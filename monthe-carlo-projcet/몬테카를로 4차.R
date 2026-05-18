############################################################
# 기침 오염 몬테카를로 (경계 발사, 앞쪽(+x)만) - 완전 실행 스크립트
# 포함 기능:
# 1) 1m x 1m, 1cm 격자(100x100)
# 2) 발사점: (x=0, y=0.5) -> 셀 중심 스냅
# 3) 공간 커널(A안): w(r,theta) = r^alpha * exp(-r/lambda) * [max(cos(theta),0)]^k
#    - +x 전방 기준: cos(theta)=dx/r, dx<=0은 0 처리(forward-only)
# 4) N ~ Poisson(mu) (기침마다 비말 수가 달라짐)
# 5) counts ~ Multinomial(N, p_ij), 오염확률맵 P_ij = 1 - exp(-mu * p_ij)
# 6) 히트맵(y축 정상), log(1+counts) 히트맵 + (색 구간별 셀 개수 표) 동시 출력
############################################################

## =========================
## 섹터 A. 유틸
## =========================

make_grid <- function(L = 1, step = 0.01) {
  dx <- step
  x <- seq(dx/2, L - dx/2, by = dx)
  y <- seq(dx/2, L - dx/2, by = dx)
  list(x = x, y = y, dx = dx, L = L)
}

snap_to_grid <- function(v, dx) {
  round((v - dx/2)/dx) * dx + dx/2
}

radial_components <- function(grid, x0, y0) {
  x <- grid$x; y <- grid$y
  dx_mat <- outer(x, y, function(xi, yj) xi - x0)
  dy_mat <- outer(x, y, function(xi, yj) yj - y0)
  r_mat  <- sqrt(dx_mat^2 + dy_mat^2)
  list(r = r_mat, dx = dx_mat, dy = dy_mat)
}

## =========================
## 섹터 B. 전방 로브(+x) & 커널
## =========================

# +x 전방 기준 전방 로브: cos(theta)=dx/r, 뒤(dx<0)는 0
forward_lobe_weight <- function(dx, r, k = 0, forward_only = TRUE) {
  eps <- 1e-12
  cos_th <- dx / pmax(r, eps)
  if (forward_only) cos_th <- pmax(cos_th, 0)
  cos_th^k
}

# A안 커널: r^alpha * exp(-r/lambda) * lobe  -> 정규화 p
kernel_prob_forward_A <- function(r, dx, lambda, k = 0, alpha = 1, forward_only = TRUE) {
  eps <- 1e-12
  radial <- pmax(r, eps)^alpha
  base   <- radial * exp(-r / lambda)
  lobe   <- forward_lobe_weight(dx, r, k = k, forward_only = forward_only)
  w <- base * lobe
  w[!is.finite(w)] <- 0
  s <- sum(w)
  if (s <= 0) stop("가중치 합이 0입니다. lambda를 키우거나 alpha/k를 줄이세요.")
  p <- w / s
  list(w = w, p = p)
}

contam_prob_map <- function(p, mu) {
  1 - exp(-mu * p)  # P(>=1) with Poisson thinning
}

simulate_once <- function(p, mu) {
  N <- rpois(1, mu)  # 매 실행마다 달라짐(재현성 고정 원하면 스크립트 맨 위에 set.seed 한 번만)
  counts_vec <- as.vector(rmultinom(1, size = N, prob = as.vector(p)))
  dim(counts_vec) <- dim(p)
  list(N = N, counts = counts_vec)
}

## =========================
## 섹터 C. lambda calibration(반경 r0 내 누적 q)
## =========================

lambda_from_fraction_box_forward_A <- function(r0, q, grid, x0, y0,
                                               k = 0, alpha = 1, forward_only = TRUE,
                                               lower = NULL, upper = NULL, max_expand = 40) {
  x0 <- snap_to_grid(x0, grid$dx)
  y0 <- snap_to_grid(y0, grid$dx)
  
  rc <- radial_components(grid, x0, y0)
  r  <- rc$r
  dx <- rc$dx
  
  if (is.null(lower)) lower <- max(1e-6, 0.2 * grid$dx)
  if (is.null(upper)) upper <- 10 * grid$L
  
  g <- function(lambda) {
    eps <- 1e-12
    rmin <- min(r)
    exp_part <- exp(-(r - rmin) / lambda)    # underflow 안정화
    radial_part <- pmax(r, eps)^alpha
    lobe <- forward_lobe_weight(dx, r, k = k, forward_only = forward_only)
    w <- radial_part * exp_part * lobe
    s <- sum(w)
    if (!is.finite(s) || s <= 0) return(NA_real_)
    p <- w / s
    sum(p[r <= r0]) - q
  }
  
  fL <- g(lower); fU <- g(upper); it <- 0
  while ((is.na(fL) || is.na(fU) || fL * fU > 0) && it < max_expand) {
    upper <- upper * 2
    fL <- g(lower); fU <- g(upper)
    it <- it + 1
  }
  if (is.na(fL)) stop("g(lower)가 NA입니다. lower를 키우거나 조건(r0,q,alpha,k)을 바꾸세요.")
  if (is.na(fU)) stop("g(upper)가 NA입니다. upper를 더 키우세요.")
  if (fL * fU > 0) stop("브래킷 내 부호 변화가 없습니다. r0, q, alpha, k를 재설정하세요.")
  
  uniroot(g, c(lower, upper))$root
}

## =========================
## 섹터 D. 플롯( y축 정상 ) + 색 구간별 셀 개수 표 동시 출력
## =========================

plot_heat <- function(M, main = "", xlab = "y (m)", ylab = "x (m)",
                      breaks = NULL, cols = NULL) {
  if (is.null(breaks)) {
    image(
      x = seq(0, 1, length.out = nrow(M)),
      y = seq(0, 1, length.out = ncol(M)),
      z = t(M),
      axes = FALSE,
      main = main,
      xlab = xlab,
      ylab = ylab
    )
  } else {
    image(
      x = seq(0, 1, length.out = nrow(M)),
      y = seq(0, 1, length.out = ncol(M)),
      z = t(M),
      axes = FALSE,
      main = main,
      xlab = xlab,
      ylab = ylab,
      breaks = breaks,
      col = cols
    )
  }
  axis(1, at = seq(0, 1, by = 0.2))
  axis(2, at = seq(0, 1, by = 0.2))
  box()
}

# log(1+counts) 히트맵 + (색 구간별 셀 개수 표) 한 번에 출력
plot_heat_with_bin_table <- function(counts, N_value, nbin = 8) {
  M <- log1p(counts)
  
  breaks <- seq(min(M), max(M), length.out = nbin + 1)
  bin_id <- cut(M, breaks = breaks, include.lowest = TRUE, right = FALSE)
  color_count_table <- table(bin_id)
  
  cols <- colorRampPalette(c("white", "orange"))(nbin)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  plot_heat(M, main = sprintf("log(1+counts), N=%d", N_value),
            breaks = breaks, cols = cols)
  
  plot.new()
  text(
    0, 1,
    paste(
      "Color-bin cell counts (log scale)\n",
      paste(names(color_count_table), as.integer(color_count_table),
            sep = " : ", collapse = "\n"),
      sep = ""
    ),
    adj = c(0, 1)
  )
  
  invisible(list(breaks = breaks, table = color_count_table))
}

## =========================
## 섹터 E. 실행(여기만 값 바꿔가며 실험)
## =========================

# (선택) 전체 재현성 고정하고 싶으면 아래 한 줄만 켜고, simulate_once 안에서 seed는 절대 쓰지 마
# set.seed(123)

grid <- make_grid(L = 1, step = 0.01)

# 발사점: +x 방향으로 분사하려면 왼쪽 경계에서 시작
x0_raw <- 0.0
y0_raw <- 0.5
x0 <- snap_to_grid(x0_raw, grid$dx)   # 0.005
y0 <- snap_to_grid(y0_raw, grid$dx)

rc <- radial_components(grid, x0, y0)
r  <- rc$r
dx <- rc$dx

# 파라미터
mu <- 500
k_lobe <- 0
alpha_radial <- 1       # A안(피크형) 강도: 1부터 시작. 피크를 더 앞으로 밀려면 alpha↑ 또는 lambda↑
r0 <- 0.5
q  <- 0.95

# lambda 역산
lambda <- lambda_from_fraction_box_forward_A(
  r0 = r0, q = q, grid = grid, x0 = x0, y0 = y0,
  k = k_lobe, alpha = alpha_radial, forward_only = TRUE
)

# 확률맵 p, 오염확률맵 P_contam
ker <- kernel_prob_forward_A(r, dx, lambda = lambda, k = k_lobe, alpha = alpha_radial, forward_only = TRUE)
p <- ker$p
P_contam <- contam_prob_map(p, mu)

# 1회 샘플
sim <- simulate_once(p, mu)
N_sampled <- sim$N
counts <- sim$counts

# 출력 1) 오염확률 히트맵
plot_heat(P_contam, main = sprintf("P[contam] (+x), λ=%.3f, μ=%d, k=%d, α=%d",
                                   lambda, mu, k_lobe, alpha_radial))

# 출력 2) log(1+counts) 히트맵 + 색 구간별 셀 개수
bin_info <- plot_heat_with_bin_table(counts, N_sampled, nbin = 8)

############################################################
# (추가) 반복 시뮬레이션으로 "시행마다 달라지는" 오염 비율 계산
# - 전체 대비 오염된 그리드 비율 (counts > 0)
# - 전방(+x) 10~40cm 구간( x in [x0+0.10, x0+0.40] ) 내 오염 비율
# - 각 시행의 N(비말 수)도 함께 저장
############################################################

T <- 200  # 기침 시행 횟수(원하는 만큼 바꾸기)

# 전방(+x) 10~40cm 구간에 해당하는 x 인덱스(한 번만 계산)
x_vals <- grid$x
x_idx  <- which(x_vals >= (x0 + 0.10) & x_vals <= (x0 + 0.40))

# 결과 저장 벡터
N_vec <- integer(T)
ratio_overall <- numeric(T)
ratio_forward <- numeric(T)
n_contam_overall <- integer(T)
n_contam_forward <- integer(T)

for (t in 1:T) {
  sim <- simulate_once(p, mu)     # 매번 N, counts가 달라짐
  N_vec[t] <- sim$N
  counts_t <- sim$counts
  
  contam_mask <- (counts_t > 0)
  
  # (1) 전체 오염 비율 및 오염 셀 수
  n_contam_overall[t] <- sum(contam_mask)
  ratio_overall[t] <- n_contam_overall[t] / length(contam_mask)
  
  # (2) 전방(+x) 10~40cm 구간 오염 비율 및 오염 셀 수
  band_mask <- contam_mask[x_idx, , drop = FALSE]
  n_contam_forward[t] <- sum(band_mask)
  ratio_forward[t] <- n_contam_forward[t] / length(band_mask)
}

# 결과 테이블
results <- data.frame(
  trial = 1:T,
  N = N_vec,
  overall_contam_cells = n_contam_overall,
  overall_contam_ratio = ratio_overall,
  forward_10_40cm_contam_cells = n_contam_forward,
  forward_10_40cm_contam_ratio = ratio_forward
)

print(head(results, 10))

cat("\n[Summary over trials]\n")
cat(sprintf("Overall contam ratio: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f\n",
            mean(ratio_overall), sd(ratio_overall),
            min(ratio_overall), max(ratio_overall)))
cat(sprintf("Forward(+x) 10~40cm ratio: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f\n",
            mean(ratio_forward), sd(ratio_forward),
            min(ratio_forward), max(ratio_forward)))

# (선택) 시행별 변동 플롯
par(mfrow=c(3,1), mar=c(4,4,2,1))
plot(results$trial, results$N, type="l", xlab="trial", ylab="N", main="Droplet count N per cough")
plot(results$trial, results$overall_contam_ratio, type="l", xlab="trial", ylab="ratio", main="Overall contaminated grid ratio (counts>0)")
plot(results$trial, results$forward_10_40cm_contam_ratio, type="l", xlab="trial", ylab="ratio", main="Forward(+x) 10~40cm contaminated grid ratio (counts>0)")
par(mfrow=c(1,1))




# 출력 3) counts 기반 요약 통계(콘솔)
summary_stats <- list(
  lambda = lambda,
  N_sampled = N_sampled,
  total_cells = length(counts),
  contaminated_cells = sum(counts > 0),
  contaminated_fraction = sum(counts > 0) / length(counts),
  max_count = max(counts),
  mean_count_all_cells = mean(counts),
  mean_count_given_contam = if (any(counts > 0)) mean(counts[counts > 0]) else 0
)
print(summary_stats)
# =========================
# 섹터 E-추가. 오염된 그리드 비율 계산(전체 + y구간 10~40cm)
# - "오염" 정의: counts > 0 (해당 셀에 비말이 1개 이상 떨어짐)
# - 전체 비율: 전체 셀 중 오염된 셀 비율
# - y축 10~40cm 구역 비율: y ∈ [0.10, 0.40]인 셀들 중 오염된 셀 비율
#   (발사점(y0) 기준 상대거리로 하고 싶으면 아래 주석 참고)
# =========================


