# Calculate AC-optimal design and AC-efficiency for the exponential model using R package ICAOD

rm(list = ls(all = TRUE))
library(ICAOD)

# True parameter values
t1_true <- 0.9144
themu <- c(0.0863, 85, 0.1)   # (theta1, theta2, mu)
E <- themu[3]
theta <- themu[1:2]

# Estimator for the defined dose level
delta_star <- theta[2] * log(E / theta[1] + 1)

# Determine tilde_c1 in equation (17) of Remark 2
c1 <- exp(delta_star / theta[2]) - 1
c2 <- -theta[1] * delta_star * exp(delta_star / theta[2]) / theta[2]^2
tilde_c1 <- c(c1, c2)

# Gradient vector of the exponential model
fx <- function(x, param) {
  c(exp(x / param[2]) - 1,
    -param[1] * x * exp(x / param[2]) / param[2]^2)
}

# Information matrix for the exponential model under SLSE
# param = c(theta1, theta2, t1)
FIM_exp_SLSE <- function(x, w, param) {
  t1 <- param[3]
  g1 <- numeric(2)
  G2 <- matrix(0, nrow = 2, ncol = 2)
  for (i in 1:length(x)) {
    fi <- fx(x[i], param)
    g1 <- g1 + w[i] * fi
    G2 <- G2 + w[i] * (fi %*% t(fi))
  }
  list(g1 = g1, InfM = G2 - t1 * (g1 %*% t(g1)))
}

# Optimality criterion function: minimizes tilde_c1^T M^-1 tilde_c1
opt2 <- function(x, w, param, fimfunc) {
  M <- fimfunc(x = x, w = w, param = param)$InfM
  t(tilde_c1) %*% solve(M) %*% tilde_c1
}

# Sensitivity function for the criterion
opt_sens2 <- function(xi_x, x, w, param, fimfunc) {
  t1 <- param[3]
  fim <- fimfunc(x = x, w = w, param = param)
  M_inv <- solve(fim$InfM)
  f_x <- fx(xi_x, param)
  term1 <- (1 - t1) * t(tilde_c1) %*% M_inv %*% f_x %*% t(f_x) %*% M_inv %*% tilde_c1
  term2 <- t1 * t(tilde_c1) %*% M_inv %*% (f_x - fim$g1) %*% t(f_x - fim$g1) %*% M_inv %*% tilde_c1
  term3 <- t(tilde_c1) %*% M_inv %*% tilde_c1
  as.numeric(term1 + term2 - term3)
}

# =========== Find locally tilde_c-optimal design for the exponential model ============
res <- locally(
  fimfunc = FIM_exp_SLSE,
  lx = 0, ux = 150,
  inipars = c(theta, t1_true),
  iter = 500,
  k = 3,                     # number of support points
  crtfunc = opt2,
  sensfunc = opt_sens2,
  ICA.control = list(checkfreq = Inf)
)

print(res)    # display results
plot(res)     # visualise equivalence theorem

# =========== AC-optimal design computations ============
# r2 is the ratio of two variance components (assumed known)
r2 <- (0.1 * (0.25 + 0.01) + 0.9 * ((4 / 72)^2 + 6 / 72^2)) / 0.001

# Function to compute AC-optimal weights for a given design
ACopt_weight <- function(x, w, param) {
  t1 <- param[3]
  fim <- FIM_exp_SLSE(x = x, w = w, param = param)$InfM
  M_inv <- solve(fim)
  rho_val <- as.numeric(sqrt(r2 * (1 - t1) * t(tilde_c1) %*% M_inv %*% tilde_c1))
  wk1 <- 1 / (1 + rho_val)
  AC_w <- c((rho_val / (1 + rho_val)) * w, wk1)
  list(wk1 = wk1, AC_weight = AC_w)
}

# tilde_c-optimal designs for different t values (from ICAOD output)
# For t = 0
despoint_t0 <- c(83.79936, 150)
wstar_t0 <- c(0.9351299, 0.06487011)
ACopt_weight(x = despoint_t0, w = wstar_t0, param = c(theta, 0))$AC_weight

# For t = 0.9144
despoint_t0.9144 <- c(0, 95.99267, 150)
wstar_t0.9144 <- c(0.3732542, 0.5468067, 0.07993916)
ACopt_weight(x = despoint_t0.9144, w = wstar_t0.9144, param = c(theta, 0.9144))$AC_weight

# =========== AC-efficiency of the AC-optimal design ============
AC_eff <- function(x_true, w_true, param_true, x_mis, w_mis, param_mis) {
  wk1_true <- ACopt_weight(x = x_true, w = w_true, param = param_true)$wk1
  wk1_mis  <- ACopt_weight(x = x_mis,  w = w_mis,  param = param_mis)$wk1
  t1 <- param_true[3]
  
  M_inv_true <- solve(FIM_exp_SLSE(x = x_true, w = w_true, param = param_true)$InfM)
  M_inv_mis  <- solve(FIM_exp_SLSE(x = x_mis,  w = w_mis,  param = param_true)$InfM)
  
  psi_true <- ((1 - t1) * t(tilde_c1) %*% M_inv_true %*% tilde_c1) / (1 - wk1_true) + 1 / (r2 * wk1_true)
  psi_mis  <- ((1 - t1) * t(tilde_c1) %*% M_inv_mis  %*% tilde_c1) / (1 - wk1_mis)  + 1 / (r2 * wk1_mis)
  
  as.numeric(psi_true / psi_mis)
}

# Efficiency comparisons: true t = 0
AC_eff(x_true = despoint_t0, w_true = wstar_t0, param_true = c(theta, 0),
       x_mis = despoint_t0,   w_mis = wstar_t0,   param_mis = c(theta, 0))

# True t = 0.9144, misspecified as t = 0
AC_eff(x_true = despoint_t0.9144, w_true = wstar_t0.9144, param_true = c(theta, 0.9144),
       x_mis = despoint_t0,        w_mis = wstar_t0,      param_mis = c(theta, 0))
