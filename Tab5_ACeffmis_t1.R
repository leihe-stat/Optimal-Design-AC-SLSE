# Calculate AC-optimal design and AC-efficiency for the Michaelis-Menten model using R package ICAOD

rm(list = ls(all = TRUE))
library(ICAOD)

# True parameter values
t1_true <- 0.8
themu <- c(0.467, 25, 0.1)   # (theta1, theta2, mu)
E <- themu[3]
theta <- themu[1:2]

# Estimator for the defined dose level
delta_star <- E * theta[2] / (theta[1] - E)

# Determine tilde_c1 in equation (17) of Remark 2
c1 <- delta_star / (theta[2] + delta_star)
c2 <- -theta[1] * delta_star / (theta[2] + delta_star)^2
tilde_c1 <- c(c1, c2)

# Gradient vector of the Michaelis-Menten model
fx <- function(x, param) {
  c(x / (param[2] + x), -param[1] * x / (param[2] + x)^2)
}

# Information matrix for the M-M model under SLSE
# param = c(theta1, theta2, t1)
FIM_MM_SLSE <- function(x, w, param) {
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

# =========== Find locally tilde_c-optimal design for the M-M model ============
res <- locally(
  fimfunc = FIM_MM_SLSE,
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
r2 <- 10

# Function to compute AC-optimal weights for a given design
ACopt_weight <- function(x, w, param) {
  t1 <- param[3]
  fim <- FIM_MM_SLSE(x = x, w = w, param = param)$InfM
  M_inv <- solve(fim)
  rho_val <- as.numeric(sqrt(r2 * (1 - t1) * t(tilde_c1) %*% M_inv %*% tilde_c1))
  wk1 <- 1 / (1 + rho_val)
  AC_w <- c((rho_val / (1 + rho_val)) * w, wk1)
  list(wk1 = wk1, AC_weight = AC_w)
}

# tilde_c-optimal designs for different t values (from ICAOD output)
# For t = 0
despoint_t0 <- c(13.76211, 150)
wstar_t0 <- c(0.9167839, 0.0832161)
ACopt_weight(x = despoint_t0, w = wstar_t0, param = c(theta, 0))$AC_weight

# For t = 0.2
despoint_t0.2 <- c(14.95025, 150)
wstar_t0.2 <- c(0.9019593, 0.0980407)
ACopt_weight(x = despoint_t0.2, w = wstar_t0.2, param = c(theta, 0.2))$AC_weight

# For t = 0.3
despoint_t0.3 <- c(15.68691, 150)
wstar_t0.3 <- c(0.8929255, 0.1070745)
ACopt_weight(x = despoint_t0.3, w = wstar_t0.3, param = c(theta, 0.3))$AC_weight

# For t = 0.4
despoint_t0.4 <- c(16.56264, 150)
wstar_t0.4 <- c(0.8823778, 0.1176222)
ACopt_weight(x = despoint_t0.4, w = wstar_t0.4, param = c(theta, 0.4))$AC_weight

# For t = 0.5
despoint_t0.5 <- c(17.63514, 150)
wstar_t0.5 <- c(0.8697781, 0.1302219)
ACopt_weight(x = despoint_t0.5, w = wstar_t0.5, param = c(theta, 0.5))$AC_weight

# For t = 0.6
despoint_t0.6 <- c(0, 18.75476, 18.75851, 150)
wstar_t0.6 <- c(0.02741868, 0.3720911, 0.4614199, 0.1390704)
ACopt_weight(x = despoint_t0.6, w = wstar_t0.6, param = c(theta, 0.6))$AC_weight

# For t = 0.7
despoint_t0.7 <- c(0, 18.75, 150)
wstar_t0.7 <- c(0.16661, 0.7142857, 0.1191042)
ACopt_weight(x = despoint_t0.7, w = wstar_t0.7, param = c(theta, 0.7))$AC_weight

# For t = 0.8
despoint_t0.8 <- c(0, 18.75, 150)
wstar_t0.8 <- c(0.2707838, 0.625, 0.1042162)
ACopt_weight(x = despoint_t0.8, w = wstar_t0.8, param = c(theta, 0.8))$AC_weight

# For t = 0.9
despoint_t0.9 <- c(0, 18.75, 150)
wstar_t0.9 <- c(0.3518078, 0.5555556, 0.09263664)
ACopt_weight(x = despoint_t0.9, w = wstar_t0.9, param = c(theta, 0.9))$AC_weight

# =========== AC-efficiency of the AC-optimal design ============
AC_eff <- function(x_true, w_true, param_true, x_mis, w_mis, param_mis) {
  wk1_true <- ACopt_weight(x = x_true, w = w_true, param = param_true)$wk1
  wk1_mis  <- ACopt_weight(x = x_mis,  w = w_mis,  param = param_mis)$wk1
  t1 <- param_true[3]
  
  M_inv_true <- solve(FIM_MM_SLSE(x = x_true, w = w_true, param = param_true)$InfM)
  M_inv_mis  <- solve(FIM_MM_SLSE(x = x_mis,  w = w_mis,  param = param_true)$InfM)
  
  psi_true <- ((1 - t1) * t(tilde_c1) %*% M_inv_true %*% tilde_c1) / (1 - wk1_true) + 1 / (r2 * wk1_true)
  psi_mis  <- ((1 - t1) * t(tilde_c1) %*% M_inv_mis  %*% tilde_c1) / (1 - wk1_mis)  + 1 / (r2 * wk1_mis)
  
  as.numeric(psi_true / psi_mis)
}

# Table 5: AC-efficiency for misspecified t values
Tab5 <- matrix(0, 4, 5)

# True t = 0.3
Tab5[1, 1] <- AC_eff( # misspecified as t = 0
  x_true = despoint_t0.3, w_true = wstar_t0.3, param_true = c(theta, 0.3),
  x_mis = despoint_t0,    w_mis = wstar_t0,    param_mis = c(theta, 0))
Tab5[1, 2] <- AC_eff( # misspecified as t = 0.2
  x_true = despoint_t0.3, w_true = wstar_t0.3, param_true = c(theta, 0.3),
  x_mis = despoint_t0.2,  w_mis = wstar_t0.2,  param_mis = c(theta, 0.2))
Tab5[1, 3] <- AC_eff( # misspecified as t = 0.4
  x_true = despoint_t0.3, w_true = wstar_t0.3, param_true = c(theta, 0.3),
  x_mis = despoint_t0.4,  w_mis = wstar_t0.4,  param_mis = c(theta, 0.4))
Tab5[1, 4] <- AC_eff( # misspecified as t = 0.6
  x_true = despoint_t0.3, w_true = wstar_t0.3, param_true = c(theta, 0.3),
  x_mis = despoint_t0.6,  w_mis = wstar_t0.6,  param_mis = c(theta, 0.6))
Tab5[1, 5] <- AC_eff( # misspecified as t = 0.8
  x_true = despoint_t0.3, w_true = wstar_t0.3, param_true = c(theta, 0.3),
  x_mis = despoint_t0.8,  w_mis = wstar_t0.8,  param_mis = c(theta, 0.8))

# True t = 0.5
Tab5[2, 1] <- AC_eff( # misspecified as t = 0
  x_true = despoint_t0.5, w_true = wstar_t0.5, param_true = c(theta, 0.5),
  x_mis = despoint_t0,    w_mis = wstar_t0,    param_mis = c(theta, 0))
Tab5[2, 2] <- AC_eff( # misspecified as t = 0.2
  x_true = despoint_t0.5, w_true = wstar_t0.5, param_true = c(theta, 0.5),
  x_mis = despoint_t0.2,  w_mis = wstar_t0.2,  param_mis = c(theta, 0.2))
Tab5[2, 3] <- AC_eff( # misspecified as t = 0.4
  x_true = despoint_t0.5, w_true = wstar_t0.5, param_true = c(theta, 0.5),
  x_mis = despoint_t0.4,  w_mis = wstar_t0.4,  param_mis = c(theta, 0.4))
Tab5[2, 4] <- AC_eff( # misspecified as t = 0.6
  x_true = despoint_t0.5, w_true = wstar_t0.5, param_true = c(theta, 0.5),
  x_mis = despoint_t0.6,  w_mis = wstar_t0.6,  param_mis = c(theta, 0.6))
Tab5[2, 5] <- AC_eff( # misspecified as t = 0.8
  x_true = despoint_t0.5, w_true = wstar_t0.5, param_true = c(theta, 0.5),
  x_mis = despoint_t0.8,  w_mis = wstar_t0.8,  param_mis = c(theta, 0.8))

# True t = 0.7
Tab5[3, 1] <- AC_eff( # misspecified as t = 0
  x_true = despoint_t0.7, w_true = wstar_t0.7, param_true = c(theta, 0.7),
  x_mis = despoint_t0,    w_mis = wstar_t0,    param_mis = c(theta, 0))
Tab5[3, 2] <- AC_eff( # misspecified as t = 0.2
  x_true = despoint_t0.7, w_true = wstar_t0.7, param_true = c(theta, 0.7),
  x_mis = despoint_t0.2,  w_mis = wstar_t0.2,  param_mis = c(theta, 0.2))
Tab5[3, 3] <- AC_eff( # misspecified as t = 0.4
  x_true = despoint_t0.7, w_true = wstar_t0.7, param_true = c(theta, 0.7),
  x_mis = despoint_t0.4,  w_mis = wstar_t0.4,  param_mis = c(theta, 0.4))
Tab5[3, 4] <- AC_eff( # misspecified as t = 0.6
  x_true = despoint_t0.7, w_true = wstar_t0.7, param_true = c(theta, 0.7),
  x_mis = despoint_t0.6,  w_mis = wstar_t0.6,  param_mis = c(theta, 0.6))
Tab5[3, 5] <- AC_eff( # misspecified as t = 0.8
  x_true = despoint_t0.7, w_true = wstar_t0.7, param_true = c(theta, 0.7),
  x_mis = despoint_t0.8,  w_mis = wstar_t0.8,  param_mis = c(theta, 0.8))

# True t = 0.9
Tab5[4, 1] <- AC_eff( # misspecified as t = 0
  x_true = despoint_t0.9, w_true = wstar_t0.9, param_true = c(theta, 0.9),
  x_mis = despoint_t0,    w_mis = wstar_t0,    param_mis = c(theta, 0))
Tab5[4, 2] <- AC_eff( # misspecified as t = 0.2
  x_true = despoint_t0.9, w_true = wstar_t0.9, param_true = c(theta, 0.9),
  x_mis = despoint_t0.2,  w_mis = wstar_t0.2,  param_mis = c(theta, 0.2))
Tab5[4, 3] <- AC_eff( # misspecified as t = 0.4
  x_true = despoint_t0.9, w_true = wstar_t0.9, param_true = c(theta, 0.9),
  x_mis = despoint_t0.4,  w_mis = wstar_t0.4,  param_mis = c(theta, 0.4))
Tab5[4, 4] <- AC_eff( # misspecified as t = 0.6
  x_true = despoint_t0.9, w_true = wstar_t0.9, param_true = c(theta, 0.9),
  x_mis = despoint_t0.6,  w_mis = wstar_t0.6,  param_mis = c(theta, 0.6))
Tab5[4, 5] <- AC_eff( # misspecified as t = 0.8
  x_true = despoint_t0.9, w_true = wstar_t0.9, param_true = c(theta, 0.9),
  x_mis = despoint_t0.8,  w_mis = wstar_t0.8,  param_mis = c(theta, 0.8))

# Display the resulting efficiency table
Tab5