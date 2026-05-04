# Calculate AC-optimal design and AC-efficiency for the Michaelis-Menten model using R package ICAOD

rm(list = ls(all = TRUE))
library(ICAOD)

# True parameter values
t1_true <- 0.9

# Parameter sets (theta1, theta2, mu) for different misspecification scenarios
themu1 <- c(0.467 * 1.2, 25 * 1.3, 0.1 * 1.4)
themu2 <- c(0.467 * 1.2, 25 * 1.0, 0.1 * 1.0)
themu3 <- c(0.467 * 1.2, 25 * 0.7, 0.1 * 0.6)
themu4 <- c(0.467 * 1.0, 25 * 1.3, 0.1 * 1.0)
themu5 <- c(0.467 * 1.0, 25 * 1.0, 0.1 * 0.6)
themu6 <- c(0.467 * 1.0, 25 * 0.7, 0.1 * 1.4)
themu7 <- c(0.467 * 0.8, 25 * 1.3, 0.1 * 0.6)
themu8 <- c(0.467 * 0.8, 25 * 1.0, 0.1 * 1.4)
themu9 <- c(0.467 * 0.8, 25 * 0.7, 0.1 * 1.0)
themu_true <- c(0.467 * 1.0, 25 * 1.0, 0.1 * 1.0)  # true parameter values

# Function to compute tilde_c1 as a function of (theta1, theta2, mu)
ftilde_c1 <- function(themu) {
  E <- themu[3]
  theta <- themu[1:2]
  # Estimator for the defined dose level
  delta_star <- E * theta[2] / (theta[1] - E)
  # Determine tilde_c1 in equation (17) of Remark 2
  c1 <- delta_star / (theta[2] + delta_star)
  c2 <- -theta[1] * delta_star / (theta[2] + delta_star)^2
  c(c1, c2)
}

# Gradient vector of the Michaelis-Menten model
fx <- function(x, param) {
  c(x / (param[2] + x), -param[1] * x / (param[2] + x)^2)
}

# Information matrix for the M-M model under SLSE (t1 is fixed globally)
FIM_MM_SLSE <- function(x, w, param) {
  # param contains (theta1, theta2); mu is not used here
  g1 <- numeric(2)
  G2 <- matrix(0, nrow = 2, ncol = 2)
  for (i in 1:length(x)) {
    fi <- fx(x[i], param)
    g1 <- g1 + w[i] * fi
    G2 <- G2 + w[i] * (fi %*% t(fi))
  }
  list(g1 = g1, InfM = G2 - t1_true * (g1 %*% t(g1)))
}

# Optimality criterion function: minimizes tilde_c1^T M^-1 tilde_c1
opt2 <- function(x, w, param, fimfunc) {
  M <- fimfunc(x = x, w = w, param = param)$InfM
  t(tilde_c1) %*% solve(M) %*% tilde_c1
}

# Sensitivity function for the criterion
opt_sens2 <- function(xi_x, x, w, param, fimfunc) {
  fim <- fimfunc(x = x, w = w, param = param)
  M_inv <- solve(fim$InfM)
  f_x <- fx(xi_x, param)
  term1 <- (1 - t1_true) * t(tilde_c1) %*% M_inv %*% f_x %*% t(f_x) %*% M_inv %*% tilde_c1
  term2 <- t1_true * t(tilde_c1) %*% M_inv %*% (f_x - fim$g1) %*% t(f_x - fim$g1) %*% M_inv %*% tilde_c1
  term3 <- t(tilde_c1) %*% M_inv %*% tilde_c1
  as.numeric(term1 + term2 - term3)
}

# =========== Find locally tilde_c-optimal design for the M-M model ============
themu <- themu_true
tilde_c1 <- ftilde_c1(themu)
res <- locally(
  fimfunc = FIM_MM_SLSE,
  lx = 0, ux = 150,
  inipars = themu[1:2],          # only (theta1, theta2) are optimized
  iter = 500,
  k = 3,                         # number of support points
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
  tilde_c1 <- ftilde_c1(param)
  fim <- FIM_MM_SLSE(x = x, w = w, param = param[1:2])$InfM
  M_inv <- solve(fim)
  rho_val <- as.numeric(sqrt(r2 * (1 - t1_true) * t(tilde_c1) %*% M_inv %*% tilde_c1))
  wk1 <- 1 / (1 + rho_val)
  AC_w <- c((rho_val / (1 + rho_val)) * w, wk1)
  list(wk1 = wk1, AC_weight = AC_w)
}

# tilde_c-optimal designs for different parameter values (from ICAOD output)
# For true parameters (t = 0.9)
despoint_true_0.9 <- c(0, 18.75, 150)
wstar_true_0.9 <- c(0.3518078, 0.5555556, 0.09263664)
ACopt_weight(x = despoint_true_0.9, w = wstar_true_0.9, param = themu_true)$AC_weight

# For Case 1
despoint_1_0.9 <- c(0, 22.67442, 150)
wstar_1_0.9 <- c(0.3662052, 0.5555556, 0.07823924)
ACopt_weight(x = despoint_1_0.9, w = wstar_1_0.9, param = themu1)$AC_weight

# For Case 2
despoint_2_0.9 <- c(0, 18.75, 150)
wstar_2_0.9 <- c(0.3420723, 0.5555556, 0.1023721)
ACopt_weight(x = despoint_2_0.9, w = wstar_2_0.9, param = themu2)$AC_weight

# For Case 3
despoint_3_0.9 <- c(0, 14.18919, 150)
wstar_3_0.9 <- c(0.3244156, 0.5555556, 0.1200288)
ACopt_weight(x = despoint_3_0.9, w = wstar_3_0.9, param = themu3)$AC_weight

# For Case 4
despoint_4_0.9 <- c(0, 22.67442, 150)
wstar_4_0.9 <- c(0.3544884, 0.5555556, 0.08995603)
ACopt_weight(x = despoint_4_0.9, w = wstar_4_0.9, param = themu4)$AC_weight

# For Case 5
despoint_5_0.9 <- c(0, 18.75, 150)
wstar_5_0.9 <- c(0.3300448, 0.5555556, 0.1143997)
ACopt_weight(x = despoint_5_0.9, w = wstar_5_0.9, param = themu5)$AC_weight

# For Case 6
despoint_6_0.9 <- c(0, 14.18919, 150)
wstar_6_0.9 <- c(0.3754471, 0.5555556, 0.06899737)
ACopt_weight(x = despoint_6_0.9, w = wstar_6_0.9, param = themu6)$AC_weight

# For Case 7
despoint_7_0.9 <- c(0, 22.67442, 150)
wstar_7_0.9 <- c(0.3392844, 0.5555556, 0.1051601)
ACopt_weight(x = despoint_7_0.9, w = wstar_7_0.9, param = themu7)$AC_weight

# For Case 8
despoint_8_0.9 <- c(0, 18.75, 150)
wstar_8_0.9 <- c(0.4134432, 0.5555556, 0.0310012)
ACopt_weight(x = despoint_8_0.9, w = wstar_8_0.9, param = themu8)$AC_weight

# For Case 9
despoint_9_0.9 <- c(0, 14.18919, 150)
wstar_9_0.9 <- c(0.3647663, 0.5555556, 0.07967818)
ACopt_weight(x = despoint_9_0.9, w = wstar_9_0.9, param = themu9)$AC_weight

# =========== AC-efficiency of the AC-optimal design ============
AC_eff <- function(x_true, w_true, param_true, x_mis, w_mis, param_mis) {
  wk1_true <- ACopt_weight(x = x_true, w = w_true, param = param_true)$wk1
  wk1_mis  <- ACopt_weight(x = x_mis,  w = w_mis,  param = param_mis)$wk1
  
  tilde_c1 <- ftilde_c1(param_true)
  
  M_inv_true <- solve(FIM_MM_SLSE(x = x_true, w = w_true, param = param_true[1:2])$InfM)
  M_inv_mis  <- solve(FIM_MM_SLSE(x = x_mis,  w = w_mis,  param = param_true[1:2])$InfM)
  
  psi_true <- ((1 - t1_true) * t(tilde_c1) %*% M_inv_true %*% tilde_c1) / (1 - wk1_true) + 1 / (r2 * wk1_true)
  psi_mis  <- ((1 - t1_true) * t(tilde_c1) %*% M_inv_mis  %*% tilde_c1) / (1 - wk1_mis)  + 1 / (r2 * wk1_mis)
  
  as.numeric(psi_true / psi_mis)
}

# Table 6: AC-efficiency when true model is themu_true but design is based on misspecified scenarios
Tab6 <- matrix(0, 10, 1)

# True t = 0.9, misspecified as Case 1
Tab6[1, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_1_0.9,    w_mis = wstar_1_0.9,    param_mis = themu1)
# Misspecified as Case 2
Tab6[2, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_2_0.9,    w_mis = wstar_2_0.9,    param_mis = themu2)
# Misspecified as Case 3
Tab6[3, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_3_0.9,    w_mis = wstar_3_0.9,    param_mis = themu3)
# Misspecified as Case 4
Tab6[4, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_4_0.9,    w_mis = wstar_4_0.9,    param_mis = themu4)
# Misspecified as Case 5
Tab6[5, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_5_0.9,    w_mis = wstar_5_0.9,    param_mis = themu5)
# Misspecified as Case 6
Tab6[6, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_6_0.9,    w_mis = wstar_6_0.9,    param_mis = themu6)
# Misspecified as Case 7
Tab6[7, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_7_0.9,    w_mis = wstar_7_0.9,    param_mis = themu7)
# Misspecified as Case 8
Tab6[8, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_8_0.9,    w_mis = wstar_8_0.9,    param_mis = themu8)
# Misspecified as Case 9
Tab6[9, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                     x_mis = despoint_9_0.9,    w_mis = wstar_9_0.9,    param_mis = themu9)
# No misspecification (true model)
Tab6[10, 1] <- AC_eff(x_true = despoint_true_0.9, w_true = wstar_true_0.9, param_true = themu_true,
                      x_mis = despoint_true_0.9,  w_mis = wstar_true_0.9,  param_mis = themu_true)

# Display the efficiency table
Tab6