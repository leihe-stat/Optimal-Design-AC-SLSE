#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All functions for evaluating finite-sample properties
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mixed error distribution (contaminated normal + scaled chi-square)
error.dist <- function(N) {
  a <- rbinom(N, 1, 0.1)
  rand.x <- numeric(N)
  for (i in 1:N) {
    if (a[i] == 1) {
      rand.x[i] <- rnorm(1, 0.5, 0.1)
    } else {
      rand.x[i] <- (rchisq(1, 3) - 7) / 72
    }
  }
  return(rand.x)
}

# Random sampling from a design
generate_random_sample <- function(xi, n1, params, mi) {
  x <- rep(xi[1, ], round(xi[2, ] * n1))
  error <- error.dist(length(x))
  y <- mi(x, params) + error
  return(list(x = x, y = y, error = error))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Optimal weighting matrix (with step control)
W <- function(x, params, mi, step = 2) {
  if (step == 1) {
    # Return 2x2 identity matrix
    return(diag(2))
  } else if (step == 2) {
    mval <- mi(x, params)
    # The constants Tmu2.est, Tmu3.est, Tmu4.est are defined globally
    constant <- 1 / (Tmu2.est * (Tmu4.est - Tmu2.est^2) - Tmu3.est^2)
    To.matrix <- matrix(c(Tmu4.est + 4 * Tmu3.est * mval + 4 * Tmu2.est * mval^2 - Tmu2.est^2,
                          -Tmu3.est - 2 * Tmu2.est * mval,
                          -Tmu3.est - 2 * Tmu2.est * mval,
                          Tmu2.est), ncol = 2)
    return(constant * To.matrix)
  } else {
    stop("step must be 1 or 2")
  }
}

# Objective function Q_n for the SLSE criterion (with step control)
Qn.m1 <- function(params, step = 2) {
  objvec <- cbind(y1.obs - m1(x1.obs, params),
                  y1.obs^2 - m1(x1.obs, params)^2 - params[3])
  objsum <- 0
  for (i in 1:length(y1.obs)) {
    objsum <- t(objvec[i, ]) %*% W(x1.obs[i], m1TRUE.est, m1, step = step) %*% objvec[i, ] + objsum
  }
  return(as.numeric(objsum))
}

# Ordinary least squares criterion (no variance parameter)
Qn.m1.OLS <- function(params) {
  objvec <- y1.obs - m1(x1.obs, params)
  obj_sum <- 0
  for (i in 1:length(y1.obs)) {
    obj_sum <- objvec[i]^2 + obj_sum
  }
  return(obj_sum)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Michaelis-Menten regression function (two parameters)
m1 <- function(x, theta) theta[1] * x / (theta[2] + x)

# Partial derivatives of m1 w.r.t. theta
Deriv.m1 <- function(x, theta) c(x / (theta[2] + x), -theta[1] * x / (theta[2] + x)^2)

# Derivative of rho function w.r.t. gamma = (theta, sigma^2)
Deriv.rho.m1 <- function(x, theta) {
  First3row <- cbind(Deriv.m1(x, theta), 2 * m1(x, theta) * Deriv.m1(x, theta))
  return(-rbind(First3row, c(0, 1)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Second derivatives of m1 w.r.t. theta (Hessian)
Deriv2.m1 <- function(x, theta) {
  matrix(c(0, -x / (theta[2] + x)^2,
           -x / (theta[2] + x)^2,
           2 * theta[1] * x / (theta[2] + x)^3), ncol = 2)
}

# Derivative of vec(rho) w.r.t. gamma
Deriv.vecrho.m1 <- function(x, theta, dim = 2) {
  D1 <- Deriv.m1(x, theta)
  D2 <- Deriv2.m1(x, theta)
  First3col <- rbind(D2, rep(0, dim), 2 * m1(x, theta) * D2 + 2 * D1 %*% t(D1), rep(0, dim))
  return(-cbind(First3col, rep(0, 2 * dim + 2)))
}

# Single observation gradient and Jacobian for Q_n(gamma)
gradient_Jacobian_single_m1 <- function(gamma1, y, x, step = 2) {
  rho <- c(y - m1(x, gamma1), y^2 - m1(x, gamma1)^2 - gamma1[3])
  gradient_single <- Deriv.rho.m1(x, gamma1) %*% W(x, m1TRUE.est, m1, step = step) %*% rho
  Jacobian_single <- Deriv.rho.m1(x, gamma1) %*% W(x, m1TRUE.est, m1, step = step) %*% t(Deriv.rho.m1(x, gamma1)) +
    kronecker(t(rho) %*% W(x, m1TRUE.est, m1, step = step), diag(rep(1, 3))) %*% Deriv.vecrho.m1(x, gamma1)
  return(list(gradient_single = gradient_single, Jacobian_single = Jacobian_single))
}

# Objective function for Newton-Raphson optimization (with step control)
opt.m1.R <- function(gamma1, step = 2) {
  objvec <- cbind(y1.obs - m1(x1.obs, gamma1), y1.obs^2 - m1(x1.obs, gamma1)^2 - gamma1[3])
  obj_sum <- 0
  gradient_sum <- rep(0, 3)
  Jacobian_sum <- matrix(0, 3, 3)
  for (i in 1:length(y1.obs)) {
    obj_sum <- t(objvec[i, ]) %*% W(x1.obs[i], m1TRUE.est, m1, step = step) %*% objvec[i, ] + obj_sum
    tmp <- gradient_Jacobian_single_m1(gamma1, y1.obs[i], x1.obs[i], step = step)
    gradient_sum <- tmp$gradient_single + gradient_sum
    Jacobian_sum <- tmp$Jacobian_single + Jacobian_sum
  }
  structure(obj_sum, gradient = 2 * gradient_sum, hessian = 2 * Jacobian_sum)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation settings
totalsize <- 60

### AC-optimal design based on the SLSE (second-order least squares)
xi1.opt<-matrix(c(0, 18.75, 150, 0.3620156, 0.5468066, 0.09117779),byrow=T, ncol = 3)
sampsize<-round(totalsize*(1-0.30390082))

### AC-optimal design based on the OLSE (ordinary least squares)
# xi1.opt <- matrix(c(13.76211, 150, 0.9167839, 0.08321611), byrow = TRUE, ncol = 2)
# sampsize <- round(totalsize * (1 - 0.17855257))

# True parameters and error moments
m1TRUE <- c(0.467, 25)
Tmu2 <- 0.02981944
Tmu3 <- 0.01370836
Tmu4 <- 0.007769021

# Number of Monte Carlo replications
N <- 1000

# Storage matrices
gamma.m1.SLS <- matrix(0, N, 3)
gamma.m1.OLS <- matrix(0, N, 2)

# Simulation loop
for (i in 1:N) {
  # Generate data
  data.m1 <- generate_random_sample(xi1.opt, sampsize, m1TRUE, m1)
  x1.obs <- data.m1$x
  y1.obs <- data.m1$y
  error.obs <- data.m1$error
  # asyest.m1 <- Rsolnp::solnp(pars = c(0.2, 1, 3),
  #                            fun = function(p) Qn.m1(p, step = 1),
  #                            LB = c(1e-6, 1e-6, 1e-6),
  #                            UB = c(10, 100, 5),
  #                            control = list(trace = 0))$par
  # m1TRUE.est <- asyest.m1[1:2]; Tmu2.est <- asyest.m1[3]; 
  # Tmu3.est <- mean(error.obs^3); Tmu4.est <-  mean(error.obs^4)
  # Use true parameters for weighting matrix (to avoid re-estimation)
  m1TRUE.est <- m1TRUE
  Tmu2.est <- Tmu2
  Tmu3.est <- Tmu3
  Tmu4.est <- Tmu4
  
  # Initial estimate for SLS via solnp (step = 2)
  asyest.m1 <- Rsolnp::solnp(pars = c(0.2, 1, 3),
                             fun = function(p) Qn.m1(p, step = 2),
                             LB = c(1e-6, 1e-6, 1e-6),
                             UB = c(Inf, Inf, Inf),
                             control = list(trace = 0))$par
  
  # Refine using Newton-Raphson
  thetamu2.m1 <- lava::NR(asyest.m1, opt.m1.R, step = 2, control = list(tol = 1e-6))$par
  gamma.m1.SLS[i, ] <- thetamu2.m1
  
  # OLS estimation (only theta)
  asyest.m1.OLS <- Rsolnp::solnp(pars = c(0.2, 1),
                                 fun = Qn.m1.OLS,
                                 LB = c(1e-6, 1e-6),
                                 UB = c(Inf, Inf),
                                 control = list(trace = 0))$par
  gamma.m1.OLS[i, ] <- asyest.m1.OLS
  
  if (i %% 100 == 0) cat("i =", i, "\n")
  flush.console()
}

# Optional: compute MSE
# apply((t(gamma.m1.SLS) - c(m1TRUE, Tmu2))^2, 1, mean)
# apply((t(gamma.m1.OLS) - m1TRUE)^2, 1, mean)

# Monte Carlo estimation of mu in active control
To.mu <- numeric(N)
for (i in 1:N) {
  mu <- 0.1 + (rchisq(totalsize - sampsize, 3) - 3) / sqrt(6000)
  To.mu[i] <- mean(mu)
}

# Estimation of target dose d* = mu * theta2 / (theta1 - mu)
To.dstar <- matrix(0, N, 2)
for (i in 1:N) {
  To.dstar[i, 1] <- To.mu[i] * gamma.m1.SLS[i, 2] / (gamma.m1.SLS[i, 1] - To.mu[i])
  To.dstar[i, 2] <- To.mu[i] * gamma.m1.OLS[i, 2] / (gamma.m1.OLS[i, 1] - To.mu[i])
}

# Standard deviations of estimated target doses
apply(To.dstar, 2, sd)