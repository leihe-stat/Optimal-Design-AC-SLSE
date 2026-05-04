# Plot of the Elfving space for the simple exponential model
# t = 0.2
rm(list = ls(all = TRUE))

# Gradient of eta(d, theta) with respect to theta
df <- function(x) -x * exp(-0.2 * x)

# tilde_f1 function (t passed as argument)
tilde_f1 <- function(x, t) {
  df_val <- df(x)
  df2 <- df_val^2
  lam1 <- (1 + df2 + sqrt(1 + (4 * t - 2) * df2 + df2^2)) / 2
  u1 <- (t + lam1 - 1) / (lam1 * sqrt(t))
  f1 <- sqrt(lam1 / (1 + u1^2 * df2)) * c(1, u1 * df_val)
  return(f1)
}

# tilde_f2 function (t passed as argument)
tilde_f2 <- function(x, t) {
  df_val <- df(x)
  df2 <- df_val^2
  lam2 <- (1 + df2 - sqrt(1 + (4 * t - 2) * df2 + df2^2)) / 2
  if (abs(lam2) < 1e-12) u2 <- 0 else u2 <- (t + lam2 - 1) / (lam2 * sqrt(t))
  f2 <- sqrt(lam2 / (1 + u2^2 * df2)) * c(1, u2 * df_val)
  return(f2)
}

# Parameter t
t <- 0.2

# ==================== Main routine ====================

# Draw the upper right part of the Elfving set
eps1 <- seq(-1, 1, length = 500)
eps2 <- sqrt(1 - eps1^2)
N <- length(eps1)
conv <- matrix(0, N, 2)

# Curve for x = 0 (support point 0)
for (i in 1:N) {
  conv[i, ] <- eps1[i] * tilde_f1(0, t) + eps2[i] * tilde_f2(0, t)
}
plot(conv[, 1], conv[, 2],
     xlab = "", ylab = "",
     xlim = c(-1, 1), ylim = c(-3, 3),
     col = "green", type = "l", lwd = 1,
     bty = "o", cex = 1, xaxs = "r", yaxs = "r")
mtext(expression(d[1]), side = 1, line = 2.5)
mtext(expression(d[2]), side = 2, line = 2.5)

# Support points other than 0 and 5
x0 <- c(seq(1, 4, by = 1), seq(6, 35, by = 3))
Nx <- length(x0)

# Upper right part (positive eps2)
for (j in 1:Nx) {
  for (i in 1:N) {
    conv[i, ] <- eps1[i] * tilde_f1(x0[j], t) + eps2[i] * tilde_f2(x0[j], t)
  }
  lines(conv[, 1], conv[, 2], col = "green")
}
# Curve for support point 5 (red)
for (i in 1:N) {
  conv[i, ] <- eps1[i] * tilde_f1(5, t) + eps2[i] * tilde_f2(5, t)
}
lines(conv[, 1], conv[, 2], col = "red", lwd = 1)

# Draw the lower left part of the Elfving set (negative eps2)
eps2 <- -sqrt(1 - eps1^2)
for (j in 1:Nx) {
  for (i in 1:N) {
    conv[i, ] <- eps1[i] * tilde_f1(x0[j], t) + eps2[i] * tilde_f2(x0[j], t)
  }
  lines(conv[, 1], conv[, 2], col = "green")
}
# Curve for support point 5, lower left part
for (i in 1:N) {
  conv[i, ] <- eps1[i] * tilde_f1(5, t) + eps2[i] * tilde_f2(5, t)
}
lines(conv[, 1], conv[, 2], col = "red", lwd = 1)

# Complete the prototype of the Elfving set
lines(c(0, 0), c(0, 3), col = "blue", lwd = 1)  # direction (0, tilde_c)

# Determine the intersection point using eps1*tilde_f1(5) + eps2*tilde_f2(5) = gamma * (0,1)^T
f11 <- tilde_f1(5, t)[1]
f12 <- tilde_f1(5, t)[2]
f21 <- tilde_f2(5, t)[1]
f22 <- tilde_f2(5, t)[2]

# Solve for eps1 and eps2 from: eps1*f11 + eps2*f21 = 0, eps1^2 + eps2^2 = 1
k21 <- f21 / f11
eps20 <- 1 / sqrt(1 + k21^2)
eps10 <- -k21 / sqrt(1 + k21^2)
gamma <- eps10 * f12 + eps20 * f22
points(0, gamma, col = "red")  # mark the intersection point

# Show the boundary of the generalized Elfving set
intercept <- gamma  # intercept of the tangent line (equals gamma)
# The red curve: x = eps1*f11 + eps2*f21, y = eps1*f12 + eps2*f22
# Slope of the tangent line obtained by differentiating parametric equations of the red curve
slope <- (f11 * f12 + f21 * f22) / (f21^2 + f11^2)
x_line <- seq(-1, 1, by = 0.1)
lines(x_line, slope * x_line + intercept, col = "black", lty = 5, lwd = 1.3)