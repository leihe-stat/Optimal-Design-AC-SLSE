# Plot of the Elfving space for the simple exponential model
# t_1 = 0.9
rm(list = ls(all = TRUE))

# Gradient of eta(d, theta) with respect to theta
df <- function(x) -x * exp(-0.2 * x)

# tilde_f1 function, with t passed as argument
tilde_f1 <- function(x, t) {
  df_val <- df(x)
  df2 <- df_val^2
  lam1 <- (1 + df2 + sqrt(1 + (4 * t - 2) * df2 + df2^2)) / 2
  u1 <- (t + lam1 - 1) / (lam1 * sqrt(t))
  f1 <- sqrt(lam1 / (1 + u1^2 * df2)) * c(1, u1 * df_val)
  return(f1)
}

# tilde_f2 function
tilde_f2 <- function(x, t) {
  df_val <- df(x)
  df2 <- df_val^2
  lam2 <- (1 + df2 - sqrt(1 + (4 * t - 2) * df2 + df2^2)) / 2
  if (abs(lam2) < 1e-12) u2 <- 0 else u2 <- (t + lam2 - 1) / (lam2 * sqrt(t))
  f2 <- sqrt(lam2 / (1 + u2^2 * df2)) * c(1, u2 * df_val)
  return(f2)
}

# Parameter t
t <- 0.9

# ==================== Main routine ====================

# Draw the upper right part of the Elfving set
eps1 <- seq(-1, 1, length = 500)
eps2 <- sqrt(1 - eps1^2)
N <- length(eps1)
conv <- matrix(0, N, 2)

# Compute and plot the curve for x = 0 (red)
for (i in 1:N) {
  conv[i, ] <- eps1[i] * tilde_f1(0, t) + eps2[i] * tilde_f2(0, t)
}
plot(conv[, 1], conv[, 2], 
     xlab = "", ylab = "", 
     xlim = c(-1, 1), ylim = c(-3, 3),
     col = "red", type = "l", lwd = 2,
     bty = "o", cex = 1, xaxs = "r", yaxs = "r")
mtext(expression(d[1]), side = 1, line = 2.5)
mtext(expression(d[2]), side = 2, line = 2.5)

# Support points except 5
x0 <- c(seq(1, 4, by = 1), seq(6, 35, by = 3))
Nx <- length(x0)

# Upper right part (positive eps2)
for (j in 1:Nx) {
  for (i in 1:N) {
    conv[i, ] <- eps1[i] * tilde_f1(x0[j], t) + eps2[i] * tilde_f2(x0[j], t)
  }
  lines(conv[, 1], conv[, 2], col = "green")
}
# Curve for x = 5 (red, thick)
for (i in 1:N) {
  conv[i, ] <- eps1[i] * tilde_f1(5, t) + eps2[i] * tilde_f2(5, t)
}
lines(conv[, 1], conv[, 2], col = "red", lwd = 1)

# Lower left part of the Elfving set (negative eps2)
eps2 <- -sqrt(1 - eps1^2)
for (j in 1:Nx) {
  for (i in 1:N) {
    conv[i, ] <- eps1[i] * tilde_f1(x0[j], t) + eps2[i] * tilde_f2(x0[j], t)
  }
  lines(conv[, 1], conv[, 2], col = "green")
}
# Again for x = 5, lower left part
for (i in 1:N) {
  conv[i, ] <- eps1[i] * tilde_f1(5, t) + eps2[i] * tilde_f2(5, t)
}
lines(conv[, 1], conv[, 2], col = "red", lwd = 1)

# Complete the prototype of the Elfving set
lines(c(0, 0), c(0, 3), col = "blue", lwd = 1)  # direction (0, tilde_c)
points(1, 0, col = "red")                       # intersection point for support point 0

# Compute intersection point for support point 5 (x-coordinate = -0.8)
f11 <- tilde_f1(5, t)[1]
f12 <- tilde_f1(5, t)[2]
f21 <- tilde_f2(5, t)[1]
f22 <- tilde_f2(5, t)[2]

# According to equation (23): x = eps12*f11 + eps22*f21 = -0.8
# The red curve: x = eps12*f11 + eps22*f21, y = eps12*f12 + eps22*f22
# Let D = f11*f22 - f21*f12, then eps12 = (x*f22 - y*f21)/D, eps22 = (y*f11 - x*f12)/D
# Substituting eps12^2 + eps22^2 = 1 and x = -0.8 gives a quadratic in y: a*y^2 + b*y + c = 0
a1 <- 0.8 * f22
a2 <- 0.8 * f12
D <- f11 * f22 - f21 * f12
a <- f21^2 + f11^2
b <- 2 * a1 * f21 + 2 * a2 * f11
c <- a1^2 + a2^2 - D^2
y1 <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a) # y-coordinate = y1
points(-0.8, y1, col = "red")  # mark intersection for support point 5

# The black dashed boundary is defined by the two red intersection points
slope <- -y1 / (1 + 0.8)        # slope of the dashed line
x_line <- seq(-1, 1, by = 0.1)
lines(x_line, slope * (x_line - 1), col = "black", lty = 5, lwd = 1.3)
points(0, slope * (0 - 1), col = "blue")  # intersection with Elfving set boundary