library(kerdec)

## Settings and samples
n <- 150                                # Sample size
l <- 5                                  # Number of columns
m <- n + 10                             # Error sample size
shape <- 5                              # X distr. shape par.
rate <- 5                               # X distr. rate par.
sd_error <- 0.2                          # std. error of error distr.
X <- rgamma(n, shape, rate)             # Uncontaminated sample
eps_panel <- matrix(rlaplace(n*l, sd = sd_error),
                    n, l)               # Panel of errors
eps <- rlaplace(m, sd = sd_error)       # Pure errors
del <- rlaplace(m, sd = sd_error) - rlaplace(m, sd = sd_error)
Y <- X + eps_panel[, 1]                 # Contaminated sample
Y_panel <- sweep(x = eps_panel, MARGIN = 1,
                 STATS = X, FUN = "+")  # Contaminated in panel



t <- seq(0, 50, .2)
## out0 <- dens_denominator(t = t, smp = eps, sigma = sd(eps),
##                          k = 1, error_dist = 2, panel_proc = 1)
out1 <- dens_denominator(t = t, smp = eps, sigma = sd_error,
                         k = 1, error_dist = 1, panel_proc = 1)
out2 <- dens_denominator(t = t, smp = eps, sigma = sd_error,
                         k = 1, error_dist = 2, panel_proc = 1)
out3 <- dens_denominator(t = t, smp = del, sigma = sd_error,
                         k = 2, error_dist = 1, panel_proc = 1)
out4 <- dens_denominator(t = t, smp = eps, sigma = sd_error,
                         k = 1, error_dist = 3, panel_proc = 1)
plot(range(t), range(c(out1, out2)), type = "n")
## lines(t, out0, col = 1, lty = 1, lwd = 3)
lines(t, out1, col = 2, lty = 1, lwd = 1.5)
lines(t, out2, col = 3, lty = 2, lwd = 1.5)
lines(t, out3, col = 4, lty = 3, lwd = 1.5)
lines(t, out4, col = "magenta", lty = 4, lwd = 1.5)
legend("topright",
       legend = c("ecf", "param", "ecf_diffs", "param_normal"),
       col = c("red", "green", "blue", "magenta"),
       lty = 1:4)


microbenchmark::microbenchmark(
    Mod(sapply(t, function(tt) mvdeconvolution::ecf(matrix(tt), matrix(eps)))),
dens_denominator(t = t, smp = eps, sigma = sd_error,
                         k = 1, error_dist = 1, panel_proc = 1),
dens_denominator(t = t, smp = eps, sigma = sd_error,
                         k = 1, error_dist = 2, panel_proc = 1),
dens_denominator(t = t, smp = del, sigma = sd_error,
                         k = 2, error_dist = 1, panel_proc = 1),
dens_denominator(t = t, smp = eps, sigma = sd_error,
            k = 1, error_dist = 3, panel_proc = 1)
)



t <- seq(0, 50, .1)
outt1 <- dens_denominator(t = t, smp = del, sigma = sd_error,
                         k = 20, error_dist = 1, panel_proc = 2)
outt2 <- dens_denominator(t = t, smp = del, sigma = sd_error,
                         k = 20, error_dist = 2, panel_proc = 2)
outt3 <- dens_denominator(t = t, smp = del, sigma = sd_error,
                         k = 20, error_dist = 3, panel_proc = 2)
plot(range(t), range(c(out1, out2)), type = "n")
## lines(t, out0, col = 1, lty = 1, lwd = 3)
lines(t, outt1, col = 2, lty = 1, lwd = 1.5)
lines(t, outt2, col = 3, lty = 2, lwd = 1.5)
lines(t, outt3, col = 4, lty = 3, lwd = 1.5)
legend("topright",
       legend = c("ecf", "param", "param_normal"),
       col = c("red", "green", "blue"),
       lty = 1:3)
