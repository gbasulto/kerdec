
library(kerdec)

## Seed
set.seed(666)

## Settings and samples
n <- 150                                # Sample size
l <- 5                                  # Number of columns
m <- n + 10                             # Error sample size
shape <- 5                              # X distr. shape par.
rate <- 5                               # X distr. rate par.
sd_error <- .2                          # std. error of error distr.
X <- rgamma(n, shape, rate)             # Uncontaminated sample
eps_panel <- matrix(rlaplace(n*l, sd = sd_error),
                    n, l)               # Panel of errors
eps <- rlaplace(m, sd = sd_error)       # Pure errors
Y <- X + eps_panel[, 1]                 # Contaminated sample
Y_panel <- sweep(x = eps_panel, MARGIN = 1,
                 STATS = X, FUN = "+")  # Contaminated in panel


## Real density
plot(function(x) dgamma(x, shape, rate), 0, 4)
#hist(Y)

## Parameters to perform estimation
lower <- -2
upper <- 7

## Case 1: normal error with known variance (misspecified
## distribution)
case1 <-
    kerdec_dens(Y, method = "NR", kernel = "triw", h0 = c(.06, 0.10),
                lower = lower, upper = upper, h = 0.08,
                error_dist = "normal",
                error_scale_par = sd_error)
with(case1, plot(x, f_vals, type = "l"))

## Case 2: Laplace error with known variance
case2 <-
    kerdec_dens(Y, method = "NR", kernel = "triw", h0 = c(.05, 0.1),
                lower = lower, upper = upper, h = 0.068,
                error_dist = "Laplace",
                error_scale_par = sd_error)
with(case2, plot(x, f_vals, type = "l"))

## Case 3: normal error with unknown variance and sample of errors
## given (misspecified distribution)
case3 <-
    kerdec_dens(Y, method = "CV", kernel = "flat",
                lower = lower, upper = upper, h = 0.2,
                error_dist = "Normal", error_smp = eps)
with(case3, plot(x, f_vals, type = "l"))

## Case 4: Laplace error with known variance and sample of errors
## given
case4 <-
    kerdec_dens(Y, method = "CV", kernel = "flat",
                lower = lower, upper = upper, h = 0.2,
                error_dist = "Laplace", error_smp = eps)
with(case4, plot(x, f_vals, type = "l"))

## Case 4: ecf to approximate errors based on sample of errors
## given
case5 <-
    kerdec_dens(Y, method = "CV", kernel = "flat",
                lower = lower, upper = upper, h = 0.2,
                error_smp = eps)
with(case5, plot(x, f_vals, type = "l"))

## Case 6: Panel data with normal errors (unknown variances)
case6 <- 
  kerdec_dens(Y_panel, method = "CV", kernel = "flat",
              lower = lower, upper = upper, h = 0.15, 
              error_dist = "normal")
with(case6, plot(x, f_vals, type = "l"))

## Case 7: Panel data with Laplace errors (unknown variances)
case7 <- 
  kerdec_dens(smp = Y_panel, method = "CV", kernel = "flat",
              lower = lower, upper = upper, h = 0.15, 
              error_dist = "laplace")
with(case7, plot(x, f_vals, type = "l"))

## Case 8: Panel data with Laplace errors (unknown variances)
case8 <- 
  kerdec_dens(smp = Y_panel, method = "CV", kernel = "flat",
              lower = lower, upper = upper, h = 0.15, 
              error_dist = "none")
with(case8, plot(x, f_vals, type = "l"))

## Case 9: Panel data with normal errors (unknown variances)
case9 <- 
  kerdec_dens(Y_panel, method = "CV", kernel = "flat",
              lower = lower, upper = upper, h = 0.15, 
              error_dist = "normal", panel_proc = "take_aver")
with(case9, plot(x, f_vals, type = "l"))

## Case 10: Panel data with Laplace errors (unknown variances)
case10 <- 
  kerdec_dens(smp = Y_panel, method = "CV", kernel = "flat",
              lower = lower, upper = upper, h = 0.15, 
              error_dist = "laplace")
with(case10, plot(x, f_vals, type = "l"), panel_proc = "take_aver")

## Case 11: Panel data with Laplace errors (unknown variances)
case11 <- 
  kerdec_dens(smp = Y_panel, method = "CV", kernel = "flat",
              lower = lower, upper = upper, h = 0.15, 
              error_dist = "none")
with(case11, plot(x, f_vals, type = "l"), panel_proc = "take_aver")



case0 <-
  kerdec_dens(Y, method = "NR", kernel = "triw",
              lower = lower, upper = upper, h = 0.07,
              h0 = c(0.03, 0.2),
              ## h0 = c(.04, 0.15),
              error_dist = "laplace",
              error_scale_par = sd_error)

with(case0, plot(x, f_vals, type = "l"))

kerdec_dens(Y, method = "NR", kernel = "triw",
            lower = lower, upper = upper, h = 0.07,
            h0 = c(0.05, 0.2),
            ## h0 = c(.04, 0.15),
            error_dist = "normal",
            error_scale_par = sd_error)
