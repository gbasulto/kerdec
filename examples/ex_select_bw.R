
library(kerdec)

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


plot(function(x) dgamma(x, shape, rate), 0, 4)

select_bw(Y,
          method = c("CV", "NR")[1],
          kernel = c("sinc", "vp", "triw", "tric",
                     "flat")[5],
          h0 = NULL,
          error_dist = c("Normal", "Laplace", "None")[1],
          error_scale_par = NULL,
          error_smp = rnorm(1000, 0, 3^2),
          resolution = 128,
          error_proc = c("all", "vs_first",
                                     "indep_pairs")[1],
          panel_proc = c("keep_first", "take_aver")[1],
          truncation_bound = NULL)

select_bw(Y,
          method = c("CV", "NR")[1],
          kernel = c("sinc", "vp", "triw", "tric",
                     "flat")[5],
          h0 = NULL,
          error_dist = c("Normal", "Laplace", "None")[2],
          error_scale_par = NULL,
          error_smp = rlaplace(1000, 0, 3^2),
          resolution = 128,
          error_proc = c("all", "vs_first",
                                     "indep_pairs")[1],
          panel_proc = c("keep_first", "take_aver")[1],
          truncation_bound = NULL)



select_bw(
    matrix(1:12, 3, 3),
          method = "nR",
          kernel = "flat",
          h0 = NULL,
          error_dist = "normal",
          error_scale_par = NULL,
          error_smp = NULL,
          resolution = 128,
          error_proc = "all",
          panel_proc = c("keep_first", "take_aver")[1],
          truncation_bound = NULL)


nn <- 100
sig <- 50
smp <- rlaplace(nn, 0, sig) - rlaplace(nn, 0, sig)

sig0 <- sd(smp)/sqrt(2); sig0

microbenchmark::microbenchmark(
    optim(par = sig0, fn = laplace_convol_loglik,
          lower = .Machine$double.eps, method = "L-BFGS-B",
          smp = smp),    
    optim(par = sig0,
          fn = laplace_convol_loglik,
          gr = dlaplace_convol_loglik,
          lower = .Machine$double.eps, method = "L-BFGS-B",
          smp = smp)
    )

plot(function(t) laplace_convol_loglik(t, smp), 0.35, 0.7)
plot(function(t) dlaplace_convol_loglik(t, smp), 0.35, 0.7)

numDeriv::grad(laplace_convol_loglik, .5, smp = smp)
dlaplace_convol_loglik(.5, smp)
