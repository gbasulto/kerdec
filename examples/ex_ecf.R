
## ---------------------------------------------------------------- ##
## ---- Different functions to find the same ---------------------- ##
## ---------------------------------------------------------------- ##

library(kerdec)

smp <- rpois(5, lambda = 100)
t <- c(-1, 1)

ecf_re_cpp(matrix(t), matrix(smp))
ecf_real(t, smp)
Re(ecf_cpp(matrix(t), matrix(smp)))


smp <- rpois(150, lambda = 1000)
t <- seq(-1, 1, length.out = 1500)
microbenchmark::microbenchmark(
    ecf_re_cpp(matrix(t), matrix(smp)),
    ecf_real(t, smp),
    Re(ecf_cpp(matrix(t), matrix(smp)))
)

## ---------------------------------------------------------------- ##
## ---- Univariate empirical characteristic function -------------- ##
## ---------------------------------------------------------------- ##

library(kerdec)
devtools::reload()
rm(list = ls())

## Parameters of Poisson distribution, sample size and grid size
lambda <- 10
n <- 150                                # Sample size
m <- 200                                # Grid size
t <- seq(-5, 5, length.out = m)         # Evaluation grid
smp <- rpois(n, lambda)                 # Random sample

## Compute empirical characteristic values and characteristic function
## values
real <- ecf_re_cpp(matrix(t), matrix(smp))
imag <- ecf_im_cpp(matrix(t), matrix(smp))
modu <- ecf_mod_cpp(matrix(t), matrix(smp))
true <- exp(lambda*(exp(1i*t) - 1))

## Make plots
                                        # Real
plot(t, real, type = "l", col = 3)
lines(t, Re(true), col = 4)
title("Real part of empirical and true characteristic functions")
legend("topleft", legend = c("ecf", "cf"), col = 3:4, lwd = 2)

                                        # Imaginary
plot(t, imag, type = "l", col = 3)
lines(t, Im(true), col = 4)
title("Imaginary part of empirical and true characteristic functions")
legend("topleft", legend = c("ecf", "cf"), col = 3:4, lwd = 2)

                                        # Modulus
plot(t, modu, type = "l", col = 3)
lines(t, Mod(true), col = 4)
lines(t, Mod(ecf_cpp(matrix(t), matrix(smp))), col = 2)
title("Modulus of empirical and true characteristic functions")
legend("topleft", legend = c("ecf", "cf"), col = 3:4, lwd = 2)
