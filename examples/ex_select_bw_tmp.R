## Case 1: normal error with known variance (misspecified
## distribution)
case1 <-
  kerdec_dens(Y, method = "CV", kernel = "sinc", h0 = c(0.08, 0.5),
              lower = lower, upper = upper, h = NULL,
              error_dist = "laplace",
              error_scale_par = sd_error, resolution = 64)

with(case1, plot(x, f_vals, type = "l"))

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "vp", h0 = c(0.08, 0.5),
              lower = lower, upper = upper, h = NULL,
              error_dist = "normal",
              error_scale_par = sd_error, resolution = 128)

with(case1, plot(x, f_vals, type = "l"))

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "triw", h0 = c(0.08, 0.5),
              lower = lower, upper = upper, h = NULL,
              error_dist = "normal",
              error_scale_par = sd_error, resolution = 128)

with(case1, plot(x, f_vals, type = "l"))

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "tric", h0 = c(0.08, 0.5),
              lower = lower, upper = upper, h = NULL,
              error_dist = "normal",
              error_scale_par = sd_error, resolution = 128)

with(case1, plot(x, f_vals, type = "l"))

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "flat", h0 = c(0.08, 0.5),
              lower = lower, upper = upper, h = NULL,
              error_dist = "normal",
              error_scale_par = sd_error, resolution = 128)

with(case1, plot(x, f_vals, type = "l"))
