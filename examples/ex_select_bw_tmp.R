## Case 1: normal error with known variance (misspecified
## distribution)
case1 <-
  kerdec_dens(Y, method = "CV", kernel = "sinc", h0 = c(.18, 0.26),
              lower = lower, upper = upper, h = 0.08,
              error_dist = "laplace",
              error_scale_par = sd_error, resolution = 64)

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "vp", h0 = c(.1, 0.5),
              lower = lower, upper = upper, h = 0.08,
              error_dist = "normal",
              error_scale_par = sd_error, resolution = 64)

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "triw", h0 = c(.09, 0.26),
              lower = lower, upper = upper, h = 0.08,
              error_dist = "normal",
              error_scale_par = sd_error, resolution = 64)

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "tric", h0 = c(.1, 0.26),
              lower = lower, upper = upper, h = 0.08,
              error_dist = "normal",
              error_scale_par = sd_error, resolution = 64)

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "flat", h0 = c(.1, 0.26),
              lower = lower, upper = upper, h = 0.08,
              error_dist = "normal",
              error_scale_par = sd_error, resolution = 64)


with(case1, plot(x, f_vals, type = "l"))
