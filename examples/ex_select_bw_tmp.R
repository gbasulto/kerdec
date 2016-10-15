## Case 1: normal error with known variance (misspecified
## distribution)
case1 <-
  kerdec_dens(Y, method = "CV", kernel = "sinc",
              lower = lower, upper = upper,
              error_dist = "laplace",
              error_scale_par = sd_error,
              # bw_interval = c(0.08, 0.5),
              resolution = 128)

aaa <- select_bw(Y, "CV", "flat", error_smp = eps,
                 error_dist = "laplace",
                 bw_interval = c(0.08, 0.5))
aaa

with(case1, plot(x_eval, f_vals, type = "l"))

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "vp",
              lower = lower, upper = upper,
              error_dist = "normal",
              error_scale_par = sd_error,
#              bw_interval = c(0.08, 0.5),
              resolution = 128)

with(case1, plot(x_eval, f_vals, type = "l"))

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "triw",
              lower = lower, upper = upper,
              error_dist = "normal",
              error_scale_par = sd_error,
              # h0 = 0.19,
              bw_interval = c(0.07, 0.2),
              resolution = 128)

with(case1, plot(x_eval, f_vals, type = "l"))

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "tric",
              lower = lower, upper = upper,
              error_dist = "normal",
              error_scale_par = sd_error,
              # bw_interval = c(0.08, 0.5), 
              resolution = 128)

with(case1, plot(x_eval, f_vals, type = "l"))

case1 <-
  kerdec_dens(Y, method = "CV", kernel = "flat",
              lower = lower, upper = upper, h = NULL,
              error_dist = "normal",
              error_scale_par = sd_error,
#              bw_interval = c(0.08, 0.5), 
              resolution = 128)

with(case1, plot(x_eval, f_vals, type = "l"))

