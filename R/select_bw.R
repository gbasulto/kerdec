### Bandwidth Selection
###
### This function provides a bandwidth for kernel denvolvolution
### density estimator. This works for several deconvolution scenarios,
### specifically, when the error distribution is known, when a sample
### of pure errors is available instead and when a contaminated sample
### is available as panel data.
###
### See the vignette for more details.
select_bw <- function(smp,
                      method = c("CV", "NR"),
                      h0 = NULL,
                      error = NULL,
                      resolution = 128,
                      error_proc = c("all", "vs_first", "indep_pairs")[1],
                      panel_proc = c("keep_first", "take_aver")[1],
                      truncation_bound = NULL,
                      kernel = 1){
    return(0)
}
                      
