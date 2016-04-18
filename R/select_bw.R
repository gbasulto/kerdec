## -------------------------------------------------------------------
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
                      method = c("CV", "NR")[1],
                      kernel = c("sinc", "vp", "triw", "tric",
                                 "flat")[5],
                      h0 = NULL,
                      error_dist = c("Normal", "Laplace", "None")[1],
                      error_scale_par = NULL,
                      error_smp = NULL,
                      resolution = 128,
                      error_proc = c("all", "vs_first",
                                     "indep_pairs")[1],
                      panel_proc = c("keep_first", "take_aver")[1],
                      truncation_bound = NULL){

    ## Check that the sample is numeric. If it is a vector, cast it to
    ## a matrix.
    if(is.numeric(smp)){
        if(is.vector(smp)) smp <- matrix(smp)
    } else{
        stop("smp must be numeric.")
    }

    n <- nrow(smp)
    if(n < 3) stop("Sample size must be of at least 3.")

    ## Convert to lower case the argument "method" and then it is
    ## implemented.
    method <- tolower(method)
    bw_methods <- c("cv", "nr")
    if(!(method %in% bw_methods)){
        msg <- paste0(c("Method ", method, " is not implemented. ",
                        "The current methods are:\n ",
                        paste0(bw_methods, collapse = "\n")))
        stop(msg)
    }
    
    return(0)
}
