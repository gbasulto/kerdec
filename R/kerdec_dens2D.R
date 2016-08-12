## -------------------------------------------------------------------

check_dims <- function(dim1, dim2){
### This is for bivariate deconvolution.

    if(dim1 == dim2) return(dim1)

    ## If dimensions are not equal, display an error.
    arg_name1 <- deparse(substitute(dim1))
    arg_name2 <- deparse(substitute(dim2))

    msg <- paste0(arg_name1, " = ", dim1, " and ",
                  arg_name2, " = ", dim2, ", but they must be equal.")
    stop(msg)

    return (NULL)
}



##' Bivariate Kernel Deconvolution Density Estimation
##'
##' This function provides a bandwidth for kernel denvolvolution
##' density estimator. This works for several deconvolution scenarios,
##' specifically, when the error distribution is known, when a sample
##' of pure errors is available instead and when a contaminated sample
##' is available as panel data.
##'
##' See the vignette for more details.
##' @param smp1 First coordinate of the sample. It is either a vector
##'     of size n or a nxl matrix (for repeated observations; l per
##'     each individual).
##' @param smp2 Second coordinate of the sample. It is either a vector
##'     of size n or a nxl matrix (for repeated observations; l per
##'     each individual).
##' @param method Method used to select the bandwidth. The currently
##'     available methods are "CV" (cross-validation) and "NR" (normal
##'     reference). Cross-validation is the default.
##' @param kernel Kernel whose Fourier transform has bounded
##'     support. There are currently 5 kernel programmed: "sinc",
##'     "triangular", "triw", "tric" and "flat". See vignette for further
##'     details and examples. Flat-top kernel is the default.
##' @param lower Lower limit where the density will be evaluated.
##' @param upper Upper limit where the density will be evaluated.
##' @param x_eval Not yet activated. grid of values where the density
##'     will be evaluated. If it is given, parameters 'lower' and
##'     'upper' will be omitted.
##' @param h Bandwidth parameter which is only required if method =
##'     NULL.
##' @param h0 Optional argument used as initial value to look for the
##'     optimal value.
##' @param error_smp1 Optional vector with the first coordinate of
##'     errors. It is necessary to approximate the error distribution
##'     if it is unknown.
##' @param error_smp2 Optional vector with the second coordinate of
##'     errors. It is necessary to approximate the error distribution
##'     if it is unknown.
##' @param error_dist Three possible values are accepted. c("Normal",
##'     "Laplace", "None").
##' @param error_scale_par1 Scale parameter matching the standard
##'     deviation of the first coordinate of the error. It is NULL by
##'     default and it is required if (and only if) error_dist is
##'     normal or Laplace and no sample of error is provided nor
##'     contaminated sample comes in panel structure.
##' @param error_scale_par2 Scale parameter matching the standard
##'     deviation of the second coordinate of the error. It is NULL by
##'     default and it is required if (and only if) error_dist is
##'     normal or Laplace and no sample of error is provided nor
##'     contaminated sample comes in panel structure.
##' @param resolution Number of points to approximate integral in
##'     inversion formula, also to estimate the density (if grid was
##'     not given).
##' @param error_proc This is required only for panel data
##'     structure. It refers to the way errors are processed. See
##'     \code{\link{process_differences}} for further details.
##' @param panel_proc It defined what will be the contaminated sample
##'     (only) for panel data. "keep_first" will use the first column
##'     as sample while "take_aver" will take the average of
##'     contaminated samples per individual.
##' @param truncation_bound It truncates the integrand in
##'     deconvolution formula when the denominator is smaller than
##'     this bound.
##' @param bw_interval Do not modify it.
##' @return A list
##' @author Guillermo Basulto-Elias
##' @export
kerdec_dens2D <- function(smp1, smp2,
                          method = c("CV", "NR")[1],
                          kernel = "flat",
                          lower = NULL, upper = NULL,
                          x_eval = NULL,
                          h = NULL, h0 = NULL,
                          error_smp1 = NULL,
                          error_smp2 = NULL,
                          error_dist = "None",
                          error_scale_par1 = NULL,
                          error_scale_par2 = NULL,
                          resolution = 128,
                          error_proc = "all",
                          panel_proc = "keep_first",
                          truncation_bound = NULL,
                          bw_interval = NULL)
{
    ## Let us first state all the implemented distributions. We will
    ## check later that the arguments are valid.
    bw_methods <- c("cv", "nr", "none")
    kernels <- c("sinc", "triangular", "triw", "tric", "flat")
    error_dists <- c("none", "laplace", "normal")
    error_procs <- c("all", "vs_first", "indep_pairs")
    panel_procs <- c("keep_first", "take_aver")

    smp1 <- check_smp(smp1)             # Check smp is num. matrix
    smp2 <- check_smp(smp2)
    
    n <- check_dims(nrow(smp1), nrow(smp2))
    k <- check_dims(ncol(smp1), ncol(smp2))
    ## Ask the sample size to be at least three and also obtain the
    ## number of repetitions.
    if(n < 3) stop("Sample size must be of at least 3.")

    method <- check_bw_method(method, bw_methods, h)
    kernel <- kernel2numeric(kernel, kernels)
    error_dist <- error_dist2numeric(error_dist, error_dists)
    error_smp1 <- check_error_smp(error_smp1)
    error_smp2 <- check_error_smp(error_smp2)
    panel_proc <- panel_proc2numeric(panel_proc, panel_procs)

    ## If data are provided in a panel structure, compute differences
    ## of errors to approximate the error distribution.
    if (k > 1) {
        error_proc <- error_proc2numeric(error_proc, error_procs)
        error_smp1 <- process_differences(smp1, error_proc)
        error_smp2 <- process_differences(smp2, error_proc)
        smp1 <- switch(panel_proc, smp1[, 1], matrix(rowMeans(smp1)))
        smp2 <- switch(panel_proc, smp2[, 1], matrix(rowMeans(smp2)))
    } 

    ## Compute error scale parameter if it was not given.
    error_scale_par1 <- compute_scale_par(error_dist, error_smp1, k,
                                         error_scale_par1)
    error_scale_par2 <- compute_scale_par(error_dist, error_smp2, k,
                                         error_scale_par2)
    
    ## Now we select the initial h0, if it was not provided.
    if(is.null(h0)){
        h0 <- ifelse(error_dist == 2,
        (5*error_scale_par1^4/n)^(1/9),
        error_scale_par1/sqrt(log(n)/2))
    }

    ## If error_smp was null, the error distribution must have been
    ## given and that was already verified above. Thus we can set
    ## other random values to it (such values will NOT be used within
    ## kerdec_dens_cpp since error_dist is either 1 or 2)
    if (is.null(error_smp1)) {
        if (!is.null(error_smp2)) {
            msg <- "Only one sample of errors was provided"
            stop(msg)
        }
        error_smp1 <- matrix(0, 5, 1)
        error_smp2 <- matrix(0, 5, 1)
    } 

    ## Define smp and error_smp, which are two-column matrices.
    smp <- cbind(smp1, smp2)
    error_smp <- cbind(error_smp1, error_smp2)
    error_scale_par <- c(error_scale_par1, error_scale_par2)

    vals <- kerdec_dens2D_cpp(smp, error_smp, h,
                              lower,
                              upper,
                              c(resolution, resolution),
                              kernel,
                              error_scale_par, k,
                              error_dist,
                              panel_proc,
                              cutoff = 999)
    
    return(list(smp = smp, error_smp = error_smp, vals = vals))
    ## ## 
    ## h_optim <-
    ##     switch(method,
    ##            nr = h_NR(h0, smp, error_smp, resolution,
    ##                      kernel, n, error_scale_par, k,
    ##                      error_dist, panel_proc, bw_interval),
    ##            cv = h_CV(h0, smp, error_smp, resolution,
    ##                      kernel, error_scale_par, k,
    ##                      error_dist, panel_proc, bw_interval),
    ##            none = NULL)
    
    ## if(is.null(h)) h <- h_optim$estimate

    ## ## Compute density values (if required).
    ## switch(is.null(lower) + is.null(upper) + 1,
    ## {
    ##     f_vals <-
    ##         kerdec_dens_cpp(smp = smp, error_smp = error_smp, h = h,
    ##                         lower = lower, upper = upper,
    ##                         resolution = resolution, ker = kernel,
    ##                         sigma = error_scale_par, k = k,
    ##                         error_dist = error_dist,
    ##                         panel_proc = panel_proc)
    ##     f_vals <- Re(f_vals)
    ##     x_eval <- seq(lower, upper, len = resolution + 1)[-resolution]
        
    ## },
    ## stop("'lower' or 'upper' arguments were not provided."),
    ## {
    ##     x_eval <- NULL
    ##     f_vals <- NULL
    ## })
    
    ## return(list(f_vals = f_vals,
    ##             x_eval= x_eval,
    ##             h = h,
    ##             h0 = h0,
    ##             h_optim = h_optim))
}

## ##' Select bandwidth
## ##'
## ##' This function computes a bandwidth with the selected method and it
## ##' also plots the corresponging function to minimize (MISE, AMISE,
## ##' CV, etc.).
## ##' @inheritParams kerdec_dens
## ##' @param bw_interval A bivariate vector with the limits where the
## ##'     function to minimize is going to be plotted.
## ##' @return A list with the initial value, the computed value and the
## ##'     result of the minimization algorithm.
## ##' @author Guillermo Basulto
## ##' @export
## select_bw <- function(smp1, smp2,
##                       method = c("CV", "NR")[1],
##                       kernel = "flat",
##                       h0 = NULL,
##                       error_smp = NULL,
##                       error_dist = "None",
##                       error_scale_par = NULL,
##                       resolution = 128,
##                       error_proc = "all",
##                       panel_proc = "keep_first",
##                       truncation_bound = NULL,
##                       bw_interval = NULL){
##     out <- kerdec_dens(smp, method, kernel, NULL, NULL, NULL, NULL,
##                        h0, error_smp, error_dist, error_scale_par,
##                        resolution, error_proc, panel_proc, NULL,
##                        bw_interval)
##     out <- out[3:5]
## }
