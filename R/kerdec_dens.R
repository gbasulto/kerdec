## -------------------------------------------------------------------

##' Find positive minimizer with the input function
##'
##' Most methods used to find optimal bandwidths require to minimize a
##' function (tipically MISE, AMISE, CV, ...). The argument ought to
##' be positive. The fastest method that I empirically found was
##' nlm. If the function is evaluated at a negative number or the
##' output is not finite or not a number, the function will return the
##' maximum possible number, so nlm can work better with the input
##' function.
##' 
##' @param f Function to be minimized. Tipically a version of CV, MISE
##'     or AMISE.
##' @param h0 Initial bandwodth value.
##' @param ... Arguments passed to f besides the first argument, which
##'     must be the bandwidth.
##' @return The output of the nlm function.
##' @author Guillermo Basulto-Elias
optimize_bw <- function(f, h0, ...){
    ## Function to be minimized
    objective <- function(bw){
        if(bw < 0) return (.Machine$double.xmax) # Check bw > 0
        val <- f(bw, ...)
        if (is.nan(val) | !is.finite(val)){ # Check val is finite
            val <- .Machine$double.xmax
        }
        return (val)
    }
    
    ## Compute optimal bandwidth
    out <- nlm(objective, h0)
    
    ## Return warning if it is likely that the function did not find
    ## the minimum.
    if (!(out$code %in% 1:3)) {
        msg <-
            paste("\n\nThe nlm function used to find the optimal",
                  "value might have not converged appropriately.",
                  "Try another initial value (argument h0) or plot",
                  "the function to be minimized by providing the",
                  "limits in the argument bw_limits.\n\n")
        warning(msg)
    }
    return(out)
}

h_NR <- function(h0, smp, error_smp, resolution, kernel, n,
                 error_scale_par, k, error_dist, panel_proc,
                 bw_interval){
    ## Normal references works only for kernels with second moment. We
    ## check that here.
    if(!(kernel %in% 3:4)){
        stop("'nr' does not work for that kernel")
    }
    
    ## (1) Assign the corresponding kernel, (2) provide an estimate
    ## for sigma_X and provide an approximation to the roughness of
    ## the second derivative.
    mu2K2 <- ifelse(kernel == 3, 6^2, (4.822182e-05)^2)
    sigY <- sd(smp)
    sigE <- error_scale_par
    if(panel_proc == 2) sigE <- sigE/sqrt(k)
    sig_hat <- sqrt(sigY^2 - sigE^2)
    R <- 0.37/(sqrt(pi)*sig_hat^5)

    ## Function to be minimized
    amise_fun <- function(bw){
        if(bw < 0) return (.Machine$double.xmax) # Check bw > 0
        out <- amise(bw, mu2K2, R, error_smp, resolution, kernel, n,
                     error_scale_par, k, error_dist, panel_proc)
        if (is.nan(out) | !is.finite(out)){
            out <- .Machine$double.xmax 
            }
        return (out)
    }
    
    ## Compute optimal bandwidth
    h_optim <- nlm(amise_fun, h0)

    if (!(h_optim$code %in% 1:2)) {
        msg <-
            paste("The NR rule might have not found the optimal",
                  "bandwidth. We recommend to provide the argument",
                  "bw_interval (a vector of size 2 with the limits)",
                  "to plot the function to minimize")
        warning(msg)
    }

    ## If required, compute and print an interval with the values of
    ## the function to be minimized.
    if(!is.null(bw_interval)){

        h_grid <- seq(from = bw_interval[1], # Grid to plot vals.
                      to = bw_interval[2],
                      length.out = 100) 
        nr_vals <-                      # Vals
            sapply(h_grid, amise_fun)
        
        h <- h_optim$estimate
    
        plot(h_grid, nr_vals, type = "l")
        abline(v = h, col = "red")
        abline(v = h0, col = "green")
        cat("h0 = ", h0, "\n")
        cat("h = ", h, "\n")
    }
    
    return (h_optim)
}

h_CV <- function(h0, smp, error_smp, resolution, kernel,
                 error_scale_par, k, error_dist, panel_proc,
                 bw_interval){
    ## Vector of differences required for CV in Youndje
    ## (2007)
    Z <- process_differences(matrix(smp, nrow = 1),
                             method = 1)
    
    ## Function to be minimized and displayed.
    cv_fun <- function(bw){
        if(bw < 0) return (.Machine$double.xmax)
        
        cv_val <- 
            CV(bw, Z, smp, error_smp, resolution, kernel,
                          error_scale_par, k, error_dist, panel_proc)
        if (is.nan(cv_val) | !is.finite(cv_val)){
            cv_val <- .Machine$double.xmax
            }
        return(cv_val)
    }
    
    h_optim <- nlm(cv_fun, h0)
    
    if(!(h_optim$code %in% 1:2)){
        msg <-
            paste("The CV might have not found the opt. bandwidth.",
                  "We recommend to provide the argument bw_interval",
                  "(a vector of size 2 with the limits) to plot the",
                  "function to minimize")
        warning(msg)
    }
    
    ## If required, compute and print an interval with the values of
    ## the function to be minimized.
    if(!is.null(bw_interval)){

        h_grid <- seq(from = bw_interval[1],
                      to = bw_interval[2],
                      length.out = 100)
        
        cv_vals <-
            sapply(h_grid, cv_fun)
        
        h <- h_optim$estimate
    
        plot(h_grid, cv_vals, type = "l")
        abline(v = h, col = "red")
        abline(v = h0, col = "green")
        cat("h0 = ", h0, "\n")
        cat("h = ", h, "\n")
        
    }

    return(h_optim)
}

error_dist2numeric <- function (error_dist, error_dists) {
### This function checks that the specified error distribution is
### programmed and it converts it (from character) to numeric,
### being the number determined by error_dists.

    ## Check that the error distribution is valid and convert it to
    ## numeric argument.
    error_dist0 <- error_dist           # Create copy to display msg
    error_dist <- tolower(error_dist)   # To lower case
    error_dist <- match(error_dist, error_dists) # To numetric value
    if (!(error_dist %in% 1:length(error_dists))) {
        msg <- paste0(c("\nError distribution '",
                        error_dist0, "' is not implemented. ",
                        "The current error distributions are:\n ",
                        paste0(error_dists, collapse = "  "),
                        "\n\n See vignette for details."))
        stop(msg)
    }
    return (error_dist)
}

error_proc2numeric <- function (error_proc, error_procs) {
### This function checks that the specified method for computing
### differences of errors in panel data distribution is programmed and
### it converts it (from character) to numeric, being the number
### determined by error_procs.

    diff_mthd <- match(error_proc, error_procs)
    if (!(diff_mthd %in% 1:length(error_procs))) {
        msg <- paste0(c("\nerror_proc '",
                            error_proc, "' is not implemented. ",
                        "The current error_procs are:\n ",
                        paste0(error_procs, collapse = "  "),
                        "\n\n See vignette for details."))
        stop(msg)
    }
    return(diff_mthd)
}

panel_proc2numeric <- function (panel_proc, panel_procs){
### This function checks that the specified method for obtaining the
### ccontaminated sample (keeping only the firt column or taking the
### average by individual) is correct and also convert it to numeric
### argument, matching the order from panel_procs argument.
    smp_mthd <- match(panel_proc, panel_procs)
    if (!(smp_mthd %in% 1:length(panel_procs))) {
            msg <- paste0(c("\npanel_proc '",
                            panel_proc, "' is not implemented. ",
                            "The current panel_procs are:\n ",
                            paste0(panel_procs, collapse = "  "),
                            "\n\n See vignette for details."))
            stop(msg)
    }
    return (smp_mthd)
}

kernel2numeric <- function (kernel, kernels){
### It checks that the kernel is valid and convert it to numeric,
### matching the number of entry in 'kernels' argument.
    kernel0 <- kernel
    if (is.character(kernel)) kernel <- match(kernel, kernels)
    if (!(kernel %in% 1:length(kernels))) {
        msg <- paste0(c("\nKernel ",
                        kernel0, " is not implemented. ",
                        "The current kernels are:\n ",
                        paste0(kernels, collapse = "  "),
                        "\n\n See vignette for details."))
        stop(msg)
    }

    return (kernel)
    }
    

check_error_smp <- function (error_smp){
### This function verifies/converts error_smp to matrix.
    
    ## Check that the error sample is numeric. If it is a vector, cast
    ## it to a matrix.
    if (is.numeric(error_smp)) {
        if (is.vector(error_smp)) error_smp <- matrix(error_smp)
    } else {
        if (!is.null(error_smp)) {
            stop("'error_smp' must be numeric.")
        }
    }
    return (error_smp)
}

check_bw_method <- function (method, bw_methods, h){
### This function checks that the bandwidth selection method is
### implemented (from 'methods' argument). Argument 'h' is usually
### null.
    method <- tolower(method)
    if(!is.null(h)) method <- "none"    # Select "none" as bandwidth
                                        # sel. if h was given.
    if(!(method %in% bw_methods)){
        msg <- paste0(c("\n Method ", method,
                        " is not implemented. ",
                        "The current methods are:\n ",
                        paste0(bw_methods, collapse = "  ")))
        stop(msg)
    }
    return (method)
}

##' Kernel Deconvolution Density Estimation
##'
##' This function provides a bandwidth for kernel denvolvolution
##' density estimator. This works for several deconvolution scenarios,
##' specifically, when the error distribution is known, when a sample
##' of pure errors is available instead and when a contaminated sample
##' is available as panel data.
##'
##' See the vignette for more details.
##' @param smp It is either a vector of size n or a nxl matrix (for
##'     repeated observations; l per each individual).
##' @param method Method used to select the bandwidth. The currently
##'     available methods are "CV" (cross-validation) and "NR" (normal
##'     reference). Cross-validation is the default.
##' @param kernel Kernel whose Fourier transform has bounded
##'     support. There are currently 5 kernel programmed: "sinc",
##'     "vp", "triw", "tric" and "flat". See vignette for further
##'     details and examples. Flat-top kernel is the default.
##' @param lower Lower limit where the density will be evaluated.
##' @param upper Upper limit where the density will be evaluated.
##' @param h Bandwidth parameter which is only required if method =
##'     NULL.
##' @param h0 Optional argument used as initial value to look for the
##'     optimal value.
##' @param error_smp Optional vector errors. It is necessary to
##'     approximate the error distribution if it is unknown.
##' @param error_dist Three possible values are accepted. c("Normal",
##'     "Laplace", "None").
##' @param error_scale_par Scale parameter matching the standard
##'     deviation. It is NULL by default and it is required if (and
##'     only if) error_dist is normal or Laplace and no sample of
##'     error is provided nor contaminated sample comes in panel
##'     structure.
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
kerdec_dens <- function(smp,
                        method = c("CV", "NR")[1],
                        kernel = "flat",
                        lower, upper, h = NULL,
                        h0 = NULL,
                        error_smp = NULL,
                        error_dist = "None",
                        error_scale_par = NULL,
                        resolution = 128,
                        error_proc = "all",
                        panel_proc = "keep_first",
                        truncation_bound = NULL,
                        bw_interval = NULL){

    ## Let us first state all the implemented distributions. We will
    ## check later that the arguments are valid.
    bw_methods <- c("cv", "nr", "none")
    kernels <- c("sinc", "vp", "triw", "tric", "flat")
    error_dists <- c("none", "laplace", "normal")
    error_procs <- c("all", "vs_first", "indep_pairs")
    panel_procs <- c("keep_first", "take_aver")

    ## Check that the sample is numeric. If it is a vector, cast it to
    ## a matrix.
    if (!is.numeric(smp)) stop("smp must be numeric.")
    if(is.vector(smp)) smp <- matrix(smp)

    ## Ask the sample size to be at least three and also obtain the
    ## number of repetitions.
    n <- nrow(smp)
    if(n < 3) stop("Sample size must be of at least 3.")
    k <- ncol(smp)

    method <- check_bw_method(method, bw_methods, h)
    kernel <- kernel2numeric(kernel, kernels)
    error_dist <- error_dist2numeric(error_dist, error_dists)
    error_smp <- check_error_smp(error_smp)
    panel_proc <- panel_proc2numeric(panel_proc, panel_procs)

    ## If data are provided in a panel structure, compute differences
    ## of errors to approximate the error distribution.
    if (k > 1) {
        error_proc <- error_proc2numeric(error_proc, error_procs)
        error_smp <- process_differences(smp, error_proc)
        smp <- switch(smp_mthd, smp[, 1], matrix(rowMeans(smp)))
    } 

    ## Compute error scale parameter if it was not given.
    error_scale_par <- compute_scale_par(error_dist, error_smp, k,
                                         error_scale_par)
    
    ## Now we select the initial h0, if it was not provided.
    if(is.null(h0)){
        h0 <- ifelse(error_dist == 2,
        (5*error_scale_par^4/n)^(1/9),
        error_scale_par/sqrt(log(n)/2))
    }
    
    ## If error_smp was null, the error distribution must have been
    ## given and that was already verified above. Thus we can set
    ## other random values to it (such values will NOT be used within
    ## kerdec_dens_cpp since error_dist is either 1 or 2)
    if(is.null(error_smp)) error_smp <- matrix(0, 5, 1)

    ## 
    h_optim <-
        switch(method,
               nr = h_NR(h0, smp, error_smp, resolution,
                         kernel, n, error_scale_par, k,
                         error_dist, panel_proc, bw_interval),
               cv = h_CV(h0, smp, error_smp, resolution,
                         kernel, error_scale_par, k,
                         error_dist, panel_proc, bw_interval),
               none = NULL)
    
    if(is.null(h)) h <- h_optim$estimate
    
    f_vals <-
        kerdec_dens_cpp(smp = smp, error_smp = error_smp, h = h,
                        lower = lower, upper = upper,
                        resolution = resolution, ker = kernel,
                        sigma = error_scale_par, k = k,
                        error_dist = error_dist,
                        panel_proc = panel_proc)

    x <- seq(lower, upper, length.out = resolution + 1)[-resolution]
    
    return(list(f_vals = Re(f_vals),
                x = x,
                h = h,
                h0 = h0,
                h_optim = h_optim))
}
