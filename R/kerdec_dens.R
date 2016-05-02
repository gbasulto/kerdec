
compute_scale_par <- function(error_dist, error_smp, k){
    if (k == 1) {                         # Pure sample of errors
        sigma <- switch(
            error_dist,
            -1,             # ecf
            sqrt(2)*mean(abs(error_smp)), # Laplace
            sqrt(mean((error_smp)^2))     # Normal
        )
    } else {
        sigma <- switch(
            error_dist,
            -1,                         #
            mle_laplace_diffs(error_smp),
            sqrt(mean((error_smp)^2)/2)
        )
    }
    return(sigma)
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
                        truncation_bound = NULL){

    ## Let us first state all the implemented distributions. We will
    ## check later that the arguments are valid.
    bw_methods <- c("cv", "nr")
    kernels <- c("sinc", "vp", "triw", "tric", "flat")
    error_dists <- c("none", "laplace", "normal")
    error_procs <- c("all", "vs_first", "indep_pairs")
    panel_procs <- c("keep_first", "take_aver")

    ## Check that the sample is numeric. If it is a vector, cast it to
    ## a matrix.
    if(is.numeric(smp)){
        if(is.vector(smp)) smp <- matrix(smp)
    } else{
        stop("smp must be numeric.")
    }

    ## Ask the sample size to be at least three and also obtain the
    ## number of repetitions.
    n <- nrow(smp)
    if(n < 3) stop("Sample size must be of at least 3.")
    k <- ncol(smp)

    ## Convert to lower case the argument "method" and then it is
    ## implemented.
    method <- tolower(method)
    if(!(method %in% bw_methods)){
        msg <- paste0(c("\n Method ", method,
                        " is not implemented. ",
                        "The current methods are:\n ",
                        paste0(bw_methods, collapse = "  ")))
        stop(msg)
    }

    ## Check 'kernel' is numeric, if not, assign it a numerica value
    ## and then verify it is in the list of implemented kernels.
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

    
    ## Check that the error sample is numeric. If it is a vector, cast
    ## it to a matrix.
    if (is.numeric(error_smp)) {
        if (is.vector(error_smp)) error_smp <- matrix(error_smp)
    } else {
        if (!is.null(error_smp)) {
            stop("'error_smp' must be numeric.")
        }
    }

    ## If data are provided in a panel structure, compute differences
    ## of errors to approximate the error distribution.
    if (k > 1) {
        diff_mthd <- match(error_proc, error_procs)
        smp_mthd <- match(panel_proc, panel_procs)
        if (!(diff_mthd %in% 1:length(error_procs))) {
            msg <- paste0(c("\nerror_proc '",
                            error_proc, "' is not implemented. ",
                            "The current error_procs are:\n ",
                            paste0(error_procs, collapse = "  "),
                            "\n\n See vignette for details."))
            stop(msg)
        }
        if (!(smp_mthd %in% 1:length(panel_procs))) {
            msg <- paste0(c("\npanel_proc '",
                            panel_proc, "' is not implemented. ",
                            "The current panel_procs are:\n ",
                            paste0(panel_procs, collapse = "  "),
                            "\n\n See vignette for details."))
            stop(msg)
        }
                                        # Rename vars. once we checked
                                        # they're valid
        panel_proc <- smp_mthd
        error_proc <- diff_mthd
                                        # Compute differences for
                                        # errors and compute the
                                        # sample to be used for decon.
        error_smp <- process_differences(smp, diff_mthd)
        smp <- switch(smp_mthd, smp[, 1], matrix(rowMeans(smp)))
    } else {
        panel_proc = 1
    }

    ## Estimate scale parameter from error_smp if required.
    if (is.null(error_smp)){
        ## Display error message if error distribution is not given
        ## and it cannot be approximated.
        if (is.null(error_scale_par) | error_dist == 1){
            stop(paste0("If a panel data structured are not given",
                        ", nor a sample of errors, then the error",
                        "distribution must be specified as well ",
                        "as its scale parameter"))
        }
    } else {
        error_scale_par <-
            compute_scale_par(error_dist, error_smp, k)
    }


    ## If error_smp was null, the error distribution must have been
    ## given and that was already verified above. Thus we can set
    ## other random values to it (such values will NOT be used within
    ## kerdec_dens_cpp since error_dist is either 1 or 2)
    if(is.null(error_smp)) error_smp <- matrix(0, 5, 1)
    
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
                h = h))
}

##
