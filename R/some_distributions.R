
##' Laplace Probability Density Function
##'
##' Evaluate a Laplace pdf in a vector
##' @param x Vector where the pdf will be evaluated
##' @param mean Mean of Laplace distribution
##' @param sd Standard deviation of Laplace distribution
##' @return A vector with pdf values
##' @author Guillermo Basulto-Elias
##' @seealso \code{\link{plaplace}} for cumulative distribution
##'     function and \code{\link{rlaplace}} for random number generation.
##' @examples
##' x <- seq(-2, 4, 0.05)
##' vals <- dlaplace(x, mean = 1, sd = 0.99)
##' plot(x, vals, type = "l")
##' @export
dlaplace <- function(x, mean = 0, sd = 1){
    b <- sd/sqrt(2)
    return(exp(-abs(x - mean)/b)/(2*b))
}

##' Laplace Cumulative Distribution Function
##'
##' Evaluate a Laplace cdf in a vector
##' @param x Vector where the cdf will be evaluated
##' @param mean Mean of Laplace distribution
##' @param sd Standard deviation of Laplace distribution
##' @return A vector with cdf values
##' @author Guillermo Basulto-Elias
##' @seealso \code{\link{dlaplace}} for probability density
##'     function and \code{\link{rlaplace}} for random number generation.
##' @examples
##' x <- seq(-2, 4, 0.05)
##' vals <- plaplace(x, mean = 1, sd = 0.99)
##' plot(x, vals, type = "l")
##' @export
plaplace <- function(x, mean = 0, sd = 1){
    b <- sd/sqrt(2)
    logic <- x < mean
    sign <- 2*logic - 1
    return((!logic) + sign*0.5*exp(sign*(x - mean)/b))
}

##' Laplace Random Generator
##'
##' Generate Laplace random sample
##' @param n Sample size
##' @param mean Mean of Laplace distribution
##' @param sd Standard deviation of Laplace distribution
##' @return A vector with random values
##' @author Guillermo Basulto-Elias
##' @seealso \code{\link{dlaplace}} for probability density
##'     function and \code{\link{plaplace}} for cumulative distribution
##'     function.
##' @examples
##' x <- seq(-2, 4, 0.05)
##' vals <- plaplace(x, mean = 1, sd = 0.99)
##' plot(x, vals, type = "l")
##' @export
rlaplace <- function(n, mean = 0, sd = 1){
    b <- sd/sqrt(2)
    out <- rexp(n, 1/b) - rexp(n, 1/b) + mean
    return(out)
}

##' Laplace Characteristic Function
##'
##' Obtain the Laplace characteristic function in a vector
##' @param t Vector where the characteristic function will be
##'     evaluated
##' @param mean Mean of Laplace distribution
##' @param sd Standard deviation of Laplace distribution
##' @return A complex vector with characteristic function values
##' @author Guillermo Basulto-Elias
##' @export
cflaplace <- function(t, mean = 0, sd = 1){
    b <- sd/sqrt(2)
    return(exp(1i*mu*t)/(1 + b^2*t^2))
}

##' Normal Characteristic Function
##'
##' Obtain the normal characteristic function in a vector
##' @param t Vector where the characteristic function will be
##'     evaluated
##' @param mean Mean of Laplace distribution
##' @param sd Standard deviation of Laplace distribution
##' @return A complex vector with characteristic function values
##' @author Guillermo Basulto-Elias
##' @export
cfnorm <- function(t, mean = 0, sd = 1){
    return(exp(1i*mu*t - 0.5*sd^2*t^2))
}

##' Log-likelihood of Laplace Convolution
##'
##' Negative value of the log-likelihood of the difference of two
##' independent Laplace random variables.
##' @param sigma Value or vector of positive values where the this
##'     function will be evaluated.
##' @param smp Sample of differences of Laplace random variables.
##' @return A vector of values of the log-likelihood.
##' @author Guillermo Basulto-Elias
laplace_convol_loglik <- function(sigma, smp){
                                        # Adapt function to receive
                                        # vectors
    if(length(sigma) > 1){
        out <- sapply(sigma, function(sig)
            laplace_convol_loglik(sig, smp))
        return(out)
    }

    n <- length(smp)
    x <- abs(smp)
    b <- sqrt(2)/sigma
    out <- n*log(b) + sum(log(1 + b*x)) - b*sum(x)
    return(-out)
}

##' Gradient of Log-likelihood of Laplace Convolution
##'
##' Negative value of the gradient of log-likelihood of the difference
##' of two independent Laplace random variables.
##' @param sigma Value or vector of positive values where the this
##'     function will be evaluated.
##' @param smp Sample of differences of Laplace random variables.
##' @return A vector of values of the gradient of the log-likelihood.
##' @author Guillermo Basulto-Elias
dlaplace_convol_loglik <- function(sigma, smp){
                                        # Adapt function to receive
                                        # vectors
    if(length(sigma) > 1){
        out <- sapply(sigma, function(sig)
            dlaplace_convol_loglik(sig, smp))
        return(out)
    }

    n <- length(smp)
    x <- sqrt(2)*abs(smp)
    s <- sigma
    
    out <- n/s + sum(x/(s^2 + s*x)) - sum(x)/s^2
    return(out)
}

##' MLE of the std. dev. for Laplace convolution
##'
##' Compute the maximum likelihood estimator for the standard
##' deviation of the difference of two independent Laplace random
##' variables.
##' @param smp Sample of differences of Laplace random variables.
##' @return The MLE of the standard deviation
##' @author Guillermo Basulto-Elias
mle_laplace_diffs <- function(smp){
    sig0 <- sd(smp)/sqrt(2)             # Use sample std. deviation as
                                        # initial value

    ## Minimize using L-BFGS-B providing the function to minimize and
    ## its gradient
    out <- optim(par = sig0,
                 fn = laplace_convol_loglik,
                 gr = dlaplace_convol_loglik,
                 lower = .Machine$double.eps, method = "L-BFGS-B",
                 smp = smp)
    return(out$par)
}

##' Compute scale parameter for error in kdde
##'
##' This function computes the maximum likelihood estimator of the
##' standard deviation for the error distribution provided as long as
##' this error is normal or Laplace and the sample is either a sample
##' or errors or a sample of differeces of errors (which arises in
##' panel data).
##'
##' This function also computes the standard deviation when the error
##' is approximated with the enpirical charcteristic function, since
##' it is useful for bandwidth selection.
##' 
##' @param error_dist An integer specifying the error distribution: 1
##'     if non-parametric, 2 if Laplace, 3 if normal.
##' @param error_smp Vector of sample of errors or differences of
##'     errors.
##' @param k An integer specifying what is being provided as sample of
##'     errors: '1' for pure sample of errors and otherwise for
##'     differences of errors.
##' @param error_scale_par NULL by default. It is is given, the
##'     function does nothing but return this value.
##' @return The (positive) standard deviation
##'
##' @examples
##' ## Difference of Laplace errors is provided:
##' sigma <- 5
##' smp <- rlaplace(n = 100, sd = sigma) - rlaplace(n = 100, sd = sigma)
##' kerdec:::compute_scale_par(2, smp, 2)
##' @author Guillermo Basulto-Elias
compute_scale_par <- function(error_dist, error_smp, k, error_scale_par = NULL){

    ## Estimate scale parameter from error_smp if required.
    if (is.null(error_smp)){
        ## Display error message if error distribution is not given
        ## and it cannot be approximated.
        if (is.null(error_scale_par) | error_dist == 1){
            stop(paste0("If a panel data structured are not given",
                        ", nor a sample of errors, then the error ",
                        "distribution must be specified as well ",
                        "as its scale parameter"))
            return(error_scale_par)
        }
    }
    
    if (k == 1) {                         # Pure sample of errors
        sigma <- switch(
            error_dist,
            sd(error_smp),             # ecf
            sqrt(2)*mean(abs(error_smp)), # Laplace
            sqrt(mean((error_smp)^2))     # Normal
        )
    } else {
        sigma <- switch(
            error_dist,
            sd(error_smp)/sqrt(2),                         #
            mle_laplace_diffs(error_smp),
            sqrt(mean((error_smp)^2)/2)
        )
    }
    return(sigma)
}

