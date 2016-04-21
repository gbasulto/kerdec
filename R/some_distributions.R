
##' Laplace Probability Density Function
##'
##' Evaluate a Laplace pdf in a vector
##' @param x Vector where the pdf will be evaluated
##' @param mean Mean of Laplace distribution
##' @param sd Standard deviation of Laplace distribution
##' @return A vector with pdf values
##' @author Guillermo Basulto-Elias
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
