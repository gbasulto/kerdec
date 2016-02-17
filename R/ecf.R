#' @useDynLib kerdec
#' @importFrom Rcpp sourceCpp
NULL

#' Real part of empirical characteristic function
#'
#'
#' Real part of empirical characteristic function of a d-dimensional
#' random variable. This function is evaluated at m vectors of size
#' d.
#'
#' This function must receive matrices or vectors. It is a wrapper
#' function that allows more general inputs.
#'
#' @param t mxd matrix where the function will be evaluated.
#' @param smp nxd matrix with sample size if size n.
#'
#' @return A vector of size m with the real part of the empirical
#' characteristic function.
#'
#' @example
#' examples/ex_ecf.R
#'
#' @export
ecf_real <- function(t, smp){

    ## Verify arguments are numeric structures.
    if(!is.numeric(t) | !is.numeric(smp)){
        stop("t and smp must be numeric")
    }

    ## Verify they are either vectors or matrices
    aux <- abs(length(dim(t)) - 1)
    if(aux != 1) stop("t and smp must be vectors or matrices")
    
    ## Create index:
    ## 1: smp vector, t vector
    ## 2: smp matrix, t vector
    ## 3: smp vector, t matrix
    ## 4: smp matrix, t matrix
    idx <- 2*is.matrix(t) + is.matrix(smp) + 1

    switch(idx,
           "1" = {
               t <- matrix(t)
               smp <- matrix(smp)
           },
           "2" = {
               ## Make t row vector
               t <- matrix(t, nrow = 1)
           },
           "3" = {
               ## Make smp row vector
               smp <- matrix(smp, nrow = 1)
           },
           "4" = {
               NULL
           })
    return(ecf_re_cpp(t, smp))
}


## if(!is.vector(t)) stop("t must be vector or matrix")
## if(!is.vector(t)) stop("t must be vector or matrix")
