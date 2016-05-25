##' Fourier transform of kernel
##'
##' Evaluate Fourier transform of kernels and product kernels.
##'
##' The Fourier transform (FT) of all the kernel functions considered
##' for this package have a bounded support, as it is usual in kernel
##' deconvolution methods. All the kernels are symmetric.
##'
##' See vignette for detailed examples.
##' @param t A vector of size n of a nxd matrix where the Fourier
##'     transform of a kernel or a d-product kernel (if t is a matrix)
##'     will be evaluated.
##' @param ker Character or number specifying the kernel. All the
##'     kernels are symmetric and they have a FT supported on [-1, 1].
##'     The possible choices
##'     are:
##'     \itemize{
##'       \item{1 or "sinc"}{ Sinc kernel}
##'       \item{2 or "triangular"}{ Kernel whose FT is a triangle}
##'       \item{3 or "triweight"}{ Kernel whose FT is proportional to
##'                               triweight kernel}
##'       \item{4 or "tricube"}{ Kernel whose FT is proportional to
##'                               tricube kernel}
##'       \item{5 or "flat-top"}{ A kernel with FT equal to one around
##'                              zero and decreasing linearly to zero
##'                              at -1/2 and 1/2 }}
##' @return A vector of size n.
##' @examples
##' plot(function(t) ft_kernel(t, "flat"), -1.5, 1.5)
##' @export
ft_kernel <- function(t, ker){
  ## If t is dataframe, convert it to matrix
  if(is.data.frame(t)) t <- as.matrix(t)
  
  ## Stop if t is not numeric
    if(!is.numeric(t)){
        stop("t must be numeric matrix or numeric vector")
    }

    ## Convert ker to numeric value if a character was provided.
    if(is.character(ker)){
        ker <- switch(ker,
                      "sinc" = 1,
                      "triangular" = 2,
                      "triweight" = 3,
                      "tricube" = 4,
                      "flat-top" = 5,
                      stop("Kernel not specified")
                      )
    }

        ## Call the appropriate function
        if(is.vector(t)) t <- matrix(t)

        return(drop(ft_kernel_cpp(t, ker)))
}

