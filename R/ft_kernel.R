##' Fourier transform of kernel
##'
##' Compute Fourier transform of product kernels.
##'
##' The Fourier transform of all the kernel functions considered for
##' this package have a bounded support, as it is usual in kernel
##' deconvolution methods.
##'
##' @param mat nxp matrix or vector of size n.
##' @param ker
##'
##' @return A vector of size n.
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

