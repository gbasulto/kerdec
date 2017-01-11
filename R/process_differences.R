##' Process differences for panel data
##'
##' Panel data allows to approximate the characteristic function of
##' the error taking differences of observations for each
##' individual. This way is not unique when there are more than two
##' replicates per individual. This function allows to do it in
##' several ways.
##'
##' @param smp n x d matrix
##' @param method Integer or word specifying method to process
##'     differences.
##'        1 or "all", all pairwise differences.
##'        2 or "vs_first", all minus first.
##'        3 or "indep_pairs", independent columns.
##' @param na.rm Logical value specifying if NAs will be removed.
##' @export
process_differences <- function(smp, method, na.rm = TRUE){
    methods <- c("all", "vs_first", "indep_pairs")
    if (!is.numeric(method)) {
        method <- match(method, methods)
    }
    
    out <- process_differences_cpp(smp, method)
    
    if (na.rm) out <- out[!is.na(out)]
    
    return (out)
}
