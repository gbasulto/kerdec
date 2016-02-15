
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

/* We will compute the real and imaginary parts of the empirical
   characteristic function as well as its modulus and the empirical
   characteristic function itself. All of this for univariate and
   multivariate cases.

   The purpose of doing it separately it is just for avoiding
   unnecesary computations or use of complex numbers.
 */

/* --------------------------------------------------------------- */
/* -------- Real part of ecf ------------------------------------- */
/* --------------------------------------------------------------- */

//' Real part of empirical characteristic function
//'
//' Real part of empirical characteristic function of a d-dimensional
//' random variable. This function is evaluated at m vectors of size
//' d.
//'
//' This function must receive matrices. Vectors or values are not
//' accepted.
//' 
//' @param t mxd matrix where the function will be evaluated.
//' @param smp nxd matrix with sample size if size n.
//'
//' @return A vector of size m with the real part of the empirical
//' characteristic function.
//'
//' @export
//[[Rcpp::export]]
arma::vec ecf_re(arma::mat t, arma::mat smp)
{
  return mean(cos(t * trans(smp)), 1);
}

/* --------------------------------------------------------------- */
/* -------- Imaginary part of ecf -------------------------------- */
/* --------------------------------------------------------------- */

//' Imaginary part of empirical characteristic function
//'
//' Imaginary part of empirical characteristic function of a
//' d-dimensional random variable. This function is evaluated at m
//' vectors of size d.
//'
//' This function must receive matrices. Vectors or values are not
//' accepted.
//' 
//' @param t mxd matrix where the function will be evaluated.
//' @param smp nxd matrix with sample size if size n.
//'
//' @return A vector of size m with the imaginary part of the
//' empirical characteristic function.
//'
//' @export
//[[Rcpp::export]]
arma::vec ecf_im(arma::mat t, arma::mat smp)
{
  return mean(sin(t * trans(smp)), 1);
}
