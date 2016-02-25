
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
//' @example
//' examples/ex_ecf.R
//[[Rcpp::export]]
arma::vec ecf_re_cpp(arma::mat t, arma::mat smp)
{
  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
    {
      Rcpp::stop("t and smp must have the same number of columns");
    }
  
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
//[[Rcpp::export]]
arma::vec ecf_im_cpp(arma::mat t, arma::mat smp)
{
  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
    {
      Rcpp::stop("t and smp must have the same number of columns");
    }

  return mean(sin(t * trans(smp)), 1);
}

/* --------------------------------------------------------------- */
/* -------- Modulus of ecf --------------------------------------- */
/* --------------------------------------------------------------- */

//' Modulus of empirical characteristic function
//'
//' Modulus of empirical characteristic function of a d-dimensional
//' random variable. This function is evaluated at m vectors of size
//' d.
//'
//' This function must receive matrices. Vectors or values are not
//' accepted.
//' 
//' @param t mxd matrix where the function will be evaluated.
//' @param smp nxd matrix with sample size if size n.
//'
//' @return A complex of size m with the modulus of the
//' empirical characteristic function.
//[[Rcpp::export]]
arma::vec ecf_mod_cpp(arma::mat t, arma::mat smp)
{
  arma::vec real, imag;
  arma::mat arg;
  
  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
    {
      Rcpp::stop("t and smp must have the same number of columns");
    }

  arg = t * trans(smp);
  real = mean(cos(arg), 1);
  imag = mean(sin(arg), 1);
  
  return sqrt(real % real + imag % imag);
}

/* --------------------------------------------------------------- */
/* -------- ecf -------------------------------------------------- */
/* --------------------------------------------------------------- */

//' Empirical characteristic function
//'
//' Empirical characteristic function of a d-dimensional
//' random variable. This function is evaluated at m vectors of size
//' d.
//'
//' This function must receive matrices. Vectors or values are not
//' accepted.
//' 
//' @param t mxd matrix where the function will be evaluated.
//' @param smp nxd matrix with sample size if size n.
//'
//' @return A complex vector of size m with the empirical
//' characteristic function evaluated at t.
//[[Rcpp::export]]
arma::cx_vec ecf_cpp(arma::mat t, arma::mat smp)
{
  arma::vec real, imag;
  arma::mat arg;
  
  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
    {
      Rcpp::stop("t and smp must have the same number of columns");
    }

  arg = t * trans(smp);
  real = mean(cos(arg), 1);
  imag = mean(sin(arg), 1);
  
  return arma::cx_vec(real, imag);
}

