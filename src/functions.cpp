
#include <RcppArmadillo.h>
#include <fourierin.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

/* 
   We start by defining functions to work with Fourier transforms of
   kernels. All of these transforms with support [-1, 1].
   Specifically:

   ft_kernel_cpp
*/

// ------------------------------------------------------------------
//             Evaluate at only one value
// ------------------------------------------------------------------
double ft_kernel_cpp(double t, int ker)
{
  double out;			// Output number
  
  t = std::abs(t);		// Only the abs. value is required.
  if(t >= 1) return 0;		// Ft of kernel must have [-1, 1] as
  // their support
  
  switch(ker)
  {
  case 1:			// Sinc kernel
    out = 1;
    break;
  case 2:			// VP kernel
    out = 1 - t;
    break;
  case 3:			// Triweight kernel
    out = std::pow(1 - t*t, 3.0);
    break;
  case 4:			// Tricube kernel
    out = std::pow(1 - t*t*t, 3.0);
    break;
  case 5:			//  Flat-top kernel
    if(t < 0.5) out = 1;
    else out = 1 - 2.0*(t - 0.5);
    break;
  default:
    Rcpp::stop("Kernel not specified.");
  }
  
  return out;
}

// ------------------------------------------------------------------
//             Evaluate at only at vectors or matrices
// ------------------------------------------------------------------
//' Fourier transforms of kernels
//'
//' See ft_kernel
//[[Rcpp::export]]
arma::vec ft_kernel_cpp(arma::mat t, int ker)
{
  int i, j, n = t.n_rows, d = t.n_cols;
  arma::vec out(n);
  
  out.ones();			// Start out with ones.
  
  for(i = 0; i < n; i++)
    for(j = 0; j < d; j++)
    {
      out(i) *= ft_kernel_cpp(t(i, j), ker);
    }
    
    return out;
}



/* 
   We will compute the real and imaginary parts of the empirical
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


/* Deconvolution formulas */

/* Deconvolution formula when a sample of pure errors is available,
   which is independent from the signal. */
arma::cx_vec kerdec_dens_pure(arma::mat smp, arma::mat error_smp,
			      arma::vec lower, arma::vec upper,
			      double cutoff, int resolution)
{

  arma::cx_vec out(smp.n_rows);
  
  out = fourierin::fourierin_1d_cpp(smp, 0, 1, 0, 1, 0, 0);
  
  return out;
}


//' Multivariate kernel deconvolution density estimator
//'
//' @export
//[[Rcpp::export]]
arma::cx_vec kerdec_dens(arma::vec smp)
{
  arma::cx_vec out(smp.n_rows);

  out = fourierin::fourierin_1d_cpp(smp, 0, 1, 0, 1, 0, 0);
  
  return out;
}
