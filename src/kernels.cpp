
#include <RcppArmadillo.h>
#include <fourierin.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//' Fourier transforms of kernels
//'
//' Explanation here.
//' @export
//[[Rcpp::export]]
double FT_ker(double t, int ker)
{
  double out;

  t = std::abs(t);

  if(t >= 1) return 0;

  switch(ker)
    {
    case 1:
      out = 1;
      break;
    case 2:
      out = 1 - t;
      break;
    case 3:
      out = std::pow(1 - t*t, 3.0);
      break;
    case 4:
      out = std::pow(1 - t*t*t, 3.0);
      break;
    case 5:
      if(t < 0.5) out = 1;
      else out = 1 - 2.0*(t - 0.5);
      break;
    default:
      Rcpp::stop("Kernel not specified.");
    }
  
   return out;
}


// //' Fourier transforms of kernels
// //'
// //' Explanation here.
// //' @export
// //[[Rcpp::export]]
// arma::vec FT_ker(arma::vec t, int ker)
// {
//   int i, n = t.n_rows;
//   vec out(n);
//   out = abs(t);

//   (out > 0).print("Greater than zero");
  
//  // case 1:
//  //   for(i = 0; i < n; i++)
//  //     { 
//  //     }
	 
//    return out;
// }


// /* We will compute the real and imaginary parts of the empirical
//    characteristic function as well as its modulus and the empirical
//    characteristic function itself. All of this for univariate and
//    multivariate cases.

//    The purpose of doing it separately it is just for avoiding
//    unnecesary computations or use of complex numbers.
//  */

// /* --------------------------------------------------------------- */
// /* -------- Real part of ecf ------------------------------------- */
// /* --------------------------------------------------------------- */

// //' Real part of empirical characteristic function
// //'
// //' Real part of empirical characteristic function of a d-dimensional
// //' random variable. This function is evaluated at m vectors of size
// //' d.
// //'
// //' This function must receive matrices. Vectors or values are not
// //' accepted.
// //' 
// //' @param t mxd matrix where the function will be evaluated.
// //' @param smp nxd matrix with sample size if size n.
// //'
// //' @return A vector of size m with the real part of the empirical
// //' characteristic function.
// //'
// //' @example
// //' examples/ex_ecf.R
// //'
// //' @export
// //[[Rcpp::export]]
// arma::vec ecf_re_cpp(arma::mat t, arma::mat smp)
// {
//   //  Display error dimensions are different in sample and
//   //  eval. points.
//   if(t.n_cols != smp.n_cols)
//     {
//       Rcpp::stop("t and smp must have the same number of columns");
//     }
  
//   return mean(cos(t * trans(smp)), 1);
// }

// /* --------------------------------------------------------------- */
// /* -------- Imaginary part of ecf -------------------------------- */
// /* --------------------------------------------------------------- */

// //' Imaginary part of empirical characteristic function
// //'
// //' Imaginary part of empirical characteristic function of a
// //' d-dimensional random variable. This function is evaluated at m
// //' vectors of size d.
// //'
// //' This function must receive matrices. Vectors or values are not
// //' accepted.
// //' 
// //' @param t mxd matrix where the function will be evaluated.
// //' @param smp nxd matrix with sample size if size n.
// //'
// //' @return A vector of size m with the imaginary part of the
// //' empirical characteristic function.
// //'
// //' @export
// //[[Rcpp::export]]
// arma::vec ecf_im_cpp(arma::mat t, arma::mat smp)
// {
//   //  Display error dimensions are different in sample and
//   //  eval. points.
//   if(t.n_cols != smp.n_cols)
//     {
//       Rcpp::stop("t and smp must have the same number of columns");
//     }

//   return mean(sin(t * trans(smp)), 1);
// }

// /* --------------------------------------------------------------- */
// /* -------- Modulus of ecf --------------------------------------- */
// /* --------------------------------------------------------------- */

// //' Modulus of empirical characteristic function
// //'
// //' Modulus of empirical characteristic function of a d-dimensional
// //' random variable. This function is evaluated at m vectors of size
// //' d.
// //'
// //' This function must receive matrices. Vectors or values are not
// //' accepted.
// //' 
// //' @param t mxd matrix where the function will be evaluated.
// //' @param smp nxd matrix with sample size if size n.
// //'
// //' @return A complex of size m with the modulus of the
// //' empirical characteristic function.
// //'
// //' @export
// //[[Rcpp::export]]
// arma::vec ecf_mod_cpp(arma::mat t, arma::mat smp)
// {
//   arma::vec real, imag;
//   arma::mat arg;
  
//   //  Display error dimensions are different in sample and
//   //  eval. points.
//   if(t.n_cols != smp.n_cols)
//     {
//       Rcpp::stop("t and smp must have the same number of columns");
//     }

//   arg = t * trans(smp);
//   real = mean(cos(arg), 1);
//   imag = mean(sin(arg), 1);
  
//   return sqrt(real % real + imag % imag);
// }
