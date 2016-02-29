
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace arma;

/* In this file we define functions to evaluate Fourier transforms of
   kernels. Such Fourier transforms must have support [-1, 1].  */

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

