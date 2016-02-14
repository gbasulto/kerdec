#include <Rcpp.h>
using namespace Rcpp;

/****************************************************************/
/***** ecf and ******************/
/****************************************************************/

double ecf_mod_1d(double t, NumericVector smp)
{
  /*
     Description: This function computes the module of the ecf from a
     sample smp at t.
   */
  int n = smp.size();
  NumericVector aux(n);
  double real_sq, imag_sq, out;

  // First find the real and imaginary parts and then the modulus.
  aux = t*smp;
  real_sq = pow(mean(cos(aux)), 2.0);
  imag_sq = pow(mean(sin(aux)), 2.0);
  out = sqrt(real_sq + imag_sq);

  return out;
}

//' Module of the empirical characteristic function
//'
//' It computes the module of the empirical characteristic function
//' evaluated on a sample
//' @export
// [[Rcpp::export]]
NumericVector ecf_mod_1d(NumericVector t, NumericVector smp)
{
  /*
     Description: This function computes the module of the ecf from a
     sample smp and it evaluates it in the vector t. It returns a
     vector of the same size as t.
   */
  int i, m = t.size();
  NumericVector out(m);

  for(i = 0; i < m; i++) out[i] = ecf_mod_1d(t[i], smp);

  return out;
}

