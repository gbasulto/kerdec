
#include <RcppArmadillo.h>
#include <fourierin.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

/* -------------------------------------------------------------------

   TABLE OF CONTENTS

   The code on this file could be separated in the following topics of
   functions:

   1. Fourier transforms of kernels.
   2. Empirical characteristic functions.
   3. Univariate kdde.

--------------------------------------------------------------------*/

/* -------------------------------------------------------------------

   FOURIER TRANSFORMS OF KERNELS AND PRODUCT KERNELS

   We start by defining functions to work with Fourier transforms of
   kernels. All of these transforms with support [-1, 1].

--------------------------------------------------------------------*/

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
    Rcpp::stop("Kernel not defined.");
  }

  return out;
}

// ------------------------------------------------------------------
//             Evaluate at only at vectors or matrices
// ------------------------------------------------------------------

//[[Rcpp::export]]
arma::vec ft_kernel_cpp(const arma::mat & t, int ker)
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

/* -------------------------------------------------------------------

   EMPIRICAL CHARACTERISTIC FUNCTIONS

   We will compute the real and imaginary parts of the empirical
   characteristic function as well as its modulus and the empirical
   characteristic function itself. All of this for univariate and
   multivariate cases.

   The purpose of doing it separately it is just for avoiding
   unnecesary computations or use of complex numbers.

--------------------------------------------------------------------*/

// ------------------------------------------------------------------
//           Real part of ecf
// ------------------------------------------------------------------

//[[Rcpp::export]]
arma::vec ecf_re_cpp(const arma::mat & t,
		     const arma::mat & smp)
{
  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
  {
    Rcpp::stop("t and smp must have the same number of columns");
  }

  return mean(cos(t * trans(smp)), 1);
}

// ------------------------------------------------------------------
//           Imaginary part of ecf
// ------------------------------------------------------------------

//[[Rcpp::export]]
arma::vec ecf_im_cpp(const arma::mat & t,
		     const arma::mat & smp)
{
  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
  {
    Rcpp::stop("t and smp must have the same number of columns");
  }

  return mean(sin(t * trans(smp)), 1);
}

// ------------------------------------------------------------------
//           Modulus of ecf
// ------------------------------------------------------------------

//[[Rcpp::export]]
arma::vec ecf_mod_cpp(const arma::mat & t, const arma::mat & smp)
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
/*                ecf                                              */
/* --------------------------------------------------------------- */

//[[Rcpp::export]]
arma::cx_vec ecf_cpp(const arma::mat & t, const arma::mat & smp)
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

/* -------------------------------------------------------------------

   UNIVARIATE KERNEL DECONVOLUTION FORMULAS

   We include here functions for univariate kdde for knwon Gaussian or
   Lpalce errors, pure samples errors and panel data.

   -----------------------------------------------------------------*/

// -------------------------------------------------------------------
// Denominator
// -------------------------------------------------------------------

arma::vec dens_denominator(const arma::vec & t,
			   const arma::vec & smp,
			   double sigma, int k,
			   int error_dist,
			   int panel_proc)
{
  arma::vec out(t.n_rows);

  // Include case where \bar{Y}_{\cdot\cdot} is being taken as the
  // sample.
  if(panel_proc == 2)
    {
      out = dens_denominator(t/k, smp, sigma, k, error_dist, 1);
      out = arma::pow(out, k);
      return out;
    }

  switch(error_dist)
  {
  case 1:			// Nonparametric case
    if (k > 1) out = sqrt(ecf_mod_cpp(t, smp));
    else out = ecf_mod_cpp(t, smp);
    break;
  case 2:			// Laplace casee
    out = 1.0/(1 + (t % t)*sigma*sigma/2.0);
    break;
  case 3:			// Normal case
    out = exp(-sigma*sigma*(t % t)/2.0);
    break;
  default:
    Rcpp::stop("Error distribution not defined.");
  }

  return out;
}

// -------------------------------------------------------------------
//    Density estimation
// -------------------------------------------------------------------

//[[Rcpp::export]]
arma::cx_vec kerdec_dens_cpp(const arma::vec & smp,
			     const arma::vec & error_smp,
			     double h,
			     double lower, double upper,
			     int resolution,
			     int ker,
			     double sigma, int k,
			     int error_dist,
			     int panel_proc,
			     double cutoff = 999)
{
  int m = resolution, i;
  arma::vec t(m), denom(m);
  arma::cx_vec fun_vals(m), out(m);

  // If no cutoff is given, it is set to the one suggested by Neumann
  // (1997).
  if(cutoff == 999) cutoff = 1/sqrt(smp.n_rows);

  // Define grid where the integrand will be evaluated.
  t = arma::linspace<arma::mat>(-1.0/h, 1.0/h - 2.0/h/m, m);

  denom = dens_denominator(t, error_smp, sigma, k, error_dist,
			   panel_proc);
  fun_vals = (ecf_cpp(t, smp) % ft_kernel_cpp(h*t, ker))/denom;

  for(i = 0; i < m; i++)
    {
      if(denom[i] < cutoff) fun_vals[i] = 0;
    }

  out = fourierin::fourierin_cx_1d_cpp(fun_vals, -1.0/h, 1.0/h,
				    lower, upper, -1.0, -1.0);

  return out;
}

// --------------------------------------------------------------------
//    Bandwidth selection: Normal reference
// --------------------------------------------------------------------


//[[Rcpp::export]]
double amise(double h,
	     double mu2K2, double R,
	     const arma::vec & error_smp,
	     int resolution,
	     int ker, int n,
	     double sigma, int k,
	     int error_dist,
	     int panel_proc)
{
  double out, delta, h4 = h*h*h*h;
  int m = resolution;
  arma::vec aux(m);

  // aux is the m-sized grid vector
  delta = 2.0/m;
  aux = arma::linspace<arma::mat>(-1.0, 1.0 - delta, m);

  // aux is now the integrand in (3.1), Delaigle & Gijbels (2004)
  aux = ft_kernel_cpp(aux, ker)/
    dens_denominator(aux/h, error_smp, sigma, k, error_dist,
		     panel_proc);
  aux = arma::abs(aux);
  aux = (aux % aux);

  // out is equal to the term with the integral
  out = delta*arma::sum(aux)/(2*datum::pi*n*h);

  // Now we add the second term
  out += h4*mu2K2*R/4.0;

  return out;
}


// --------------------------------------------------------------------
//    Cross-validation bandwidth
// --------------------------------------------------------------------

//[[Rcpp::export]]
double CV(double h, const arma::vec & Z, const arma::vec & smp,
	  const arma::vec & error_smp, int resolution, int ker,
	  double sigma, int k, int error_dist, int panel_proc)
{
  /* 
     This function in the one involved in formula (1.7) from Youndje &
     Wells (2007), required by the cross-validation formula.
  */
  
  int m = resolution;
  arma::vec t(m), denom(m), fun_vals(m), cv_aux(m), kernel_vals(m);
  double st, fvals, out, delta = 2.0/h/m;	// Gridsize
  
  // Define grid where the integrand will be evaluated.
  t = arma::linspace<arma::mat>(-1.0/h, 1.0/h - delta, m);

  // ST_hat has the denominator squared.
  denom = dens_denominator(t, error_smp, sigma, k, error_dist,
			   panel_proc);
  denom = denom % denom;

  // kernel values are common for ST and the integral in (1.7) after
  // using Parseval's identity.
  kernel_vals = ft_kernel_cpp(h*t, ker);

  // ST integrand values
  fun_vals = (ecf_re_cpp(t, Z) % kernel_vals)/denom;

  // ST_hat in formula (1.7)
  st = sum(fun_vals)*delta/(2*datum::pi);
 
  // Integral of sqrd. f using Parseval's identity.
  cv_aux = ecf_mod_cpp(t, smp) % kernel_vals;
  cv_aux = cv_aux % cv_aux;	// square numerator
  cv_aux = cv_aux/denom;	// ... And divide by sq. denominator.
  fvals = sum(cv_aux)*delta/(2*datum::pi);

  out = (fvals - 2.0*st);
  
  return out;
}

// -------------------------------------------------------------------
// Pure error sample (independent from contaminated sample)
// -------------------------------------------------------------------

//' @export
//[[Rcpp::export]]
arma::cx_vec kerdec_dens_pure_1d_cpp(const arma::vec & smp,
				     const arma::vec & error_smp,
				     double h,
				     double lower, double upper,
				     int resolution,
				     int ker,
				     double cutoff = 999)
{
  int m = resolution, i;
  arma::vec t(m), denom(m);
  arma::cx_vec fun_vals(m), out(m);

  // If no cutoff is given, it is set to the one suggested by Neumann
  // (1997).
  if(cutoff == 999) cutoff = 1/sqrt(smp.n_rows);

  // Define grid where the integrand will be evaluated.
  t = arma::linspace<arma::mat>(-1.0/h, 1.0/h - 2.0/h/m, m);

  denom = ecf_mod_cpp(t, error_smp);
  fun_vals = ecf_cpp(t, smp) % ft_kernel_cpp(h*t, ker)/denom;

  for(i = 0; i < m; i++)
    {
      if(denom[i] < cutoff) fun_vals[i] = 0;
    }

  out = fourierin::fourierin_cx_1d_cpp(fun_vals, -1/h, 1/h,
				    lower, upper, -1.0, -1.0);

  return out;
}

//' Process differences for panel data
//'
//' Panel data allows to approximate the characteristic function of
//' the error taking differences of observations for each
//' individual. This way is not unique when there are more than two
//' replicates per individual. This function allows to do it in
//' several ways.
//'
//' @param smp n x d matrix
//' @param method Integer specifying method to process differences.
//'        1, all pairwise differences.
//'        2, all minus first.
//'        3, independent columns.
//' @export
//[[Rcpp::export]]
arma::vec process_differences(const arma::mat & smp, int method)
{
  // n is the sample size, d the dimension, l the output vector size
  // and idx index tro fill out output vector. The rest are indices
  // for loops.
  int n, d, i, j, k, l, idx;

  n = smp.n_rows;
  d = smp.n_cols;

  if(d == 1)
    {
      Rcpp::stop("smp must be a matrix with at least two columns.");
    }

  // Select vector size based on the selected method.
  switch(method)
    {
    case 1:			// All pairwise differences
      l = n*d*(d - 1)/2;
     break;
    case 2:			// All versus first.
      l = n*(d - 1);
      break;
    case 3:			// Group columns by pairs.
      l  = n*(d/2);
      break;
    default:
      Rcpp::stop("Differences method not defined.");
    }

  arma::vec out(l);
  idx = 0;

  switch(method)
    {
    case 1:			// All pairwise differences
      for(i = 0; i < n; i++)
	for(j = d - 1; j > 0; j--)
	  for(k = j - 1; k >= 0; k--)
	    {
	      out(idx) = smp(i, j) - smp(i, k);
	      idx++;
	    }
      break;
    case 2:			// All versus first.
      for(i = 0; i < n; i++)
	for(j = d - 1; j > 0; j--)
	  {
	    out(idx) = smp(i, j) - smp(i, 0);
	    idx++;
	  }
      break;
    case 3:			// Group columns by pairs.
      for(i = 0; i < n; i++)
	for(j = d/2 - 1; j >= 0; j--)
	  {
	    out(idx) = smp(i, 2*j + 1) - smp(i, 2*j);
	    idx++;
	  }
      break;
    default:
      Rcpp::stop("Differences method not defined.");
    }

  return out;
}

//' @export
//[[Rcpp::export]]
arma::vec error_cf_approx(const arma::vec & t,
			  const arma::mat & smp,
			  int diff_method)
{
  int k = smp.n_cols;
  arma::vec out(t);

  out = ecf_mod_cpp(t/k, process_differences(smp, diff_method));
  out = arma::pow(out, k/2.0);

  return out;

}

//' @export
//[[Rcpp::export]]
arma::cx_vec kerdec_dens_panel_1d_cpp(const arma::mat & smp,
				      double h,
				      double lower, double upper,
				      int resolution,
				      int ker,
				      double cutoff = 999,
				      int diff_processing = 1)
{
  // n in the sample size and l is the number of repetitions.
  int m = resolution, i, n;
  arma::vec t(m), denom(m);
  arma::cx_vec fun_vals(m), out(m);

  n = smp.n_rows;

  // If no cutoff is given, it is set to the one suggested by Neumann
  // (1997).
  if(cutoff == 999) cutoff = 1/sqrt(n);

  // Define grid where the integrand will be evaluated.
  t = arma::linspace<arma::mat>(-1.0/h, 1.0/h - 2.0/h/m, m);

  // Compute approximation to cf of averaged errors
  denom = error_cf_approx(t, smp, diff_processing);

  // Find values taking averaged obs. by row.
  fun_vals = ecf_cpp(t, mean(smp, 1)) % ft_kernel_cpp(h*t, ker)/denom;

  // Set integrand value to zero if the denominator is small.
  for(i = 0; i < m; i++)
    {
      if(denom[i] < cutoff) fun_vals[i] = 0;
    }

  out = fourierin::fourierin_cx_1d_cpp(fun_vals, -1/h, 1/h,
  				    lower, upper, -1.0, -1.0);

  return out;
}


