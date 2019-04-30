
#include <RcppArmadillo.h>
#include <fourierin.h>
#include <empichar.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace empichar;

/* --------------------------------------------------------------------

   TABLE OF CONTENTS

   The code on this file could be separated in the following topics of
   functions:

   1. Fourier transforms of kernels.
   2. Empirical characteristic functions.
   3. Univariate kdde.

---------------------------------------------------------------------*/

/* --------------------------------------------------------------------

   FOURIER TRANSFORMS OF KERNELS AND PRODUCT KERNELS

   We start by defining functions to work with Fourier transforms of
   kernels. All of these transforms with support [-1, 1].

---------------------------------------------------------------------*/

// --------------------------------------------------------------------
// Evaluate at only one value
// --------------------------------------------------------------------
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
  // t is nxd, where n is the number of points to evaluate and d the
  // dimension
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

   UNIVARIATE KERNEL DECONVOLUTION FORMULAS

   We include here functions for univariate kdde for knwon Gaussian or
   Laplace errors, pure samples errors and panel data.

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
//    Process differences for panel data
// -------------------------------------------------------------------


//[[Rcpp::export]]
arma::vec process_differences_cpp(const arma::mat & smp, int method)
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

  // Truncate integrand if denominator is very small.
  for(i = 0; i < m; i++)
    {
      if(denom[i] < cutoff) fun_vals[i] = 0;
    }

  out = fourierin::fourierin_cx_1d_cpp(fun_vals, -1.0/h, 1.0/h,
				    lower, upper, -1.0, -1.0);

  return out;
}

//[[Rcpp::export]]
arma::cx_vec kerdec_dens_nonreg_cpp(const arma::vec & smp,
				    const arma::vec & error_smp,
				    double h,
				    const arma::vec & x_eval,
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

  // Truncate integrand if denominator is very small.
  for(i = 0; i < m; i++)
    {
      if(denom[i] < cutoff) fun_vals[i] = 0;
    }

  out = fourierin::fourierin_cx_1d_nonregular_cpp(fun_vals,
						  -1.0/h, 1.0/h,
						  x_eval, resolution,
						  -1.0, -1.0);

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
// Cross-validation bandwidth
// --------------------------------------------------------------------

//[[Rcpp::export]]
double CV(double h, const arma::vec & Z, const arma::vec & smp,
	  const arma::vec & error_smp, int resolution, int ker,
	  double sigma, int k, int error_dist, int panel_proc,
	  double cutoff = 999)
{
  /* 
     This function in the one involved in formula (1.7) from Youndje &
     Wells (2007), required by the cross-validation formula.
  */
  
  int m = resolution;
  arma::vec t(m), denom(m), fun_vals(m), cv_aux(m), kernel_vals(m);
  double st, fvals, out, delta;

  // Gridsize
  delta = 2.0/h/m;		
  
  // If no cutoff is given, it is set to the one suggested by Neumann
  // (1997).
  if(cutoff == 999) cutoff = 1/sqrt(smp.n_rows);

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
  st = sum(fun_vals)*delta/(2.0*datum::pi);
 
  // Integral of sqrd. f using Parseval's identity.
  cv_aux = ecf_mod_cpp(t, smp) % kernel_vals;
  cv_aux = cv_aux % cv_aux;	// square numerator
  cv_aux = cv_aux/denom;	// ... And divide by sq. denominator.

  // Add sq. integral values.
  fvals = sum(cv_aux)*delta/(2.0*datum::pi);
  
  out = fvals - 2.0*st;
  
  return out;
}


/* -------------------------------------------------------------------

   BIVARIATE KERNEL DECONVOLUTION FORMULAS

   We include here functions for bivariate kdde for knwon Gaussian or
   Laplace errors, pure samples errors and panel data.

   -----------------------------------------------------------------*/

// -------------------------------------------------------------------
// Bivariate Denominator
// -------------------------------------------------------------------

arma::vec dens_denominator2D(const arma::mat & t,
			     const arma::mat & smp,
			     const arma::vec & sigma,
			     int k, int error_dist,
			     int panel_proc)
{
  int i, n = t.n_rows;
  arma::vec out(n);
  double aux1, aux2;
  
  // Include case where \bar{Y}_{\cdot\cdot} is being taken as the
  // sample.
  if(panel_proc == 2)
    {
      out = dens_denominator2D(t/k, smp, sigma, k, error_dist, 1);
      out = arma::pow(out, k);
      return out;
    }
  
  switch(error_dist)
    {
    case 1:			// Nonparametric case
      if (k > 1) out = sqrt(ecf_mod_cpp(t, smp));
      else out = ecf_mod_cpp(t, smp);
      break;
    case 2:			// Laplace case
      for (i = 0; i < n; i++)
	{
	  aux1 = 1.0 + t(i, 0)*t(i, 0)*sigma(0)*sigma(0)/2.0;
	  aux2 = 1.0 + t(i, 1)*t(i, 1)*sigma(1)*sigma(1)/2.0;
	  out(i) = 1.0/(aux1*aux2);
	}
      break;
    case 3:			// Normal case
      out = (t % t)*(sigma % sigma);
      out = exp(-out/2.0);
      break;
    default:
      Rcpp::stop("Error distribution not defined.");
  }
  
  return out;
}

// -------------------------------------------------------------------
//    Bivariate Density estimation
// -------------------------------------------------------------------

//[[Rcpp::export]]
arma::cx_mat kerdec_dens2D_cpp(const arma::mat & smp,
			       const arma::mat & error_smp,
			       double h,
			       const arma::vec & lower,
			       const arma::vec & upper,
			       const arma::vec & resolution,
			       int ker,
			       const arma::vec & sigma, int k,
			       int error_dist,
			       int panel_proc,
			       double cutoff = 999)
{
  // resolution.print("resolution = ");
  // Rprintf("Flag 1, %d.\n", resolution(0));
  int m1 = resolution(0), m2 = resolution(1), i, j;
  arma::vec t1(m1), t2(m2), aux_col(m1), denom(m1);
  arma::mat integrand(m1, m2), t_temp(m1, 2);
  arma::cx_mat out(m1, m2);

  // If no cutoff is given, it is set to the one suggested by Neumann
  // (1997).
  if (cutoff == 999) cutoff = 1/sqrt(smp.n_rows);

  // Define grid where the integrand will be evaluated.
  t1 = arma::linspace<arma::mat>(-1.0/h, 1.0/h - 2.0/h/m1, m1);
  t2 = arma::linspace<arma::mat>(-1.0/h, 1.0/h - 2.0/h/m2, m2);

  // Compute denominator column by column. denom is m1 x m2
  for (i = 0; i < m2; i++)
    {
      // Fill a column vector with the i-th entry of t2.
      aux_col.fill(t2(i));
      // aux_mat is m1 x 2 matrix with grid values.
      t_temp = arma::join_horiz(t1, aux_col);
      // Compute denominator
      denom = dens_denominator2D(t_temp, error_smp, sigma, k, error_dist,
				 panel_proc);
      // We will fill out this m1 sized complex vector
      out.col(i) = 
	(ecf_cpp(t_temp, smp) % ft_kernel_cpp(t_temp, ker))/
	denom;
      // Truncate function if denominator is very small
      for (j = 0; j < m1; j++)
      	{
	  if (denom(i) < cutoff) out(j, i) = 0;
      	}
    }
  
  out = fourierin::fourierin_cx_2d_cpp(out, -ones(2)/h, ones(2)/h,
				       lower, upper, -1.0, -1.0);

  return out;
}

// --------------------------------------------------------------------
// Cross-validation bandwidth
// --------------------------------------------------------------------

//[[Rcpp::export]]
double CV2D(double h, const arma::vec & Z, const arma::vec & smp,
	    const arma::vec & error_smp, int resolution, int ker,
	    double sigma, int k, int error_dist, int panel_proc,
	    double cutoff = 999)
{
  /* 
     This function in the one involved in formula (1.7) from Youndje &
     Wells (2007), required by the cross-validation formula.
  */
  
  int m = resolution;
  arma::vec t(m), denom(m), fun_vals(m), cv_aux(m), kernel_vals(m);
  double st, fvals, out, delta;
  
  // Gridsize
  delta = 2.0/h/m;		
  
  // If no cutoff is given, it is set to the one suggested by Neumann
  // (1997).
  if(cutoff == 999) cutoff = 1/sqrt(smp.n_rows);
  
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
  st = sum(fun_vals)*delta/(2.0*datum::pi);
 
  // Integral of sqrd. f using Parseval's identity.
  cv_aux = ecf_mod_cpp(t, smp) % kernel_vals;
  cv_aux = cv_aux % cv_aux;	// square numerator
  cv_aux = cv_aux/denom;	// ... And divide by sq. denominator.

  // Add sq. integral values.
  fvals = sum(cv_aux)*delta/(2.0*datum::pi);
  
  out = fvals - 2.0*st;
  
  return out;
}



// // --------------------------------------------------------------------
// //    Bandwidth selection: Normal reference
// // --------------------------------------------------------------------


// //[[Rcpp::export]]
// double amise(double h,
// 	     double mu2K2, double R,
// 	     const arma::vec & error_smp,
// 	     int resolution,
// 	     int ker, int n,
// 	     double sigma, int k,
// 	     int error_dist,
// 	     int panel_proc)
// {
//   double out, delta, h4 = h*h*h*h;
//   int m = resolution;
//   arma::vec aux(m);

//   // aux is the m-sized grid vector
//   delta = 2.0/m;
//   aux = arma::linspace<arma::mat>(-1.0, 1.0 - delta, m);

//   // aux is now the integrand in (3.1), Delaigle & Gijbels (2004)
//   aux = ft_kernel_cpp(aux, ker)/
//     dens_denominator(aux/h, error_smp, sigma, k, error_dist,
// 		     panel_proc);
//   aux = arma::abs(aux);
//   aux = (aux % aux);

//   // out is equal to the term with the integral
//   out = delta*arma::sum(aux)/(2*datum::pi*n*h);

//   // Now we add the second term
//   out += h4*mu2K2*R/4.0;

//   return out;
// }


