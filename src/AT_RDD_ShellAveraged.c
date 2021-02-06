/**
 * @brief Radial Dose Distribution models
 */

/*
 *    AT_RDD_ShellAveraged.c
 *    ========
 *
 *    Created on: 05.04.2010
 *    Creator: kongruencja
 *
 *    Copyright 2006, 2010 The libamtrack team
 *
 *    This file is part of the AmTrack program (libamtrack.sourceforge.net).
 *
 *    AmTrack is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    AmTrack is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with AmTrack (file: copying.txt).
 *    If not, see <http://www.gnu.org/licenses/>
 */

#include "AT_RDD_ShellAveraged.h"

#include <gsl/gsl_sf_psi.h>

/* --------------------------------------------------- SHELL AVERAGE DOSE ---------------------------------------------------*/

 double   AT_RDD_Katz_LinearER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy){

  // Dav(r1,r2) = 2 * coeff * ( log(r2/r1) - (r2 - r1)/rmax ) / ((r2/rmax)^2 - (r1/rmax)^2)
  return 2.0 * Katz_point_coeff_Gy * ( log(r2_m/r1_m) - (r2_m - r1_m)/max_electron_range_m ) / (gsl_pow_2(r2_m/max_electron_range_m) - gsl_pow_2(r1_m/max_electron_range_m));
}

double   beta0(double x, double a, double iteration_threshold) {
    if (x < 0.99) {
        // This iterative algorithms computes the value of the beta function very accurately
        // but converges very slowly when the value of x is very close to 1.
        double ret = 0.0;
        double xi = 1.0;
        long i = 0;
        while (xi > iteration_threshold) {
            ret += xi / (a + i);
            xi *= x;
            i += 1;
        }
        ret *= pow(x, a);
        return ret;
    }
    else {
        // Series expansion taken from https://www.wolframalpha.com/input/?i=int+-1%2Fx+%281+-+x%29%5E%28a-1%29+dx
        // This can be obtained by Taylor expansion of the calculated integral + log(x) at x = 0.
        // Adding more terms didn't improve the accuracy.
        const double psi_a = gsl_sf_psi(a);
        const double psi_a_1 = psi_a + 1/a;
        const double psi_a_2 = psi_a_1 + 1/(a+1);
        const double psi_a_3 = psi_a_2 + 1/(a+2);
        const double a2 = a*a;
        return -psi_a - log(1.0-x) - M_EULER +
               a*(1.0-x)*(psi_a - psi_a_1 + 1.0) -
               0.25 * gsl_pow_2(1-x) * (
                       a * (a - 4*a*psi_a_1 + 2*a*psi_a_2 + 2 * (a-1)*psi_a + 2*psi_a_2 - 3)) +
               (1.0/18.0)*a*gsl_pow_3(1-x) * (
                       a2 + 9*a2 * psi_a_2 - 3*a2*psi_a_3 + 3*(a2 - 3*a + 2)*psi_a - 6*a -
                       9*(a-1)*a*psi_a_1 + 9*a*psi_a_2 - 9*a*psi_a_3 - 6*psi_a_3 + 11);
    }
}

double   AT_RDD_Katz_PowerLawER_DaverageKernel(  const double x1,
    const double x2,
    const double alpha){

  // this threshold was selected to get the overall error < 1e-10;
  const double iteration_threshold = 1e-13;
  const double ch1 = beta0(1-x1, 1+1/alpha, iteration_threshold);
  const double ch2 = beta0(1-x2, 1+1/alpha, iteration_threshold);
  return 2.0 * (ch1 - ch2) / (gsl_pow_2(x2) - gsl_pow_2(x1));
}


double   AT_RDD_Katz_PowerLawER_DaverageKernel_approx(  const double x1,
    const double x2,
    const double alpha){

   // F1 = x1 / alpha^2 ( (x1 / 4alpha) * (1/alpha - 1) - 1 ) + log(x1)
   // F2 = x2 / alpha^2 ( (x2 / 4alpha) * (1/alpha - 1) - 1 ) + log(x2)
   // kernel =  2/(x2^2 - x1^2) * (F2 - F1)
   const double F1 = (x1 / gsl_pow_2(alpha)) * ( (x1 / (4.0 * alpha)) * (1.0/alpha - 1.0) - 1.0 ) + log(x1);
   const double F2 = (x2 / gsl_pow_2(alpha)) * ( (x2 / (4.0 * alpha)) * (1.0/alpha - 1.0) - 1.0 ) + log(x2);
   return 2.0 * (F2 - F1) / (gsl_pow_2(x2) - gsl_pow_2(x1));
}


double   AT_RDD_Katz_PowerLawER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double alpha,
    const double Katz_point_coeff_Gy){

  //Dav(r1,r2) = coeff * kernel_av( x1, x2 )
  const double x1 = r1_m / max_electron_range_m;
  const double x2 = r2_m / max_electron_range_m;
  return Katz_point_coeff_Gy * (1.0/alpha) * AT_RDD_Katz_PowerLawER_DaverageKernel(x1, x2, alpha);
}


double AT_RDD_Cucinotta_Ddelta_average_integrand_m(  double r_m,
    void * params){

  AT_RDD_Cucinotta_Ddelta_average_integrand_m_parameters * RDD_parameters = (AT_RDD_Cucinotta_Ddelta_average_integrand_m_parameters*)(params);
  double max_electron_range_m  =  RDD_parameters->max_electron_range_m;
  double beta                  =  RDD_parameters->beta;
  double res                   =  1.0/r_m;
  res                         *=  AT_RDD_Cucinotta_f_shortRange( r_m, beta);
  res                         *=  AT_RDD_Cucinotta_f_longRange( r_m, max_electron_range_m);
  return res;
}


double   AT_RDD_Cucinotta_Ddelta_average_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double beta,
    const double Katz_point_coeff_Gy){

  // fL(r) = exp( -(r/(0.37rmax))^2 )
  // fS(r) = 1.0/( r0/r + 0.6 + 1.7 beta + 1.1 beta^2)           [here r0 = 10^(-9) [m]]
  // Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 Ddelta(r) 2 pi r dr
  //  thus: ...
  // Dav(r1,r2) = 2 * coeff/ ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 fS(r) * fL(r) * 1/r dr

  if( (r2_m > max_electron_range_m) || (r1_m > max_electron_range_m) || (r1_m > r2_m)){
#ifndef NDEBUG
    printf("wrong parameters given to AT_RDD_Cucinotta_Ddelta_average_Gy\n");
#endif
    return 0.0;
  }

  gsl_set_error_handler_off(); // TODO improve error handling

  double delta_average_integral;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (10000);
  gsl_function F;
  F.function = &AT_RDD_Cucinotta_Ddelta_average_integrand_m;
  AT_RDD_Cucinotta_Ddelta_average_integrand_m_parameters RDD_parameters;
  RDD_parameters.max_electron_range_m = max_electron_range_m;
  RDD_parameters.beta = beta;
  F.params = (void*)(&RDD_parameters);
  int status = gsl_integration_qags (&F, r1_m, r2_m, 1e-11, 1e-7, 10000, w1, &delta_average_integral, &error);
  //printf("integral = %g , error = %g, status = %d\n", delta_average_integral, error, status);
  if (status > 0){
#ifndef NDEBUG
    printf("integration error %d in AT_RDD_Cucinotta_Ddelta_average_Gy\n", status);
#endif
    delta_average_integral = -1.0;
  }
  gsl_integration_workspace_free (w1);

  return (2.0 * Katz_point_coeff_Gy / (gsl_pow_2(r2_m/max_electron_range_m) - gsl_pow_2(r1_m/max_electron_range_m))) * delta_average_integral ;
}


double   AT_RDD_Cucinotta_Dexc_average_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double beta,
    const double Katz_point_coeff_Gy){

  if( (r2_m > max_electron_range_m) || (r1_m > max_electron_range_m) || (r1_m > r2_m) || (r1_m < 0.0) ){
#ifndef NDEBUG
    printf("wrong parameters given to AT_RDD_Cucinotta_Dexc_average_Gy\n");
#endif
    return 0.0;
  }

  //  Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 Dexc(r) 2 pi r dr
  //  Dav(r1,r2) = 2 coeff/ ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 exp( - r / 2d ) * 1/r dr
  //
  // here:
  // integral = \int_r1^r2 exp( - r / 2d ) * 1/r * dr = Ei( -r2 / 2d ) - Ei( -r1 / 2d )
  // Ei is the exponential integral function
  const double d_m = beta * AT_RDD_Cucinotta_C_dm_wr;
  double FA = gsl_sf_expint_Ei( -0.5 * r1_m / d_m );
  double FB = 0;
  if( r2_m / d_m < 100.0 ){ // if x < -50 then Ei(x) < 1e-24, so we can assume than FB = 0
    FB = gsl_sf_expint_Ei( -0.5 * r2_m / d_m );
  }
  return (2.0 * Katz_point_coeff_Gy / (gsl_pow_2(r2_m/max_electron_range_m) - gsl_pow_2(r1_m/max_electron_range_m))) * (FB - FA) ;
}


double   AT_RDD_Cucinotta_Cnorm( const double r_min_m,
    const double max_electron_range_m,
    const double beta,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double Katz_point_coeff_Gy){

  const double LETfactor_Gy = LET_J_m * M_1_PI / (material_density_kg_m3 * (gsl_pow_2(max_electron_range_m) - gsl_pow_2(r_min_m)));
  const double Ddelta_average_Gy = AT_RDD_Cucinotta_Ddelta_average_Gy(r_min_m, max_electron_range_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
  const double Dexc_average_Gy = AT_RDD_Cucinotta_Dexc_average_Gy(r_min_m, max_electron_range_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
  if( (LETfactor_Gy > 0.0) && (Ddelta_average_Gy > 0.0) && (Dexc_average_Gy > 0.0)){
    return (LETfactor_Gy - Ddelta_average_Gy) / Dexc_average_Gy;
  } else {
#ifndef NDEBUG
    printf("problem in AT_RDD_Cucinotta_Cnorm\n");
#endif
    return 0.0;
  }
}


double   AT_RDD_Geiss_average_Gy(  const double r1_m,
    const double r2_m,
    const double a0_m,
    const double max_electron_range_m,
    const double norm_Gy){

  if( (r2_m > max_electron_range_m) || (r1_m > max_electron_range_m) || (r1_m > r2_m) || (r1_m < 0.0) ){
#ifndef NDEBUG
    printf("wrong parameters given to AT_RDD_Geiss_average_Gy\n");
#endif
    return 0.0;
  }

  // Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr

  if( r2_m <= a0_m ){
    return norm_Gy;
  } else if ( (r1_m < a0_m ) && ( r2_m > a0_m ) ) {
    return (norm_Gy / (gsl_pow_2(r2_m) - gsl_pow_2(r1_m)) ) * ( gsl_pow_2(a0_m) - gsl_pow_2(r1_m) + 2.0 * gsl_pow_2(a0_m) * log( r2_m / a0_m ));
  } else if ( r1_m >= a0_m ){
    return (norm_Gy / (gsl_pow_2(r2_m) - gsl_pow_2(r1_m)) ) * 2.0 * gsl_pow_2(a0_m) * log( r2_m / r1_m );
  }
  return 0.0;
}


/* --------------------------------------------------- dEdx IN OUTER SHELL ---------------------------------------------------*/


double   AT_RDD_Katz_LinearER_dEdx_J_m(  const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double Katz_point_coeff_Gy){

// dEdx = rho \int_a0^rmax r D(r) dr =
//      = rho * (pi rmax^2 - pi a0^2) * D_av(a0,rmax)
  return material_density_kg_m3 * M_PI * \
          (gsl_pow_2(max_electron_range_m) - gsl_pow_2(a0_m)) * \
          AT_RDD_Katz_LinearER_Daverage_Gy(a0_m, max_electron_range_m, max_electron_range_m, Katz_point_coeff_Gy);
}


double   AT_RDD_Katz_PowerLawER_dEdx_J_m(  const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double alpha,
    const double Katz_point_coeff_Gy){

  // dEdx = rho \int_a0^rmax r D(r) dr =
  //      = rho * (pi rmax^2 - pi a0^2) D_av(a0,rmax)
  return material_density_kg_m3 * M_PI * \
  (gsl_pow_2(max_electron_range_m) - gsl_pow_2(a0_m)) * \
  AT_RDD_Katz_PowerLawER_Daverage_Gy(a0_m, max_electron_range_m, max_electron_range_m, alpha, Katz_point_coeff_Gy);
}


/* --------------------------------------------------- SITE RDD ---------------------------------------------------*/


double   AT_RDD_Katz_LinearER_DSite_Gy( const double r_m,
    const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy){

  //Dsite(r) = 1 / (rho pi a0^2) * (LET - dEdx)   for r < a0
  //Dsite(r) = D(r)                               for r >= a0
  if( r_m < a0_m ){
    return M_1_PI * (LET_J_m - dEdx_J_m)/ (material_density_kg_m3 * gsl_pow_2(a0_m));
  } else {
    return AT_RDD_Katz_LinearER_Dpoint_Gy(r_m, max_electron_range_m, Katz_point_coeff_Gy);
  }

}


double   AT_RDD_Katz_PowerLawER_DSite_Gy( const double r_m,
    const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double alpha,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy){
  //Dsite(r) = 1 / (rho pi a0^2) * (LET - dEdx)   for r < a0
  //Dsite(r) = D(r)                               for r >= a0
  if( r_m <= a0_m ){
    return M_1_PI * (LET_J_m - dEdx_J_m)/ (material_density_kg_m3 * gsl_pow_2(a0_m));
  } else {
    return AT_RDD_Katz_PowerLawER_Dpoint_Gy(r_m, alpha, max_electron_range_m, Katz_point_coeff_Gy);
  }
}


 double  AT_RDD_KatzSite_Gy( const double r_m,
    const double r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy){

  if (r_m >= r_min_m && r_m <= max_electron_range_m){
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){
      return AT_RDD_Katz_PowerLawER_DSite_Gy(r_m, a0_m, max_electron_range_m, density_kg_m3, alpha, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
    } else if (er_model == ER_ButtsKatz){
      return AT_RDD_Katz_LinearER_DSite_Gy(r_m, a0_m, max_electron_range_m, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
    }
  }
  return 0.;
}


double  AT_inverse_RDD_KatzSite_m( const double D_Gy,
    const double r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double d_max_Gy,
    const double Katz_point_coeff_Gy){

  const double D_crit_Gy = AT_RDD_KatzPoint_Gy(a0_m, r_min_m, max_electron_range_m, er_model, alpha, Katz_point_coeff_Gy);
  if( D_Gy > d_max_Gy )
    return 0.0;
  if( D_Gy <= D_crit_Gy )
    return AT_inverse_RDD_KatzPoint_m(D_Gy, r_min_m, max_electron_range_m, er_model, alpha, Katz_point_coeff_Gy);
  else
    return a0_m;
}

