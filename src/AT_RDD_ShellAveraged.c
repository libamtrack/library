/**
 * @file
 * @brief Radial Dose Distribution models
 */

/*
 *    AT_RDD_ShellAveraged.c
 *    ========
 *
 *    Created on: 05.04.2010
 *    Author: kongruencja
 *
 *    Copyright 2006, 2009 Steffen Greilich / the libamtrack team
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

/* --------------------------------------------------- SHELL AVERAGE DOSE ---------------------------------------------------*/

inline double   AT_RDD_Katz_LinearER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy){

  // Dav(r1,r2) = 2 * coeff * ( log(r2/r1) - (r2 - r1)/rmax ) / ((r2/rmax)^2 - (r1/rmax)^2)
  return 2.0 * Katz_point_coeff_Gy * ( log(r2_m/r1_m) - (r2_m - r1_m)/max_electron_range_m ) / (gsl_pow_2(r2_m/max_electron_range_m) - gsl_pow_2(r1_m/max_electron_range_m));
}


double   AT_RDD_Katz_PowerLawER_DaverageKernel(  const double x1,
    const double x2,
    const double alpha){

  // C1 = (1-x1)^(1/alpha) ((x1-1)/x1)^(-1/alpha)
  // HGF1 =  _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x1)
  // C2 = (1-x2)^(1/alpha) ((x2-1)/x2)^(-1/alpha)
  // HGF2 =  _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x2)
  // kernel =  2/(x2^2 - x1^2) * (C2* HGF2 - C1 * HGF1)
  const double alpha_inv = 1.0/alpha;
  const double C1 = pow( 1.0-x1, alpha_inv ) * pow( (x1-1.0)/x1, -alpha_inv );
  const double HGF1 = gsl_sf_hyperg_2F1( -alpha_inv, -alpha_inv, 1.0-alpha_inv, 1.0/x1 );
  const double C2 = pow( 1.0-x2, alpha_inv ) * pow( (x2-1.0)/x2, -alpha_inv );
  const double HGF2 = gsl_sf_hyperg_2F1( -alpha_inv, -alpha_inv, 1.0-alpha_inv, 1.0/x2 );
  return 2.0 * (C2 * HGF2 - C1 * HGF1) / (gsl_pow_2(x2) - gsl_pow_2(x1));
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
  return Katz_point_coeff_Gy * AT_RDD_Katz_PowerLawER_DaverageKernel_approx(x1, x2, alpha);
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


inline double  AT_RDD_KatzSite_Gy( const double r_m,
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

