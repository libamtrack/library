/**
 * @file
 * @brief Radial Dose Distribution models
 */

/*
 *    AT_RDD_Simple.c
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

#include "AT_RDD_Simple.h"

/* --------------------------------------------------- TEST RDD ---------------------------------------------------*/

inline double  AT_RDD_Test_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double norm_constant_Gy){

  if( r_m >= r_min_m && r_m <= r_max_m ){
    return norm_constant_Gy;
  } else {
    return 0.0;
  }
}


inline double  AT_inverse_RDD_Test_m( const double D_Gy,
    const double r_max_m){
  if( D_Gy > 0)
    return r_max_m;
  return 0.0;
}

/* --------------------------------------------------- GEISS RDD ---------------------------------------------------*/


inline double  AT_RDD_Geiss_Gy( const double r_m,
    const double r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const double norm_constant_Gy){

  if ( (r_min_m <= r_m) && (r_m <= a0_m)){
    return norm_constant_Gy;
  }
  if ((a0_m < r_m) && (r_m <= max_electron_range_m)){
    return norm_constant_Gy * gsl_pow_2(a0_m / r_m);
  }
  return 0.;
}


inline double  AT_inverse_RDD_Geiss_m( const double D_Gy,
    const double d_min_Gy,
    const double d_max_Gy,
    const double a0_m,
    const double norm_constant_Gy){
  if ((d_min_Gy <= D_Gy) && (D_Gy <= d_max_Gy)){
    return a0_m * sqrt(norm_constant_Gy / D_Gy);
  }
  return 0.0;
}

/* --------------------------------------------------- KATZ RDD ---------------------------------------------------*/


inline double AT_RDD_Katz_coeff_Gy( const double Z_eff,
    const double beta,
    const double density_kg_m3,
    const double electron_density_m3,
    const double max_electron_range_m){

  const double C_J_m = electron_density_m3 * AT_Katz_C1_J_m2;
  // (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
  return C_J_m * 0.5 * M_1_PI * gsl_pow_2(Z_eff/beta) / (density_kg_m3 * gsl_pow_2(max_electron_range_m));
}


inline double AT_RDD_Katz_coeff_Gy_general(     const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    er_model){

  double beta   =  AT_beta_from_E_single( E_MeV_u );

  // Get densities
  double   density_g_cm3, electron_density_m3;
  AT_get_material_data(
      material_no,
      &density_g_cm3,
      &electron_density_m3,
      NULL, NULL, NULL, NULL, NULL, NULL);
  double density_kg_m3      =  density_g_cm3 * 1000.0;

  const double max_electron_range_m  =  AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);

  const long   Z                     =  AT_Z_from_particle_no_single(particle_no);
  const double Z_eff                 =  AT_effective_charge_from_beta_single(beta, Z);
  const double Katz_point_coeff_Gy   =  AT_RDD_Katz_coeff_Gy(Z_eff, beta, density_kg_m3, electron_density_m3, max_electron_range_m);

  return Katz_point_coeff_Gy;
}


inline double   AT_RDD_Katz_LinearER_Dpoint_Gy(        const double r_m,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy){

  // D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/rmax^2 * rmax/r * (rmax/r - 1.)
  return Katz_point_coeff_Gy * max_electron_range_m/r_m * (max_electron_range_m/r_m - 1.0);
}


inline double AT_RDD_Katz_PowerLawER_Dpoint_Gy(const double r_m,
    const double alpha,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy){

  // D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/r^2 * 1/alpha * (1 - r/rmax)^(1/alpha)
  if( r_m > max_electron_range_m )
    return 0.0;
  else
    return Katz_point_coeff_Gy * (1.0/alpha) * gsl_pow_2(max_electron_range_m/r_m) * pow(1.0 - r_m/max_electron_range_m, 1.0 / (alpha));
}


inline double  AT_RDD_KatzPoint_Gy( const double r_m,
    const double r_min_m,
    const double max_electron_range_m,
    const long   er_model,
    const double alpha,
    const double Katz_point_coeff_Gy){

  if (r_m >= r_min_m && r_m <= max_electron_range_m){
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // use power law form of RDD for power law ER models (new Katz style)
      return AT_RDD_Katz_PowerLawER_Dpoint_Gy(r_m, alpha, max_electron_range_m, Katz_point_coeff_Gy);
    } else if (er_model == ER_ButtsKatz){ // linear law form of RDD for other ER models (old Katz style)
      return AT_RDD_Katz_LinearER_Dpoint_Gy(r_m, max_electron_range_m, Katz_point_coeff_Gy);
    }
  }
  return 0.;
}


double  AT_inverse_RDD_KatzPoint_LinearER_m( const double D_Gy,
    const double r_min_m,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy){
  double x = 0.5 * (1.0 + sqrt(1.0 + 4.0 * D_Gy / Katz_point_coeff_Gy));
  return max_electron_range_m / x;
}


double AT_inverse_RDD_KatzPoint_PowerLawER_solver_function_Gy( const double r_m , void * params ){
  AT_inverse_RDD_KatzPoint_PowerLawER_parameters* RDD_parameters = (AT_inverse_RDD_KatzPoint_PowerLawER_parameters*)(params);

  if( r_m < (RDD_parameters->r_min_m)){
    return -(RDD_parameters->D_Gy);
  } else {
    return AT_RDD_Katz_PowerLawER_Dpoint_Gy(r_m,
        RDD_parameters->alpha,
        RDD_parameters->max_electron_range_m,
        RDD_parameters->Katz_point_coeff_Gy) - (RDD_parameters->D_Gy);
  }
}


double  AT_inverse_RDD_KatzPoint_m( const double D_Gy,
    const double r_min_m,
    const double max_electron_range_m,
    const long   er_model,
    const double alpha,
    const double Katz_point_coeff_Gy){

  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
    double  solver_accuracy  =  1e-13;

    AT_inverse_RDD_KatzPoint_PowerLawER_parameters RDD_parameters;

    RDD_parameters.D_Gy                 =  D_Gy;
    RDD_parameters.r_min_m              =  r_min_m;
    RDD_parameters.max_electron_range_m =  max_electron_range_m;
    RDD_parameters.alpha                =  alpha;
    RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;
    double r_m =  zriddr(AT_inverse_RDD_KatzPoint_PowerLawER_solver_function_Gy,
            (void*)(&RDD_parameters),
            r_min_m,
            max_electron_range_m,
            solver_accuracy);

    return r_m;

  } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
    return AT_inverse_RDD_KatzPoint_LinearER_m(D_Gy, r_min_m, max_electron_range_m, Katz_point_coeff_Gy);
  }

  return 0.0;
}


/* --------------------------------------------------- CUCINOTTA RDD ---------------------------------------------------*/


inline double   AT_RDD_Cucinotta_f_shortRange( const double r_m,
    const double beta){

  // fS(r) = 1.0/( r0/r + 0.6 + 1.7 beta + 1.1 beta^2)           [here r0 = 10^(-9) [m]]
  return 1.0/(1e-9/r_m + 0.6 + 1.7 * beta + 1.1 * gsl_pow_2(beta));
}


inline double   AT_RDD_Cucinotta_f_longRange( const double r_m,
    const double max_electron_range_m){

  //fL(r) = exp( -(r/(0.37rmax))^2 )
  return exp( - gsl_pow_2( r_m / (0.37 * max_electron_range_m ) ) );
}


inline double AT_RDD_Cucinotta_Ddelta_Gy( const double r_m,
    const double max_electron_range_m,
    const double beta,
    const double Katz_point_coeff_Gy){

  // Ddelta(r) = coeff * fS(r) * fL(r) * rmax^2/r^2
  return Katz_point_coeff_Gy * AT_RDD_Cucinotta_f_longRange(r_m, max_electron_range_m) * AT_RDD_Cucinotta_f_shortRange(r_m, beta) * gsl_pow_2(max_electron_range_m/r_m);
}


double AT_RDD_Cucinotta_Ddelta_average_integrand_m(  double r_m,
    void * params){

  double max_electron_range_m  =  ((double*)params)[0];
  double beta                  =  ((double*)params)[1];
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
    printf("wrong parameters given to AT_RDD_Cucinotta_Ddelta_average_Gy\n");
    return 0.0;
  }

  gsl_set_error_handler_off(); // TODO improve error handling

  double delta_average_integral;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (10000);
  gsl_function F;
  F.function = &AT_RDD_Cucinotta_Ddelta_average_integrand_m;
  double params[] = {max_electron_range_m, beta};
  F.params = params;
  int status = gsl_integration_qags (&F, r1_m, r2_m, 1e-11, 1e-7, 10000, w1, &delta_average_integral, &error);
  //printf("integral = %g , error = %g, status = %d\n", delta_average_integral, error, status);
  if (status > 0){
    printf("integration error %d in AT_RDD_Cucinotta_Ddelta_average_Gy\n", status);
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

  if( (r2_m > max_electron_range_m) || (r1_m > max_electron_range_m) || (r1_m > r2_m) || (r1_m <= 0.0) ){
    printf("wrong parameters given to AT_RDD_Cucinotta_Dexc_average_Gy\n");
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
    printf("problem in AT_RDD_Cucinotta_Cnorm\n");
    return 0.0;
  }
}


double   AT_RDD_Cucinotta_Dexc_Gy( const double r_m,
    const double max_electron_range_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy){

  // Dexc(r) = C exp( - r / 2d ) / r^2            [where d = (beta/2) * (hbar * c / wr) and wr = 13eV ]
  // Dexc(r) = Cnorm * coeff * exp( - r / 2d ) * (rmax/r)^2
  const double d_m = beta * AT_RDD_Cucinotta_C_dm_wr;
  return C_norm * Katz_point_coeff_Gy * exp( - 0.5 * r_m / d_m ) * gsl_pow_2(max_electron_range_m/r_m);
}


inline double   AT_RDD_CucinottaPoint_Gy( const double r_m,
    const double r_min_m,
    const double max_electron_range_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy){

  if( (r_m >= r_min_m) && (r_m <= max_electron_range_m) ){
    // D(r)    = Dexc(r) + Ddelta(r)
    return AT_RDD_Cucinotta_Dexc_Gy(r_m, max_electron_range_m, beta, C_norm, Katz_point_coeff_Gy) + AT_RDD_Cucinotta_Ddelta_Gy(r_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
  } else {
    return 0.0;
  }
}


double AT_inverse_RDD_Cucinotta_solver_function_Gy( const double r_m , void * params ){
  AT_inverse_RDD_Cucinotta_parameters* RDD_parameters = (AT_inverse_RDD_Cucinotta_parameters*)(params);

  if( r_m < (RDD_parameters->r_min_m)){
    return -(RDD_parameters->D_Gy);
  } else {
    return AT_RDD_CucinottaPoint_Gy(r_m,
        RDD_parameters->r_min_m,
        RDD_parameters->max_electron_range_m,
        RDD_parameters->beta,
        RDD_parameters->C_norm,
        RDD_parameters->Katz_point_coeff_Gy) - (RDD_parameters->D_Gy);
  }
}


double  AT_inverse_RDD_Cucinotta_m( const double D_Gy,
    const double r_min_m,
    const double max_electron_range_m,
    const long   er_model,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy){

  if( er_model == ER_Tabata ){
    double  solver_accuracy  =  1e-13;

    AT_inverse_RDD_Cucinotta_parameters RDD_parameters;

    RDD_parameters.D_Gy                 =  D_Gy;
    RDD_parameters.r_min_m              =  r_min_m;
    RDD_parameters.max_electron_range_m =  max_electron_range_m;
    RDD_parameters.beta                 =  beta;
    RDD_parameters.C_norm               =  C_norm;
    RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;

    double r_m =  zriddr(AT_inverse_RDD_Cucinotta_solver_function_Gy,
            (void*)(&RDD_parameters),
            r_min_m,
            max_electron_range_m,
            solver_accuracy);

    return r_m;
  }

  return 0.0;
}

