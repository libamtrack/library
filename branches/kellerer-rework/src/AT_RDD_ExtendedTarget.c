/**
 * @brief ...
 */

/*
 *    AT_RDD_ExtendedTarget.c
 *    ===========================
 *
 *    Created on: 01.03.2010
 *    Creator: greilich
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

#include "AT_RDD_ExtendedTarget.h"


double geometryFunctionPhi(         const double r_m,
    const double a0_m,
    const double t_m)
{
  /* Phi(r,a0,t) = 2 arctan( sqrt( (a0^2 - (r-t)^2) / ((t + r)^2  - a0^2) ) )   for t > | r - a0 |  */
  /* Phi(r,a0,t) = pi                                                          for t <= | r - a0 |  */
  double res = 0.;
  double factor = 0.;
  gsl_complex carg, cres;
  if (t_m <= fabs (r_m - a0_m)){
    if (r_m >= a0_m)
      res = 0.0;
    else
      res = M_PI;
  } else {
    factor  = gsl_pow_2 (a0_m) - gsl_pow_2 (r_m - t_m);
    factor /= gsl_pow_2 (t_m + r_m) - gsl_pow_2 (a0_m);
    GSL_SET_COMPLEX (&carg, sqrt (factor), 0.);
    cres = gsl_complex_arctan (carg);
    res = 2.0 * GSL_REAL (cres);
  }
  return res;
}

/* --------------------------------------------------- KATZ EXT TARGET RDD ---------------------------------------------------*/

double AT_RDD_ExtendedTarget_KatzPoint_integrand_Gy(
    double t_m,
    void* params){

  AT_RDD_ExtendedTarget_KatzPoint_parameters * RDD_parameters = (AT_RDD_ExtendedTarget_KatzPoint_parameters*)params;
  const double  r_m                   =  RDD_parameters->r_m;
  const double  a0_m                  =  RDD_parameters->a0_m;

  const double  D_Gy                  =  AT_RDD_KatzPoint_Gy( t_m,
      RDD_parameters->KatzPoint_r_min_m,
      RDD_parameters->max_electron_range_m,
      RDD_parameters->er_model,
      RDD_parameters->alpha,
      RDD_parameters->Katz_point_coeff_Gy);

  return 2.0 * t_m * D_Gy * geometryFunctionPhi(r_m, a0_m, t_m);
}


double AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration(
    const double  r_m,
    const double  a0_m,
    const long    er_model,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  alpha,
    const double  Katz_point_coeff_Gy){

  double low_lim_m = 0.0;
  if( r_m > a0_m ){
    low_lim_m = GSL_MAX(r_m - a0_m,0.0);
  }
  gsl_set_error_handler_off();

  double ext_integral_Gy;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &AT_RDD_ExtendedTarget_KatzPoint_integrand_Gy;

  AT_RDD_ExtendedTarget_KatzPoint_parameters RDD_parameters;

  RDD_parameters.r_m                  =  r_m;
  RDD_parameters.a0_m                 =  a0_m;
  RDD_parameters.er_model             =  er_model;
  RDD_parameters.max_electron_range_m =  max_electron_range_m;
  RDD_parameters.alpha                =  alpha;
  RDD_parameters.KatzPoint_r_min_m    =  KatzPoint_r_min_m;
  RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;

  F.params = (void*)(&RDD_parameters);
  int status = gsl_integration_qags (&F, low_lim_m, r_m + a0_m, 0, 1e-5, 1000, w1, &ext_integral_Gy, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
    printf("Error in AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration\n");
    ext_integral_Gy = -1.0;
  }
  gsl_integration_workspace_free (w1);

  ext_integral_Gy *= M_1_PI / gsl_pow_2(a0_m);

  return ext_integral_Gy;
}


double AT_RDD_ExtendedTarget_KatzPoint_Gy(
    const double  r_m,
    const double  a0_m,
    const long    er_model,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  alpha,
    const double  Katz_plateau_Gy,
    const double  Katz_point_coeff_Gy){

  double D_Gy   =  0.0;

  const double r_max_m = GSL_MIN(a0_m, max_electron_range_m);

  /* plateau region */
  if( (r_m <=  0.01 * a0_m) && (r_m >= 0.0)){
    D_Gy  =  Katz_plateau_Gy;
    if( max_electron_range_m < a0_m ){
      D_Gy *= (gsl_pow_2(r_max_m / a0_m));
    }
  }

  /* intermediate */
  if( (r_m <  100.0 * a0_m) && (r_m >  0.01 * a0_m) ){
    D_Gy  =  AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration(r_m, a0_m, er_model, KatzPoint_r_min_m, max_electron_range_m, alpha, Katz_point_coeff_Gy);
  }

  /* far 1/r^2 region */
  if( (r_m >=  100.0 * a0_m) && (r_m <= max_electron_range_m) ){
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
      D_Gy  =  AT_RDD_Katz_PowerLawER_Dpoint_Gy( r_m, alpha, max_electron_range_m, Katz_point_coeff_Gy );
    } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
      D_Gy  =  AT_RDD_Katz_LinearER_Dpoint_Gy( r_m, max_electron_range_m, Katz_point_coeff_Gy );
    }
  }

  return D_Gy;
}


double AT_inverse_RDD_ExtendedTarget_KatzPoint_solver_function_Gy( const double r_m , void * params ){
  AT_inverse_RDD_ExtendedTarget_KatzPoint_parameters* RDD_parameters = (AT_inverse_RDD_ExtendedTarget_KatzPoint_parameters*)(params);

  return AT_RDD_ExtendedTarget_KatzPoint_Gy(r_m,
      RDD_parameters->a0_m,
      RDD_parameters->er_model,
      RDD_parameters->Katz_point_r_min_m,
      RDD_parameters->max_electron_range_m,
      RDD_parameters->alpha,
      RDD_parameters->Katz_plateau_Gy,
      RDD_parameters->Katz_point_coeff_Gy) - RDD_parameters->D_Gy;
}


double  AT_inverse_RDD_ExtendedTarget_KatzPoint_m( const double D_Gy,
    const double KatzPoint_r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double Katz_plateau_Gy,
    const double Katz_point_coeff_Gy){

  const double  solver_accuracy  =  1e-13;

  AT_inverse_RDD_ExtendedTarget_KatzPoint_parameters RDD_parameters;

  RDD_parameters.D_Gy                 =  D_Gy;
  RDD_parameters.a0_m                 =  a0_m;
  RDD_parameters.er_model             =  er_model;
  RDD_parameters.max_electron_range_m =  max_electron_range_m;
  RDD_parameters.alpha                =  alpha;
  RDD_parameters.Katz_plateau_Gy      =  Katz_plateau_Gy;
  RDD_parameters.Katz_point_r_min_m   =  KatzPoint_r_min_m;
  RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;

  double r_m =  zriddr(AT_inverse_RDD_ExtendedTarget_KatzPoint_solver_function_Gy,
          (void*)(&RDD_parameters),
          0.0,
          max_electron_range_m + a0_m,
          solver_accuracy);

  return r_m;
}

/* --------------------------------------------------- CUCINOTTA EXT TARGET RDD ---------------------------------------------------*/


double AT_RDD_ExtendedTarget_CucinottaPoint_integrand_Gy(
    double t_m,
    void* params){

  AT_RDD_ExtendedTarget_CucinottPoint_parameters * RDD_parameters = (AT_RDD_ExtendedTarget_CucinottPoint_parameters*)params;
  const double  r_m                   =  RDD_parameters->r_m;
  const double  a0_m                  =  RDD_parameters->a0_m;

  const double  D_Gy                  =  AT_RDD_CucinottaPoint_Gy( t_m,
      RDD_parameters->KatzPoint_r_min_m,
      RDD_parameters->max_electron_range_m,
      RDD_parameters->beta,
      RDD_parameters->C_norm,
      RDD_parameters->Katz_point_coeff_Gy);

  return 2.0 * t_m * D_Gy * geometryFunctionPhi(r_m, a0_m, t_m);
}


double AT_RDD_ExtendedTarget_CucinottaPoint_Gy_by_integration(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  Katz_point_coeff_Gy,
    const double  C_norm){

  double low_lim_m = 0.0;
  if( r_m > a0_m ){
    low_lim_m = GSL_MAX(r_m - a0_m,0.0);
  }
  gsl_set_error_handler_off();

  double ext_integral_Gy;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &AT_RDD_ExtendedTarget_CucinottaPoint_integrand_Gy;

  AT_RDD_ExtendedTarget_CucinottPoint_parameters RDD_parameters;

  RDD_parameters.r_m                  =  r_m;
  RDD_parameters.a0_m                 =  a0_m;
  RDD_parameters.KatzPoint_r_min_m    =  KatzPoint_r_min_m;
  RDD_parameters.max_electron_range_m =  max_electron_range_m;
  RDD_parameters.beta                 =  beta;
  RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;
  RDD_parameters.C_norm               =  C_norm;

  F.params = (void*)(&RDD_parameters);
  int status = gsl_integration_qags (&F, low_lim_m, r_m + a0_m, 0, 1e-5, 1000, w1, &ext_integral_Gy, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
    printf("Error in AT_RDD_ExtendedTarget_CucinottaPoint_Gy_by_integration\n");
    ext_integral_Gy = -1.0;
  }
  gsl_integration_workspace_free (w1);

  ext_integral_Gy *= M_1_PI / gsl_pow_2(a0_m);

  return ext_integral_Gy;
}


double AT_RDD_ExtendedTarget_CucinottaPoint_Gy(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  Katz_point_coeff_Gy,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy){

  double D_Gy   =  0.0;

  const double r_max_m = GSL_MIN(a0_m, max_electron_range_m);
  if( (r_m <=  0.01 * a0_m) && (r_m >= 0.0)){
    D_Gy  =  Cucinotta_plateau_Gy;
    if( max_electron_range_m < a0_m ){
      D_Gy *= (gsl_pow_2(r_max_m / a0_m));
    }
  }

  /* intermediate */
  if( (r_m <  100.0 * a0_m) && (r_m >  0.01 * a0_m) ){
    D_Gy  =  AT_RDD_ExtendedTarget_CucinottaPoint_Gy_by_integration( r_m, a0_m, KatzPoint_r_min_m, max_electron_range_m, beta, Katz_point_coeff_Gy, C_norm);
  }

  /* far 1/r^2 region */
  if( (r_m >=  100.0 * a0_m) && (r_m <= max_electron_range_m) ){
    D_Gy = AT_RDD_CucinottaPoint_Gy( r_m, KatzPoint_r_min_m, max_electron_range_m, beta, C_norm, Katz_point_coeff_Gy);
  }

  return D_Gy;
}


double AT_inverse_RDD_ExtendedTarget_CucinottaPoint_solver_function_Gy( const double r_m , void * params ){
  AT_inverse_RDD_ExtendedTarget_CucinottaPoint_parameters* RDD_parameters = (AT_inverse_RDD_ExtendedTarget_CucinottaPoint_parameters*)(params);

  return AT_RDD_ExtendedTarget_CucinottaPoint_Gy(r_m,
      RDD_parameters->a0_m,
      RDD_parameters->KatzPoint_r_min_m,
      RDD_parameters->max_electron_range_m,
      RDD_parameters->beta,
      RDD_parameters->Katz_point_coeff_Gy,
      RDD_parameters->C_norm,
      RDD_parameters->Cucinotta_plateau_Gy) - RDD_parameters->D_Gy;
}


double  AT_inverse_RDD_ExtendedTarget_CucinottaPoint_m( const double D_Gy,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  Katz_point_coeff_Gy,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy){

  const double  solver_accuracy  =  1e-13;

  AT_inverse_RDD_ExtendedTarget_CucinottaPoint_parameters RDD_parameters;

  RDD_parameters.D_Gy                 =  D_Gy;
  RDD_parameters.a0_m                 =  a0_m;
  RDD_parameters.KatzPoint_r_min_m    =  KatzPoint_r_min_m;
  RDD_parameters.max_electron_range_m =  max_electron_range_m;
  RDD_parameters.beta                 =  beta;
  RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;
  RDD_parameters.C_norm               =  C_norm;
  RDD_parameters.Cucinotta_plateau_Gy =  Cucinotta_plateau_Gy;

  double r_m =  zriddr(AT_inverse_RDD_ExtendedTarget_CucinottaPoint_solver_function_Gy,
          (void*)(&RDD_parameters),
          0.0,
          max_electron_range_m+a0_m,
          solver_accuracy);

  return r_m;
}

