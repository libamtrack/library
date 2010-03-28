/**
* @file
* @brief ...
*/

/*
 *    AT_RDD_ExtendedTarget.c
 *    ===========================
 *
 *    Created on: 01.03.2010
 *    Author: greilich
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

#include "AT_RDD_ExtendedTarget.h"

///////////////////////////////////////////// SITE /////////////////////////////////////////////
void AT_RDD_Site_Gy( const long  n,
    const float   r_m[],
    const double  a0_m,
    /* radiation field parameters */
    const double  E_MeV_u,
    const long    particle_no,
    /* detector parameters */
    const long    material_no,
    /* radial dose distribution model */
    const long    rdd_model,
    const float   rdd_parameter[],
    /* electron range model */
    const long    er_model,
    const float   er_parameter[],
    float         D_RDD_Gy[]){

  D_RDD_Gy[0] = 0.0;
}


///////////////////////////////////////////// EXTENDED TARGET /////////////////////////////////////////////

double          geometryFunctionPhi(         const double r_m,
    const double a0_m,
    const double t_m)
{
  // Phi(r,a0,t) = 2 arctan( sqrt( (a0^2 - (r-t)^2) / ((t + r)^2  - a0^2) ) )   for t > | r - a0 |
  // Phi(r,a0,t) = pi                                                          for t <= | r - a0 |
  double res = 0.;
  double factor = 0.;
  gsl_complex carg, cres;
  if (t_m <= fabs (r_m - a0_m)){
    if (r_m >= a0_m)
      res = 0.0;
    else
      res = M_PI;
  } else {
    factor = gsl_pow_2 (a0_m) - gsl_pow_2 (r_m - t_m);
    factor /= gsl_pow_2 (t_m + r_m) - gsl_pow_2 (a0_m);
    GSL_SET_COMPLEX (&carg, sqrt (factor), 0.);
    cres = gsl_complex_arctan (carg);
    res = 2.0 * GSL_REAL (cres);
  }
  return res;
}


double AT_RDD_ExtendedTarget_KatzPoint_integrand_Gy(
    double t_m,
    void* params){

  AT_RDD_ExtendedTarget_KatzPoint_parameters * RDD_parameters = (AT_RDD_ExtendedTarget_KatzPoint_parameters*)params;
  const double  r_m                   =  RDD_parameters->r_m;
  const double  a0_m                  =  RDD_parameters->a0_m;
  const long    er_model              =  RDD_parameters->er_model;
  const double  r_min_m               =  RDD_parameters->KatzPoint_r_min_m;
  const double  r_max_m               =  RDD_parameters->r_max_m;
  const double  alpha                 =  RDD_parameters->alpha;
  const double  Katz_point_coeff_Gy   =  RDD_parameters->Katz_point_coeff_Gy;


  const double  D_Gy                  = AT_RDD_KatzPoint_Gy( t_m, r_min_m, r_max_m, er_model, alpha, Katz_point_coeff_Gy);

  return 2.0 * t_m * D_Gy * geometryFunctionPhi(r_m, a0_m, t_m);
}


double AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration(
    const double  r_m,
    const double  a0_m,
    const long    er_model,
    const double  KatzPoint_r_min_m,
    const double  r_max_m,
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
  RDD_parameters.r_max_m              =  r_max_m;
  RDD_parameters.alpha                =  alpha;
  RDD_parameters.KatzPoint_r_min_m    =  KatzPoint_r_min_m;
  RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;

  F.params = (void*)(&RDD_parameters);
  int status = gsl_integration_qags (&F, low_lim_m, r_m+a0_m, 0, 1e-5, 1000, w1, &ext_integral_Gy, &error);
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

  // plateau region
  if( (r_m <=  0.01 * a0_m) && (r_m >= 0.0)){
    D_Gy  =  Katz_plateau_Gy;
    if( max_electron_range_m < a0_m ){
      D_Gy *= (gsl_pow_2(r_max_m / a0_m));
    }
  }

  // intermediate
  if( (r_m <  100.0 * a0_m) && (r_m >  0.01 * a0_m) ){
    D_Gy  =  AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration(r_m, a0_m, er_model, KatzPoint_r_min_m, max_electron_range_m, alpha, Katz_point_coeff_Gy);
  }

  // far 1/r^2 region
  if( (r_m >=  100.0 * a0_m) && (r_m <= max_electron_range_m) ){
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
      D_Gy  =  AT_RDD_Katz_PowerLawER_Dpoint_Gy( r_m, alpha, max_electron_range_m, Katz_point_coeff_Gy );
    } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
      D_Gy  =  AT_RDD_Katz_LinearER_Dpoint_Gy( r_m, max_electron_range_m, Katz_point_coeff_Gy );
    }
  }

  return D_Gy;
}


float AT_inverse_RDD_ExtendedTarget_KatzPoint_solver_function_Gy( const float r_m , void * params ){
  AT_inverse_RDD_ExtendedTarget_KatzPoint_parameters* RDD_parameters = (AT_inverse_RDD_ExtendedTarget_KatzPoint_parameters*)(params);

  const double  D_Gy                  =  RDD_parameters->D_Gy;
  const double  a0_m                  =  RDD_parameters->a0_m;
  const long    er_model              =  RDD_parameters->er_model;
  const double  r_max_m               =  RDD_parameters->r_max_m;
  const double  alpha                 =  RDD_parameters->alpha;
  const double  Katz_plateau_Gy       =  RDD_parameters->Katz_plateau_Gy;
  const double  Katz_point_r_min_m    =  RDD_parameters->Katz_point_r_min_m;
  const double  Katz_point_coeff_Gy   =  RDD_parameters->Katz_point_coeff_Gy;

  return (float)AT_RDD_ExtendedTarget_KatzPoint_Gy((double)r_m, a0_m, er_model, Katz_point_r_min_m, r_max_m, alpha, Katz_plateau_Gy, Katz_point_coeff_Gy) - (float)D_Gy;
}


double  AT_inverse_RDD_ExtendedTarget_KatzPoint_m( const double D_Gy,
    const double KatzPoint_r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double Katz_plateau_Gy,
    const double Katz_point_coeff_Gy){

  float  solver_accuracy  =  1e-13f;

  AT_inverse_RDD_ExtendedTarget_KatzPoint_parameters RDD_parameters;

  RDD_parameters.D_Gy                 =  D_Gy;
  RDD_parameters.a0_m                 =  a0_m;
  RDD_parameters.er_model             =  er_model;
  RDD_parameters.r_max_m              =  max_electron_range_m;
  RDD_parameters.alpha                =  alpha;
  RDD_parameters.Katz_plateau_Gy      =  Katz_plateau_Gy;
  RDD_parameters.Katz_point_r_min_m   =  KatzPoint_r_min_m;
  RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;

  double r_m =  zriddr(AT_inverse_RDD_ExtendedTarget_KatzPoint_solver_function_Gy,
          (void*)(&RDD_parameters),
          0.0f,
          (float)max_electron_range_m+(float)a0_m,
          solver_accuracy);

  return r_m;
}

void AT_RDD_ExtendedTarget_Gy( const long  n,
    const float  r_m[],
    const double a0_m,
    /* radiation field parameters */
    const double E_MeV_u,
    const long   particle_no,
    /* detector parameters */
    const long   material_no,
    /* radial dose distribution model */
    const long   rdd_model,
    const float  rdd_parameter[],
    /* electron range model */
    const long   er_model,
    const float  er_parameter[],
    float        D_RDD_Gy[]){

  // Get material data
  double   density_g_cm3, density_kg_m3, electron_density_m3;
  AT_get_material_data(
      material_no,
      &density_g_cm3,
      &electron_density_m3,
      NULL, NULL, NULL, NULL, NULL, NULL);
  density_kg_m3      =  density_g_cm3 * 1000.0;

  // Get beta, Z and Zeff
  double beta    =  AT_beta_from_E_single(E_MeV_u);
  long   Z       =  AT_Z_from_particle_no_single(particle_no);
  double Z_eff   =  AT_effective_charge_from_beta_single(beta, Z);

  // Get the maximum electron range
  const double  max_electron_range_m      =  AT_max_electron_range_m(E_MeV_u, (int)material_no, (int)er_model);

  const double C_J_m                      =  AT_RDD_Katz_C_J_m(electron_density_m3);

  double Katz_point_coeff_Gy =  AT_RDD_Katz_coeff_Gy(C_J_m, Z_eff, beta, density_kg_m3, max_electron_range_m);
  double r_max_m = GSL_MIN(a0_m, max_electron_range_m);

  double alpha           =  0.0;
  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
    alpha               =  1.667;
    double wmax_MeV     =  AT_max_E_transfer_MeV_single(E_MeV_u);
    if(wmax_MeV <= 1e-3){
      alpha             =  1.079;
    }
  }

  /********************************************************
  ****************** LOOP OVER DISTANCES *****************
  *******************************************************/
  long     i;

  double r_min_m = (double)rdd_parameter[0];

  if( rdd_model == RDD_KatzPoint){ // TODO write comments

    double KatzPoint_r_min_m = (double)rdd_parameter[0];

    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD

      double Katz_PowerLawER_plateau_Gy = AT_RDD_Katz_PowerLawER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, alpha, Katz_point_coeff_Gy );

      for( i = 0 ; i < n ; i++){
          D_RDD_Gy[i]  =  (float)AT_RDD_ExtendedTarget_KatzPoint_Gy((double)r_m[i], a0_m, er_model, KatzPoint_r_min_m, max_electron_range_m, alpha, Katz_PowerLawER_plateau_Gy, Katz_point_coeff_Gy);
      }
      // end if ER_Waligorski, ER_Edmund
    } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD

      double Katz_LinearLawER_plateau_Gy = AT_RDD_Katz_LinearER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, Katz_point_coeff_Gy );

      for( i = 0 ; i < n ; i++){
            D_RDD_Gy[i]  =  (float)AT_RDD_ExtendedTarget_KatzPoint_Gy((double)r_m[i], a0_m, er_model, KatzPoint_r_min_m, max_electron_range_m, alpha, Katz_LinearLawER_plateau_Gy, Katz_point_coeff_Gy);
      }

      // end if ER_ButtsKatz
    } else {
      for( i = 0 ; i < n ; i++){
        D_RDD_Gy[i] = 0.0f;
      }
    }

  } else { // RDD other that RDD_KatzPoint

    for( i = 0 ; i < n ; i++){
      D_RDD_Gy[i] = AT_RDD_ExtendedTarget_integrate_Gy(r_m[i], a0_m, r_min_m, r_max_m, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, er_parameter);
    }
  }

}


double AT_RDD_Katz_ext_integrand_Gy(  double t_m,
    void * params){

  AT_RDD_ExtendedTarget_parameters* RDD_parameters = (AT_RDD_ExtendedTarget_parameters*)params;

  const long   n              =  1;
  const float  r_m            =  RDD_parameters->r_m;
  const float  a0_m           =  RDD_parameters->a0_m;
  const float  E_MeV_u        =  RDD_parameters->E_MeV_u;
  const long   particle_no    =  RDD_parameters->particle_no;
  const long   material_no    =  RDD_parameters->material_no;
  const long   rdd_model      =  RDD_parameters->rdd_model;
  const float* rdd_parameter  =  RDD_parameters->rdd_parameter;
  const long   er_model       =  RDD_parameters->er_model;
  const float* er_parameter   =  RDD_parameters->er_parameter;

  float D_Gy;
  float t_m_float = (float)t_m;
  AT_D_RDD_Gy( n, &t_m_float, E_MeV_u,  particle_no, material_no,  rdd_model,  rdd_parameter,  er_model,  er_parameter,  &D_Gy);

  return 2.0 * t_m * (double)D_Gy * geometryFunctionPhi(r_m, a0_m, t_m);
}

double AT_RDD_ExtendedTarget_integrate_Gy(  const double r_m,
    const double a0_m,
    const double r_min_m,
    const double r_max_m,
    const double E_MeV_u,
    const long   particle_no,
    /* detector parameters */
    const long   material_no,
    /* radial dose distribution model */
    const long   rdd_model,
    const float  rdd_parameter[],
    /* electron range model */
    const long   er_model,
    const float  er_parameter[]){

  double low_lim_m = 0.0;
  if( r_m > a0_m ){
    low_lim_m = GSL_MAX(r_m - a0_m,r_min_m);
  }
  gsl_set_error_handler_off();

  double ext_integral_Gy;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &AT_RDD_Katz_ext_integrand_Gy;

  AT_RDD_ExtendedTarget_parameters RDD_parameters;

  RDD_parameters.r_m           =  r_m;
  RDD_parameters.a0_m          =  a0_m;
  RDD_parameters.E_MeV_u       =  (float)E_MeV_u;
  RDD_parameters.particle_no   =  particle_no;
  RDD_parameters.material_no   =  material_no;
  RDD_parameters.rdd_model     =  rdd_model;
  RDD_parameters.rdd_parameter =  (float*)rdd_parameter;
  RDD_parameters.er_model      =  er_model;
  RDD_parameters.er_parameter  =  (float*)er_parameter;

  F.params = (void*)(&RDD_parameters);
  int status = gsl_integration_qags (&F, low_lim_m, r_m+a0_m, 0, 1e-5, 1000, w1, &ext_integral_Gy, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
    printf("Error in AT_RDD_ExtendedTarget_integrate_Gy\n");
    ext_integral_Gy = -1.0;
  }
  gsl_integration_workspace_free (w1);

  ext_integral_Gy *= M_1_PI / gsl_pow_2(a0_m);

  return ext_integral_Gy;
}
