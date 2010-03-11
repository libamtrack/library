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
void AT_RDD_Site_Gy( const long*  n,
    const float*  r_m,
    const float*  a0_m,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const long*   particle_no,
    /* detector parameters */
    const long*   material_no,
    /* radial dose distribution model */
    const long*   rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const long*   er_model,
    const float*  er_parameter,
    float*        D_RDD_Gy){

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

void AT_RDD_ExtendedTarget_Gy( const long  n,
    const float* r_m,
    const float  a0_m,
    /* radiation field parameters */
    const float  E_MeV_u,
    const long   particle_no,
    /* detector parameters */
    const long   material_no,
    /* radial dose distribution model */
    const long   rdd_model,
    const float* rdd_parameter,
    /* electron range model */
    const long   er_model,
    const float* er_parameter,
    float*       D_RDD_Gy){

  // Get material data
  double   density_g_cm3, density_kg_m3, electron_density_m3;
  AT_get_material_data(
              material_no,
              &density_g_cm3,
              &electron_density_m3,
              NULL, NULL, NULL, NULL, NULL, NULL);
  density_kg_m3      =  density_g_cm3 * 1000.0f;

  // Get beta, Z and Zeff
  float  beta   =  0.0;
  long  Z       =  0;
  float  Z_eff  =  0.0;
  long  n_tmp    = 1;
  AT_beta_from_E(  &n_tmp,
                &E_MeV_u,
                &beta);

  AT_Z_from_particle_no(  &n_tmp,
                &particle_no,
                &Z);
  AT_effective_charge_from_beta(  &n_tmp,
                  &beta,
                  &Z,
                  &Z_eff);

  // Get the maximum electron range
  float  max_electron_range_m      =  (double)AT_max_electron_range_m((double)E_MeV_u, (int)material_no, (int)er_model);

  const double C_J_m               =  AT_RDD_Katz_C_J_m((double)electron_density_m3);

  // Get LET
  float  LET_MeV_cm2_g  =  0.0f;
  AT_LET_MeV_cm2_g(  &n_tmp,
            &E_MeV_u,
            &particle_no,
            &material_no,
            &LET_MeV_cm2_g);

  float  LET_J_m        =  LET_MeV_cm2_g * density_g_cm3; // [MeV / cm]
  LET_J_m              *=  100.0f;       // [MeV / m]
  LET_J_m              *=  MeV_to_J;     // [J/m]

  double Katz_point_coeff_Gy =  AT_RDD_Katz_coeff_Gy(C_J_m,(double)Z_eff,(double)beta,(double)density_kg_m3,(double)max_electron_range_m);
  float r_max_m = GSL_MIN((double)a0_m, max_electron_range_m);

  float alpha           =  0.0f;
  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
    alpha               =  1.667f;
    float wmax_MeV      =  0.0f;
    AT_max_E_transfer_MeV(&n_tmp,&E_MeV_u,&wmax_MeV);
    if(wmax_MeV <= 1e-3){
      alpha             =  1.079f;
    }
  }

  /********************************************************
   ****************** LOOP OVER DISTANCES *****************
   *******************************************************/
  long     i;

  float r_min_m = rdd_parameter[0];

  for( i = 0 ; i < n ; i++){

    D_RDD_Gy[i] = 0.0f;

    if( rdd_model == RDD_KatzPoint){ // TODO write comments

      if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD

        if( (r_m[i] <=  0.01 * a0_m) && (r_m[i] >=  r_min_m)){
          D_RDD_Gy[i] = (float)AT_RDD_Katz_PowerLawER_Daverage_Gy( (double)r_min_m, (double)r_max_m, (double)max_electron_range_m, (double)alpha, Katz_point_coeff_Gy );
          if( max_electron_range_m < a0_m ){
            D_RDD_Gy[i] *= (gsl_pow_2(r_max_m) - gsl_pow_2(r_min_m));
            D_RDD_Gy[i] /= (gsl_pow_2(a0_m) - gsl_pow_2(r_min_m));
          }
        }

        if( (r_m[i] >=  100.0 * a0_m) && (r_m[i] <= max_electron_range_m) ){
          D_RDD_Gy[i] = (float)AT_RDD_Katz_PowerLawER_Dpoint_Gy( (double)r_m[i], (double)alpha, (double)max_electron_range_m, Katz_point_coeff_Gy );
        }

        if( (r_m[i] <  100.0 * a0_m) && (r_m[i] >  0.01 * a0_m) ){
          D_RDD_Gy[i] = AT_RDD_ExtendedTarget_integrate_Gy(r_m[i],a0_m,r_min_m,r_max_m,E_MeV_u,particle_no,material_no,rdd_model,rdd_parameter,er_model,er_parameter);
        }

        // end if ER_Waligorski, ER_Edmund
      } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD

        if( (r_m[i] <=  0.01 * a0_m) && (r_m[i] >=  r_min_m)){
          D_RDD_Gy[i]  =  (float)AT_RDD_Katz_LinearER_Daverage_Gy( (double)r_min_m, (double)r_max_m, (double)max_electron_range_m, Katz_point_coeff_Gy );
          if( max_electron_range_m < a0_m ){
            D_RDD_Gy[i] *= (gsl_pow_2(r_max_m) - gsl_pow_2(r_min_m));
            D_RDD_Gy[i] /= (gsl_pow_2(a0_m) - gsl_pow_2(r_min_m));
          }
        }

        if( (r_m[i] >=  100.0 * a0_m) && (r_m[i] <= max_electron_range_m) ){
          D_RDD_Gy[i] = (float)AT_RDD_Katz_LinearER_Dpoint_Gy( (double)r_m[i], (double)max_electron_range_m, Katz_point_coeff_Gy );
        }

        if( (r_m[i] <  100.0 * a0_m) && (r_m[i] >  0.01 * a0_m) ){
          D_RDD_Gy[i] = AT_RDD_ExtendedTarget_integrate_Gy(r_m[i],a0_m,r_min_m,r_max_m,E_MeV_u,particle_no,material_no,rdd_model,rdd_parameter,er_model,er_parameter);
        }

        // end if ER_ButtsKatz
      } else {
        D_RDD_Gy[i] = 0.0f;
      }

    } else { // RDD other that RDD_KatzPoint
      D_RDD_Gy[i] = AT_RDD_ExtendedTarget_integrate_Gy(r_m[i],a0_m,r_min_m,r_max_m,E_MeV_u,particle_no,material_no,rdd_model,rdd_parameter,er_model,er_parameter);
    }

  }// end loop over r_m

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
  AT_D_RDD_Gy( &n, &t_m_float, &E_MeV_u,  &particle_no, &material_no,  &rdd_model,  rdd_parameter,  &er_model,  er_parameter,  &D_Gy);

  return 2.0 * t_m * (double)D_Gy * geometryFunctionPhi(r_m,a0_m,t_m);
}

double AT_RDD_ExtendedTarget_integrate_Gy(  const double r_m,
    const double a0_m,
    const double r_min_m,
    const double r_max_m,
    const float  E_MeV_u,
    const long   particle_no,
    /* detector parameters */
    const long   material_no,
    /* radial dose distribution model */
    const long   rdd_model,
    const float* rdd_parameter,
    /* electron range model */
    const long   er_model,
    const float* er_parameter){

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
  RDD_parameters.E_MeV_u       =  E_MeV_u;
  RDD_parameters.particle_no   =  particle_no;
  RDD_parameters.material_no   =  material_no;
  RDD_parameters.rdd_model     =  rdd_model;
  RDD_parameters.rdd_parameter =  rdd_parameter;
  RDD_parameters.er_model      =  er_model;
  RDD_parameters.er_parameter  =  er_parameter;

  F.params = (void*)(&RDD_parameters);
  int status = gsl_integration_qags (&F, low_lim_m, r_m+a0_m, 0, 1e-5, 1000, w1, &ext_integral_Gy, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
    printf("Error\n");
    ext_integral_Gy = -1.0;
  }
  gsl_integration_workspace_free (w1);

  ext_integral_Gy *= M_1_PI / gsl_pow_2(a0_m);

  return ext_integral_Gy;
}
