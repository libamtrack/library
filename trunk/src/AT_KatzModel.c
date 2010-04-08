/**
 * @file
 * @brief ...
 */

/*
 *    AT_KatzModel.c
 *    ===========================
 *
 *    Created on: 01.03.2010
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

#include "AT_KatzModel.h"


double AT_KatzModel_KatzExtTarget_inactivation_probability(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const long    er_model,
    const double  alpha,
    const double  Katz_plateau_Gy,
    const double  Katz_point_coeff_Gy,
    const double  D0_characteristic_dose_Gy,
    const double  c_hittedness,
    const double  m_number_of_targets){

  if( (r_m >= 0.0) && (r_m <= a0_m + max_electron_range_m ) ){
    double        D_Gy = AT_RDD_ExtendedTarget_KatzPoint_Gy( r_m, a0_m, er_model, KatzPoint_r_min_m, max_electron_range_m, alpha, Katz_plateau_Gy, Katz_point_coeff_Gy);

    const long    gamma_model                =  GR_GeneralTarget;
    const double  gamma_parameters[5]        =  { 1.0, D0_characteristic_dose_Gy, c_hittedness, m_number_of_targets, 0.0};
    double        inactivation_probability;
    const long    n_tmp                      =  1;
    AT_gamma_response(  n_tmp,
        &D_Gy,
        gamma_model,
        gamma_parameters,
        &inactivation_probability);

    return inactivation_probability;
  } else {
    return 0.0;
  }
}


int AT_KatzModel_inactivation_probability(
    const long    n,
    const double  r_m[],
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameters[],
    const long    er_model,
    const double  gamma_parameters[5],
    double        inactivation_probability[]){

  const double max_electron_range_m  =  AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);

  const double a0_m                  =  AT_RDD_a0_m( max_electron_range_m, rdd_model, rdd_parameters );

  const double KatzPoint_r_min_m     =  AT_RDD_r_min_m( max_electron_range_m, rdd_model, rdd_parameters );

  const double Katz_point_coeff_Gy   =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);

  const double r_max_m               =  GSL_MIN(a0_m, max_electron_range_m);

  double Katz_plateau_Gy  =  0.0;
  double alpha                       =  0.0;
  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
    alpha            =  AT_ER_PowerLaw_alpha(E_MeV_u);
    Katz_plateau_Gy  =  AT_RDD_Katz_PowerLawER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, alpha, Katz_point_coeff_Gy );
  } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
    Katz_plateau_Gy  =  AT_RDD_Katz_LinearER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, Katz_point_coeff_Gy );
  }

  const double D0_characteristic_dose_Gy  =  gamma_parameters[1];
  const double c_hittedness               =  gamma_parameters[2];
  const double m_number_of_targets        =  gamma_parameters[3];

  if( rdd_model == RDD_KatzExtTarget ){
    long i;
    for( i = 0 ; i < n ; i++){
      inactivation_probability[i] = AT_KatzModel_KatzExtTarget_inactivation_probability(
          r_m[i],
          a0_m,
          KatzPoint_r_min_m,
          max_electron_range_m,
          er_model,
          alpha,
          Katz_plateau_Gy,
          Katz_point_coeff_Gy,
          D0_characteristic_dose_Gy,
          c_hittedness,
          m_number_of_targets);
    }
  }

  return 0;
}


double AT_KatzModel_KatzExtTarget_inactivation_cross_section_integrand_m(
    double t_m,
    void* params){

  AT_KatzModel_KatzExtTarget_inactivation_probability_parameters * inact_prob_parameters = (AT_KatzModel_KatzExtTarget_inactivation_probability_parameters*)params;

  const double  inactivation_probability   =  AT_KatzModel_KatzExtTarget_inactivation_probability( t_m,
      inact_prob_parameters->a0_m,
      inact_prob_parameters->KatzPoint_r_min_m,
      inact_prob_parameters->max_electron_range_m,
      inact_prob_parameters->er_model,
      inact_prob_parameters->alpha,
      inact_prob_parameters->Katz_plateau_Gy,
      inact_prob_parameters->Katz_point_coeff_Gy,
      inact_prob_parameters->D0_characteristic_dose_Gy,
      inact_prob_parameters->c_hittedness,
      inact_prob_parameters->m_number_of_targets);

  return inactivation_probability * t_m;
}


double AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2(
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const long    er_model,
    const double  alpha,
    const double  Katz_plateau_Gy,
    const double  Katz_point_coeff_Gy,
    const double  D0_characteristic_dose_Gy,
    const double  c_hittedness,
    const double  m_number_of_targets){

  double low_lim_m = 0.0;
  gsl_set_error_handler_off();

  double integral_m2;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &AT_KatzModel_KatzExtTarget_inactivation_cross_section_integrand_m;

  AT_KatzModel_KatzExtTarget_inactivation_probability_parameters inact_prob_parameters;

  inact_prob_parameters.a0_m                       =  a0_m;
  inact_prob_parameters.KatzPoint_r_min_m          =  KatzPoint_r_min_m;
  inact_prob_parameters.max_electron_range_m       =  max_electron_range_m;
  inact_prob_parameters.er_model                   =  er_model;
  inact_prob_parameters.alpha                      =  alpha;
  inact_prob_parameters.Katz_plateau_Gy            =  Katz_plateau_Gy;
  inact_prob_parameters.Katz_point_coeff_Gy        =  Katz_point_coeff_Gy;
  inact_prob_parameters.D0_characteristic_dose_Gy  =  D0_characteristic_dose_Gy;
  inact_prob_parameters.c_hittedness               =  c_hittedness;
  inact_prob_parameters.m_number_of_targets        =  m_number_of_targets;

  F.params = (void*)(&inact_prob_parameters);
  int status = gsl_integration_qags (&F, low_lim_m, max_electron_range_m + a0_m, 0, 1e-5, 1000, w1, &integral_m2, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
    printf("Error in AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2\n");
    integral_m2 = -1.0;
  }
  gsl_integration_workspace_free (w1);

  return 2.0 * M_PI * integral_m2;
}


double AT_KatzModel_inactivation_cross_section_m2(
    const double E_MeV_u,
    const long   particle_no,
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double gamma_parameters[5]){
  return 0.0;
}


double AT_D_RDD_Gy_int( double  r_m,
    void*   params){
  double D_Gy;
  long  n_tmp              = 1;
  AT_P_RDD_parameters* par = (AT_P_RDD_parameters*)params;
  double fr_m               = r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  n_tmp,
      &fr_m,
      *(par->E_MeV_u),
      *(par->particle_no),
      *(par->material_no),
      *(par->rdd_model),
      par->rdd_parameters,
      *(par->er_model),
      par->er_parameters,
      &D_Gy);

  return (2.0 * M_PI * r_m * D_Gy);
}


double AT_sI_int( double  r_m,
    void*   params){
  double  P = AT_P_RDD(r_m, params);
  return (r_m * P);
}


double AT_P_RDD( double  r_m,
    void*   params)
{
  double  D_Gy;
  long   n_tmp              = 1;
  AT_P_RDD_parameters* par  = (AT_P_RDD_parameters*)params;
  double  fr_m               = r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  n_tmp,
      &fr_m,
      *(par->E_MeV_u),
      *(par->particle_no),
      *(par->material_no),
      *(par->rdd_model),
      par->rdd_parameters,
      *(par->er_model),
      par->er_parameters,
      &D_Gy);

  long gamma_model = GR_GeneralTarget;
  double P;
  AT_gamma_response(  n_tmp,
      &D_Gy,
      gamma_model,
      par->gamma_parameters,
      // return
      &P);
  return P;
}
