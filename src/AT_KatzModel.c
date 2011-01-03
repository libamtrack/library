/**
 * @brief TODO
 */

/*
 *    AT_KatzModel.c
 *    ===========================
 *
 *    Created on: 01.03.2010
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
        false,
        &inactivation_probability);

    return inactivation_probability;
  } else {
    return 0.0;
  }
}


double AT_KatzModel_CucinottaExtTarget_inactivation_probability(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy,
    const double  KatzPoint_point_coeff_Gy,
    const double  D0_characteristic_dose_Gy,
    const double  c_hittedness,
    const double  m_number_of_targets){

  if( (r_m >= 0.0) && (r_m <= a0_m + max_electron_range_m ) ){
    double        D_Gy = AT_RDD_ExtendedTarget_CucinottaPoint_Gy( r_m, a0_m, KatzPoint_r_min_m, max_electron_range_m, beta, KatzPoint_point_coeff_Gy, C_norm, Cucinotta_plateau_Gy);

    const long    gamma_model                =  GR_GeneralTarget;
    const double  gamma_parameters[5]        =  { 1.0, D0_characteristic_dose_Gy, c_hittedness, m_number_of_targets, 0.0};
    double        inactivation_probability;
    const long    n_tmp                      =  1;
    AT_gamma_response(  n_tmp,
        &D_Gy,
        gamma_model,
        gamma_parameters,
        false,
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
    const double  gamma_parameters[],
    double        inactivation_probability[]){

  const double max_electron_range_m  =  AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);

  const double a0_m                  =  AT_RDD_a0_m( max_electron_range_m, rdd_model, rdd_parameters );

  const double KatzPoint_r_min_m     =  AT_RDD_r_min_m( max_electron_range_m, rdd_model, rdd_parameters );

  const double Katz_point_coeff_Gy   =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);

  const double r_max_m               =  GSL_MIN(a0_m, max_electron_range_m);

  const double D0_characteristic_dose_Gy  =  gamma_parameters[1];
  const double c_hittedness               =  gamma_parameters[2];
  const double m_number_of_targets        =  gamma_parameters[3];

  if( rdd_model == RDD_KatzExtTarget ){
    double Katz_plateau_Gy  =  0.0;
    double alpha                       =  0.0;
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
      alpha            =  AT_ER_PowerLaw_alpha(E_MeV_u);
      Katz_plateau_Gy  =  AT_RDD_Katz_PowerLawER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, alpha, Katz_point_coeff_Gy );
    } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
      Katz_plateau_Gy  =  AT_RDD_Katz_LinearER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, Katz_point_coeff_Gy );
    }

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

  if( rdd_model == RDD_CucinottaExtTarget ){
    const double  density_g_cm3        =  AT_density_g_cm3_from_material_no( material_no );
    const double  density_kg_m3        =  density_g_cm3 * 1000.0;
    const double  LET_MeV_cm2_g        =  AT_LET_MeV_cm2_g_single(E_MeV_u, particle_no, material_no);
    const double  LET_J_m              =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
    const double  beta                 =  AT_beta_from_E_single( E_MeV_u );
    const double  C_norm               =  AT_RDD_Cucinotta_Cnorm(KatzPoint_r_min_m, max_electron_range_m, beta, density_kg_m3, LET_J_m, Katz_point_coeff_Gy);
    double Cucinotta_plateau_Gy        =  AT_RDD_Cucinotta_Ddelta_average_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
    Cucinotta_plateau_Gy              +=  C_norm * AT_RDD_Cucinotta_Dexc_average_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);

    long i;
    for( i = 0 ; i < n ; i++){
      inactivation_probability[i] = AT_KatzModel_CucinottaExtTarget_inactivation_probability(
          r_m[i],
          a0_m,
          KatzPoint_r_min_m,
          max_electron_range_m,
          beta,
          C_norm,
          Cucinotta_plateau_Gy,
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

  assert( params != NULL );
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

  double low_lim_m = 0.01*a0_m;
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
  int status = gsl_integration_qag (&F, low_lim_m, max_electron_range_m + a0_m, 0, 1e-4, 1000, GSL_INTEG_GAUSS21, w1, &integral_m2, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
    printf("Error in AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2: er_model = %ld, integration from %g to %g [m] + %g [m]\n", er_model, low_lim_m, max_electron_range_m, a0_m);
    integral_m2 = 0.0;
  }
  gsl_integration_workspace_free (w1);

  const long    gamma_model                =  GR_GeneralTarget;
  const double  gamma_parameters[5]        =  { 1.0, D0_characteristic_dose_Gy, c_hittedness, m_number_of_targets, 0.0};
  double        inactivation_probability_plateau;
  const long    n_tmp                      =  1;
  AT_gamma_response(  n_tmp,
      &Katz_plateau_Gy,
      gamma_model,
      gamma_parameters,
      false,
      &inactivation_probability_plateau);

  return 2.0 * M_PI * integral_m2 + M_PI * gsl_pow_2(0.01*a0_m) * inactivation_probability_plateau;
}


double AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_integrand_m(
    double t_m,
    void* params){

  assert( params != NULL );
  AT_KatzModel_CucinottaExtTarget_inactivation_probability_parameters * inact_prob_parameters = (AT_KatzModel_CucinottaExtTarget_inactivation_probability_parameters*)params;

  const double  inactivation_probability   =  AT_KatzModel_CucinottaExtTarget_inactivation_probability( t_m,
      inact_prob_parameters->a0_m,
      inact_prob_parameters->KatzPoint_r_min_m,
      inact_prob_parameters->max_electron_range_m,
      inact_prob_parameters->beta,
      inact_prob_parameters->C_norm,
      inact_prob_parameters->Cucinotta_plateau_Gy,
      inact_prob_parameters->KatzPoint_coeff_Gy,
      inact_prob_parameters->D0_characteristic_dose_Gy,
      inact_prob_parameters->c_hittedness,
      inact_prob_parameters->m_number_of_targets);

  return inactivation_probability * t_m;
}


double AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2(
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy,
    const double  KatzPoint_point_coeff_Gy,
    const double  D0_characteristic_dose_Gy,
    const double  c_hittedness,
    const double  m_number_of_targets){

  double low_lim_m = 0.01*a0_m;
  gsl_set_error_handler_off();

  double integral_m2;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_integrand_m;

  AT_KatzModel_CucinottaExtTarget_inactivation_probability_parameters inact_prob_parameters;

  inact_prob_parameters.a0_m                       =  a0_m;
  inact_prob_parameters.KatzPoint_r_min_m          =  KatzPoint_r_min_m;
  inact_prob_parameters.max_electron_range_m       =  max_electron_range_m;
  inact_prob_parameters.beta                       =  beta;
  inact_prob_parameters.C_norm                     =  C_norm;
  inact_prob_parameters.Cucinotta_plateau_Gy       =  Cucinotta_plateau_Gy;
  inact_prob_parameters.KatzPoint_coeff_Gy         =  KatzPoint_point_coeff_Gy;
  inact_prob_parameters.D0_characteristic_dose_Gy  =  D0_characteristic_dose_Gy;
  inact_prob_parameters.c_hittedness               =  c_hittedness;
  inact_prob_parameters.m_number_of_targets        =  m_number_of_targets;

  F.params = (void*)(&inact_prob_parameters);
  int status = gsl_integration_qag (&F, low_lim_m, max_electron_range_m + a0_m, 0, 1e-4, 1000, GSL_INTEG_GAUSS21, w1, &integral_m2, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
    printf("Error in AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2: integration from %g to %g [m] + %g [m]\n", low_lim_m, max_electron_range_m, a0_m);
    integral_m2 = 0.0;
  }
  gsl_integration_workspace_free (w1);

  const long    gamma_model                =  GR_GeneralTarget;
  const double  gamma_parameters[5]        =  { 1.0, D0_characteristic_dose_Gy, c_hittedness, m_number_of_targets, 0.0};
  double        inactivation_probability_plateau;
  const long    n_tmp                      =  1;
  AT_gamma_response(  n_tmp,
      &Cucinotta_plateau_Gy,
      gamma_model,
      gamma_parameters,
      false,
      &inactivation_probability_plateau);

 return 2.0 * M_PI * integral_m2 + M_PI * gsl_pow_2(0.01*a0_m) * inactivation_probability_plateau;

}


int AT_KatzModel_inactivation_cross_section_m2(
    const long   n,
    const double E_MeV_u[],
    const long   particle_no,
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double gamma_parameters[],
    double inactivation_cross_section_m2[]){

  const double D0_characteristic_dose_Gy  =  gamma_parameters[1];
  const double c_hittedness               =  gamma_parameters[2];
  const double m_number_of_targets        =  gamma_parameters[3];

  if( rdd_model == RDD_KatzExtTarget ){
    long i;
    for( i = 0 ; i < n ; i++){

      const double max_electron_range_m  =  AT_max_electron_range_m( E_MeV_u[i], (int)material_no, (int)er_model);
      const double a0_m                  =  rdd_parameters[1];
      const double KatzPoint_r_min_m     =  AT_RDD_r_min_m( max_electron_range_m, rdd_model, rdd_parameters );
      const double Katz_point_coeff_Gy   =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u[i], particle_no, material_no, er_model);
      const double r_max_m               =  GSL_MIN(a0_m, max_electron_range_m);

      double Katz_plateau_Gy  =  0.0;
      double alpha                       =  0.0;
      if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
        alpha            =  AT_ER_PowerLaw_alpha(E_MeV_u[i]);
        Katz_plateau_Gy  =  AT_RDD_Katz_PowerLawER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, alpha, Katz_point_coeff_Gy );
      } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
        Katz_plateau_Gy  =  AT_RDD_Katz_LinearER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, Katz_point_coeff_Gy );
      }

      inactivation_cross_section_m2[i] = AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2(
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

  if( rdd_model == RDD_CucinottaExtTarget ){
    long i;
    const double  density_g_cm3        =  AT_density_g_cm3_from_material_no( material_no );
    const double  density_kg_m3        =  density_g_cm3 * 1000.0;

    for( i = 0 ; i < n ; i++){

      const double max_electron_range_m  =  AT_max_electron_range_m( E_MeV_u[i], (int)material_no, (int)er_model);
      const double a0_m                  =  rdd_parameters[1]; // AT_RDD_a0_m( max_electron_range_m, rdd_model, rdd_parameters );
      const double KatzPoint_r_min_m     =  AT_RDD_r_min_m( max_electron_range_m, rdd_model, rdd_parameters );
      const double Katz_point_coeff_Gy   =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u[i], particle_no, material_no, er_model);
      const double r_max_m               =  GSL_MIN(a0_m, max_electron_range_m);

      const double  LET_MeV_cm2_g        =  AT_LET_MeV_cm2_g_single(E_MeV_u[i], particle_no, material_no);
      const double  LET_J_m              =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
      const double  beta                 =  AT_beta_from_E_single( E_MeV_u[i] );
      const double  C_norm               =  AT_RDD_Cucinotta_Cnorm(KatzPoint_r_min_m, max_electron_range_m, beta, density_kg_m3, LET_J_m, Katz_point_coeff_Gy);
      double Cucinotta_plateau_Gy        =  AT_RDD_Cucinotta_Ddelta_average_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
      Cucinotta_plateau_Gy              +=  C_norm * AT_RDD_Cucinotta_Dexc_average_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);

      inactivation_cross_section_m2[i] = AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2(
          a0_m,
          KatzPoint_r_min_m,
          max_electron_range_m,
          beta,
          C_norm,
          Cucinotta_plateau_Gy,
          Katz_point_coeff_Gy,
          D0_characteristic_dose_Gy,
          c_hittedness,
          m_number_of_targets);
    }
  }


  return 0;
}


/* TODO implement old Katz with kappa and sigma instead of track-width here */
double AT_KatzModel_single_field_survival(
    const double fluence_cm2,
	const double E_MeV_u,
    const long   particle_no,
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double D0_characteristic_dose_Gy,
    const double m_number_of_targets,
    const double sigma0_m2){

	/* some useful variables */
	double dose_Gy = AT_dose_Gy_from_fluence_cm2_single( E_MeV_u, particle_no, fluence_cm2, material_no); /* fluence + LET -> dose */

	assert( sigma0_m2 > 0);

	/* single particle inactivation cross section calculation */
	double inactivation_cross_section_m2 = 0.0;
	double gamma_parameters[5] = {1.,D0_characteristic_dose_Gy,1.,m_number_of_targets,0.};
	AT_KatzModel_inactivation_cross_section_m2(
	    1,
	    &E_MeV_u,
	    particle_no,
	    material_no,
	    rdd_model,
	    rdd_parameters,
	    er_model,
	    gamma_parameters,
	    &inactivation_cross_section_m2);    /* here we use D0, m and a0 */
	printf("inactivation_cross_section = %g [m2]\n", inactivation_cross_section_m2);

	/* fraction of dose delivered in ion kill mode */
	double ion_kill_mode_fraction = inactivation_cross_section_m2 / sigma0_m2;
	if( ion_kill_mode_fraction > 1.0) ion_kill_mode_fraction = 1.0;
	printf("ion_kill_mode_fraction = %g\n", ion_kill_mode_fraction);

	double gamma_kill_dose = (1. - ion_kill_mode_fraction) * dose_Gy;


	double gamma_kill_mode_survival;
	double ion_kill_mode_survival;

	/* ion kill mode survival gives exponential SF curve */
	ion_kill_mode_survival = exp( - inactivation_cross_section_m2 * 1e4 * fluence_cm2);  /* depends on D0, m and a0; not on kappa !*/

	assert( D0_characteristic_dose_Gy > 0 );

	/* gamma kill mode survival gives shouldered SF curve */
	if( ion_kill_mode_fraction > 0.98 ){
		gamma_kill_mode_survival = 1.0;
		/* this if is quite artificial, for ion_kill_mode_fraction greater than 0.98
		 * gamma_kill_mode_survival is so close to 1 that it does not make sense to calculate it
		 */
	} else {
		gamma_kill_mode_survival = 1. - pow( 1. - exp( - gamma_kill_dose / D0_characteristic_dose_Gy), m_number_of_targets);  /* depends on D0, m and kappa; not on a0 ! */
	}

	printf("ion_kill_mode_survival = %g , gamma_kill_mode_survival = %g[m2]\n", ion_kill_mode_survival, gamma_kill_mode_survival);

	/* finally survival as a product of ion and gamma kill modes */
	return ion_kill_mode_survival * gamma_kill_mode_survival;
}


double AT_D_RDD_Gy_int( double  r_m,
    void*   params){
  double D_Gy;
  long  n_tmp              = 1;

  assert( params != NULL );

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
      &D_Gy);

  return (2.0 * M_PI * r_m * D_Gy);
}


double AT_sI_int( double  r_m,
    void*   params){
  assert( params != NULL );
  double  P = AT_P_RDD(r_m, params);
  return (r_m * P);
}


double AT_P_RDD( double  r_m,
    void*   params)
{
  double  D_Gy;
  long   n_tmp              = 1;
  assert( params != NULL );
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
      &D_Gy);

  long gamma_model = GR_GeneralTarget;
  double P;
  AT_gamma_response(  n_tmp,
      &D_Gy,
      gamma_model,
      par->gamma_parameters,
      false,
      // return
      &P);
  return P;
}
