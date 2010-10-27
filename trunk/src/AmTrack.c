/**
 * @brief main file holding the amorphous track methods for RE/RBE calculation
 */

/*
 *    AmTrack.c
 *    =========
 *
 *    Created on: 28.07.2009
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

#include "AmTrack.h"
#include <math.h>

void AT_run_CPPSC_method(  const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameters[],
    const long    er_model,
    const long    gamma_model,
    const double  gamma_parameters[],
    long          N2, // TODO investigate if this can be changed inside
    const double  fluence_factor,
    const bool    write_output,
    const bool    shrink_tails,
    const double  shrink_tails_under,
    const bool    adjust_N2,
    const bool    lethal_events_mode,
    double*       relative_efficiency,
    double*       d_check,
    double*       S_HCP,
    double*       S_gamma,
    double*       mean_number_of_tracks_contrib,
    double*       start_number_of_tracks_contrib,
    long*         n_convolutions,
    double*		  lower_Jensen_bound,
    double*       upper_Jensen_bound)
{
  long i;


  long     n_bins_f1 = AT_n_bins_for_singe_impact_local_dose_distrib(  number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      N2);

  double*  f1_parameters      =  (double*)calloc(AT_SC_F1_PARAMETERS_SINGLE_LENGTH * number_of_field_components, sizeof(double));

  AT_RDD_f1_parameters_mixed_field( number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      f1_parameters
  );

  double*  f1_d_Gy       =  (double*)calloc(n_bins_f1, sizeof(double));
  double*  f1_dd_Gy      =  (double*)calloc(n_bins_f1, sizeof(double));
  double*  f1            =  (double*)calloc(n_bins_f1, sizeof(double));

  AT_single_impact_local_dose_distrib(  number_of_field_components,
      E_MeV_u,
      particle_no,
      fluence_cm2_or_dose_Gy,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      N2,
      n_bins_f1,
      f1_parameters,
      f1_d_Gy,
      f1_dd_Gy,
      f1);

  double*  fluence_cm2    =  (double*)calloc(number_of_field_components, sizeof(double));

  if(fluence_cm2_or_dose_Gy[0] < 0){
    double*  dose_Gy      =  (double*)calloc(number_of_field_components, sizeof(double));
    for (i = 0; i < number_of_field_components; i++){
      dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    // convert dose to fluence
    AT_fluence_cm2(  number_of_field_components,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        fluence_cm2);
    free( dose_Gy );
  }else{
    for (i = 0; i < number_of_field_components; i++){
      fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
    }
  }

  *mean_number_of_tracks_contrib  =       AT_total_u(     number_of_field_components,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      er_model);

  free( fluence_cm2 );

  long      n_bins_f;

  AT_n_bins_for_low_fluence_local_dose_distribution(  *mean_number_of_tracks_contrib,
      fluence_factor,
      N2,
      n_bins_f1,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      // from here: return values
      &n_bins_f,
      start_number_of_tracks_contrib,
      n_convolutions);

  double*  f_d_Gy       =  (double*)calloc(n_bins_f, sizeof(double));
  double*  f_dd_Gy      =  (double*)calloc(n_bins_f, sizeof(double));
  double*  f            =  (double*)calloc(n_bins_f, sizeof(double));
  double*  fdd          =  (double*)calloc(n_bins_f, sizeof(double));
  double*  dfdd         =  (double*)calloc(n_bins_f, sizeof(double));
  double   f0           =  0.0;

  AT_low_fluence_local_dose_distribution(  n_bins_f1,
      N2,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      n_bins_f,
      f_d_Gy,
      f_dd_Gy,
      f);

  AT_SuccessiveConvolutions(  *mean_number_of_tracks_contrib,
      n_bins_f,
      &N2,
      // input + return values
      &n_bins_f1,
      f_d_Gy,
      f_dd_Gy,
      f,
      // return values
      &f0,
      fdd,
      dfdd,
      d_check,
      write_output,
      shrink_tails,
      shrink_tails_under,
      adjust_N2);

  long     n_bins_f_used  	= n_bins_f1;
  double*  S            	= (double*)calloc(n_bins_f_used, sizeof(double));

  AT_get_response_distribution_from_dose_distribution(  n_bins_f_used,
      f_d_Gy,
      f,
      gamma_model,
      gamma_parameters,
      lethal_events_mode,
      S);

  *S_HCP = AT_get_ion_response_from_response_distribution( n_bins_f_used,
		  f_dd_Gy,
		  f,
		  S);

  *S_gamma = AT_get_gamma_response_for_average_dose( n_bins_f_used,
		  f_d_Gy,
		  f_dd_Gy,
		  f,
		  gamma_model,
		  gamma_parameters,
		  lethal_events_mode);

  *relative_efficiency = *S_HCP / *S_gamma;


  /* Get zero-dose reponse and its ln */
  double		s0 = 0.0, log_s0 = 0.0;
  const long 	number_of_bins = 1;
  const double	d0 = 0.0;

  AT_gamma_response(  number_of_bins,
      &d0,
      gamma_model,
      gamma_parameters,
      false,
      // return
      &s0);

  if(s0 > 0){
	  log_s0 = log(s0);
  }

  /* Compute lower bound */
  *lower_Jensen_bound = f0*log_s0;
  for (i = 0; i < n_bins_f_used; i++){
	  double	log_s	= 0.0;
	  if (S[i] > 0){
		  log_s 	= log(S[i]);
	  }
	  *lower_Jensen_bound 	+= f[i] * log_s * f_dd_Gy[i];
  }
  *lower_Jensen_bound	=	exp(*lower_Jensen_bound);

  /* Compute upper bound */
  *upper_Jensen_bound = f0*s0;
  for (i = 0; i < n_bins_f_used; i++){
	  *upper_Jensen_bound 	+= f[i] * S[i] * f_dd_Gy[i];
  }

  free(f1_parameters);
  free(f1_d_Gy);
  free(f1_dd_Gy);
  free(f1);
  free(f_d_Gy);
  free(f_dd_Gy);
  free(f);
  free(fdd);
  free(dfdd);
  free(S);
}


void AT_run_IGK_method(  const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameters[],
    const long    er_model,
    const long    gamma_model,
    const double  gamma_parameters[],
    const double  saturation_cross_section_factor,
    const bool    write_output,
    double*       relative_efficiency,
    double*       S_HCP,
    double*       S_gamma,
    double*       sI_cm2,
    double*       gamma_dose_Gy,
    double*       P_I,
    double*       P_g)
{
  long   i, j;
  long   n_tmp      = 1;

  /* Clear return values */
  *relative_efficiency    = 0.0;
  *S_HCP                  = 0.0;
  *S_gamma                = 0.0;
  *sI_cm2                 = 0.0;
  *gamma_dose_Gy          = 0.0;
  *P_I                    = 0.0;
  *P_g                    = 0.0;

  /* convert dose to fluence or vice versa
   * depending on given dose/fluence value
   */
  double*  fluence_cm2    =  (double*)calloc(number_of_field_components, sizeof(double));
  double*  dose_Gy        =  (double*)calloc(number_of_field_components, sizeof(double));

  if(fluence_cm2_or_dose_Gy[0] < 0){
    for (i = 0; i < number_of_field_components; i++){
      dose_Gy[i]               = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    AT_fluence_cm2(  number_of_field_components,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        fluence_cm2);
  }else{
    for (i = 0; i < number_of_field_components; i++){
      fluence_cm2[i]          = fluence_cm2_or_dose_Gy[i];
    }
    AT_dose_Gy_from_fluence_cm2(  number_of_field_components,
        E_MeV_u,
        particle_no,
        fluence_cm2,
        material_no,
        dose_Gy);
  }

  /* Get total dose, total fluence */
  double total_dose_Gy        = 0.0;
  double total_fluence_cm2    = 0.0;

  for (i = 0; i < number_of_field_components; i++){
	  total_dose_Gy                +=  dose_Gy[i];
	  total_fluence_cm2            +=  fluence_cm2[i];
  }
  free( dose_Gy );


  /* Get normalized fluences and dose contributions
   * of each mixed field component
   */
  double u_single;
  double*  norm_fluence          =  (double*)calloc(number_of_field_components, sizeof(double));
  double*  dose_contribution_Gy  =  (double*)calloc(number_of_field_components, sizeof(double));

  /* TODO: Replace by explicit functions */
  for (i = 0; i < number_of_field_components; i++){
	  double LET_MeV_cm2_g              = AT_LET_MeV_cm2_g_single(E_MeV_u[i], particle_no[i], material_no);
	  double single_impact_fluence_cm2  = AT_single_impact_fluence_cm2_single(E_MeV_u[i], material_no, er_model);
	  norm_fluence[i]                   =  fluence_cm2[i] / total_fluence_cm2;
	  u_single                          =  fluence_cm2[i] / single_impact_fluence_cm2;
	  double single_impact_dose_Gy      = AT_single_impact_dose_Gy_single(LET_MeV_cm2_g, single_impact_fluence_cm2);
	  dose_contribution_Gy[i]           =  u_single * single_impact_dose_Gy;
  }
  free( fluence_cm2 );


  /* Get accumulated normalized fluence
   * this is used as pdf for later
   * sampling of particle type from
   * mixed field
   */
  // TODO do we really need accu_fluence ? it is not needed anywhere else
  double*  accu_fluence          =  (double*)calloc(number_of_field_components, sizeof(double));
  accu_fluence[0]                =  norm_fluence[0];
  if(number_of_field_components > 1){
    for (i = 1; i < number_of_field_components; i++){
      accu_fluence[i]                 += accu_fluence[i-1] + norm_fluence[i];
    }
  }
  free(accu_fluence);

  /* Open output file */
  //TODO rename KatseMitGlatse to something more reasonable
  FILE*    output_file = NULL;
  if( write_output ){
	  output_file          =  fopen("KatseMitGlatse.log","w");
	  if (output_file == NULL) return;                      // File error

	  fprintf(output_file, "##############################################################\n");
	  fprintf(output_file, "##############################################################\n");
	  fprintf(output_file, "This is SGP efficiency Katz, version(2009/10/08).\n");
	  fprintf(output_file, "##############################################################\n");
	  fprintf(output_file, "\n\n\n");
  }

  /* Check whether the general hit/target
   * gamma response model is used as the
   * Katz model is inheritely linked to it
   */
  if (gamma_model != GR_GeneralTarget ||
      rdd_model   == RDD_Test){
	  if( write_output ){
		  fprintf(output_file, "##############################################################\n");
		  fprintf(output_file, "Sorry, no IGK with other than the general hit-target model\n");
		  fprintf(output_file, "or with test RDD\n");
		  fprintf(output_file, "Please choose models accordingly. Exiting now...\n");
		  fprintf(output_file, "##############################################################\n");
	  }
	  return;
  }

  /* Browse given gamma parameters until terminating zero
   * and count the components (1 component = 4 parameters)
   */
  long   n_components         = 0;
  long   n_gamma_parameters   = 0;
  while  (gamma_parameters[n_gamma_parameters] != 0){
	  n_gamma_parameters           += 4;
	  n_components                 += 1;
  }

  /* Initialize variables */
  double   sI_m2               = 0.0;
  double   cross_section_ratio = 0.0;
  double   gamma_contribution  = 0.0;
  *S_HCP                       = 0.0;

  /* Create and fill the structure which
   * is needed for the RDD integration
   * by the GSL integration routine
   */
  AT_P_RDD_parameters* params;
  params                       = (AT_P_RDD_parameters*)calloc(1,sizeof(AT_P_RDD_parameters));
  params->E_MeV_u              = (double*)E_MeV_u;
  params->particle_no          = (long*)particle_no;
  params->material_no          = (long*)(&material_no);
  params->rdd_model            = (long*)(&rdd_model);
  params->rdd_parameters       = (double*)rdd_parameters;
  params->er_model             = (long*)(&er_model);
  params->gamma_parameters[0]  = 1; // No multiple components
  params->gamma_parameters[4]  = 0;

  /* Main loop though all the components
   * For each component the IGK will be individually applied
   * and later the results will be averaged fluence-weighted
   * This is due to the fact that IGK works on monoenergetic
   * situations and cannot consider true interaction
   * of particles from mixed fields
   */
  for(i = 0; i < n_components; i++){
      /**************************************************/
	  /* 1. Do integration of RDD for ion cross-section */
      /**************************************************/

	  /* Copy gamma parameters for current component
	   * into integration structure
	   */for (j = 1; j < 4; j++){
		  params->gamma_parameters[j]   = gamma_parameters[i*4 + j];
	  }

	  /* Initialize GSL integration workspace */
	  gsl_set_error_handler_off();
	  gsl_integration_workspace *w1   = gsl_integration_workspace_alloc (10000);
	  gsl_function F;
	  F.function                      = &AT_sI_int;
	  F.params                        = (void*)params;

	  /* Set integration limits */
	  double   lower_lim_m            = 0.0;
	  if(rdd_model == RDD_KatzPoint){
		  lower_lim_m                     = rdd_parameters[0];
	  }
	  // TODO energy is an array of size n, why do we calculate upper_lim_m only from one energy value ?
	  double   upper_lim_m            = AT_max_electron_range_m( *E_MeV_u, (int)material_no, (int)er_model);
	  double error;

	  /* Perform integration */
	  int status      = gsl_integration_qags (        &F,
			  lower_lim_m,
			  upper_lim_m,
			  1e-20,
			  1e-20,
			  10000,
			  w1,
			  &sI_m2,
			  &error);
	  if (status == GSL_EROUND || status == GSL_ESING){
		  printf("Error in integration (cross section calculation) - IGK\n");
	  }

	  /* Transform ion cross-section and close integration workspace */
	  sI_m2           *= 2.0 * M_PI;
	  *sI_cm2          = sI_m2 * 10000.0;
	  gsl_integration_workspace_free (w1);

	  // TODO: INTERCEPT Katz point RDD for m / c detectors here!

	  /* Get saturation cross-section from target size
	   * and sat.-cross-section factor
	   */
	  double   s0_m2   = 0.0;
	  double   a0_m    = 0.0;
	  if(rdd_model == RDD_KatzExtTarget){
		  a0_m              = rdd_parameters[1];
	  }
	  if(rdd_model == RDD_Geiss ||
			  rdd_model == RDD_KatzSite){
		  a0_m              =  rdd_parameters[0];
	  }
	  s0_m2            = saturation_cross_section_factor * M_PI * gsl_pow_2(a0_m);

	  /* Compute Ion-kill and gamma-kill probabilities
	   * The use them to compute HCP response of component
	   */
	  double   S_HCP_component, gamma_D_Gy;

	  double   fluence_cm2  = norm_fluence[0] * total_fluence_cm2;    // norm. fluence for particle i * total_fluence
	  double   D_Gy         = dose_contribution_Gy[0];                // dose by particle i

	  cross_section_ratio   = sI_m2 / s0_m2;
	  if( (cross_section_ratio < 1) & (cross_section_ratio >= 0) & (params->gamma_parameters[2] > 1)){
		  *P_I                = exp(-1.0 * (*sI_cm2) * fluence_cm2); // prob of being activated by ion kill mode
		  gamma_contribution  = 1.0 - cross_section_ratio;
		  gamma_D_Gy          = gamma_contribution * D_Gy;
		  AT_gamma_response(  n_tmp,
				  &gamma_D_Gy,
				  gamma_model,
				  params->gamma_parameters,
				  false,
				  // return
				  P_g);
		  *P_g            = 1.0 - *P_g;                                       // prob of being activated by gamma kill mode
		  S_HCP_component = gamma_parameters[i*4] * (1.0 - (*P_I) * (*P_g));  // activation prob, weighted by S0 for ith component
	  }else{
		  *P_I             = 1.0 - exp(-1.0 * (*sI_cm2) * fluence_cm2);     // prob of being activated by ion kill mode
		  S_HCP_component = gamma_parameters[i*4] * (*P_I);                // activation prob, weighted by S0 for ith component
	  }

	  /* Add to total HCP response and process next component */
	  *S_HCP += S_HCP_component;

  }

  /* Get gamma response for relative efficiency and gamma-kill dose*/
  AT_gamma_response(  n_tmp,
		  &total_dose_Gy,
		  gamma_model,
		  gamma_parameters,
		  false,
		  // return
		  S_gamma);

  *relative_efficiency      =       *S_HCP / *S_gamma;
  *gamma_dose_Gy            =       gamma_contribution * dose_contribution_Gy[0];

  /* Free allocated space, close output file and exit */
  free(norm_fluence);
  free(dose_contribution_Gy);
  free(params);

  if(write_output){
	  fclose(output_file);
  }
}


void AT_run_SPISS_method(  const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameters[],
    const long    er_model,
    const long    gamma_model,            // TODO do we really use gamma response here ?
    const double  gamma_parameters[],
    const long    n_runs,
    const long    N2,
    const double  fluence_factor,
    const int     write_output,
    const long    importance_sampling,
    double        results[])
{
  printf("\n############################################################\n");
  printf("\n############################################################\n");
  printf("This is AmTrack - SPISS algorithm\n");
  printf("\n");

  FILE*    output_file = NULL;
  if( write_output ){
    output_file    =  fopen("SPISS.log","w");
    if (output_file == NULL) return;                      // File error
  }

  // The histogram initialization and handling has been adapted to CPPSC
  // although some features are not used here
  long    n_bins_f1 = AT_n_bins_for_singe_impact_local_dose_distrib(          number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      N2);

  double* f1_parameters        = (double*)calloc(AT_SC_F1_PARAMETERS_SINGLE_LENGTH * number_of_field_components, sizeof(double));

  AT_RDD_f1_parameters_mixed_field( number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      f1_parameters
  );

  double*  f1_d_Gy                                      =  (double*)calloc(n_bins_f1, sizeof(double));
  double*  f1_dd_Gy                                     =  (double*)calloc(n_bins_f1, sizeof(double));
  double*  f1                                           =  (double*)calloc(n_bins_f1, sizeof(double));

  AT_single_impact_local_dose_distrib(           number_of_field_components,
      E_MeV_u,
      particle_no,
      fluence_cm2_or_dose_Gy,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      N2,
      n_bins_f1,
      f1_parameters,
      // from here: return values
      f1_d_Gy,
      f1_dd_Gy,
      f1);

  double*  fluence_cm2    =  (double*)calloc(number_of_field_components, sizeof(double));
  double*  dose_Gy        =  (double*)calloc(number_of_field_components, sizeof(double));

  long i;
  if(fluence_cm2_or_dose_Gy[0] < 0){
    for (i = 0; i < number_of_field_components; i++){
      dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    AT_fluence_cm2(  number_of_field_components,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        fluence_cm2);
  }else{
    for (i = 0; i < number_of_field_components; i++){
      fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
    }
    AT_dose_Gy_from_fluence_cm2(  number_of_field_components,
        E_MeV_u,
        particle_no,
        fluence_cm2,
        material_no,
        dose_Gy);
  }
  double*  norm_fluence                                 =  (double*)calloc(number_of_field_components, sizeof(double));

  // Normalize fluence vector
  AT_normalize(    number_of_field_components,
                fluence_cm2,
                norm_fluence);

  const double u  =       AT_total_u(     number_of_field_components,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      er_model);

  free( fluence_cm2 );
  free( dose_Gy );

  double*  accu_fluence                                 =  (double*)calloc(number_of_field_components, sizeof(double));

  // Get accumulated normalized fluence for later sampling of particle type
  accu_fluence[0]   =   norm_fluence[0];

  if(number_of_field_components > 1){
    for (i = 1; i < number_of_field_components; i++){
      accu_fluence[i] +=  accu_fluence[i-1] + norm_fluence[i];
    }
  }

  free(norm_fluence);

  long     n_bins_f;
  double   u_start;
  long     n_convolutions;

  AT_n_bins_for_low_fluence_local_dose_distribution(   u,
      fluence_factor,
      N2,
      n_bins_f1,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      // from here: return values
      &n_bins_f,
      &u_start,
      &n_convolutions);

  double*  f_d_Gy               =  (double*)calloc(n_bins_f, sizeof(double));
  double*  f_dd_Gy              =  (double*)calloc(n_bins_f, sizeof(double));
  double*  f                    =  (double*)calloc(n_bins_f, sizeof(double));
  double   f0                   =  0.0;

  AT_low_fluence_local_dose_distribution(  n_bins_f1,
      N2,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      n_bins_f,
      f_d_Gy,
      f_dd_Gy,
      f);

  // We are only interested in f_d_Gy and f_dd_Gy, so clear f
  for (i = 0; i < n_bins_f; i++){
    f[i] = 0;
  }

  if(importance_sampling){
    printf("\n");
    printf("Importance sampling chosen. Biasing function G(r)=r^%ld\n", importance_sampling);
  }else{
    printf("\n");
    printf("No importance sampling chosen.\n");
  }

  // init RNG
  gsl_rng * rng1   = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rng1, 12345678);

  long  act_number_particles;
  double  d_Gy;
  double  d_j_Gy;
  double  weight;
  double  r_m;
  long    n_tmp = 1;
  double  F;
  long    bin_no;
  double  max_bin_Gy    = log10(f_d_Gy[n_bins_f-1]);
  double  min_bin_Gy    = log10(f_d_Gy[0]);
  double  dd_bin_Gy     = (max_bin_Gy - min_bin_Gy) / (double)n_bins_f;

  // Do n_runs runs
  for (i = 0; i < n_runs; i++){
    // Get actual number particles for this run from Poisson generator
    act_number_particles          =   (long)gsl_ran_poisson(rng1, u);
    // Reset local dose for run i
    d_Gy                = 0.0;
    // Reset weight for importance sampling
    weight              = 1.0;
    // Add n individual doses according to their distribution
    long j;
    for (j = 0; j < act_number_particles; j++){
      // (1) draw random number 0..1 and sample particle type
      F                 = gsl_rng_uniform (rng1);
      long k;
      for (k = 0; k < number_of_field_components; k++){
        if (accu_fluence[k] >= F){
          break;
        }
      }

      // (2) draw again random number 0..1 for radius sampling
      F = gsl_rng_uniform (rng1);

      // (3) Apply importance sampling / weighting
      if (importance_sampling){
        weight  *= importance_sampling * pow(F, importance_sampling - 1.0);
        F        = pow(F, importance_sampling);
      }

      // (4) get dose d_Gy[j](r_max * F)
      r_m        = f1_parameters[k*9 + 2] * sqrt(F); // r_max for particle type k * 0..1
      AT_D_RDD_Gy(      n_tmp,
          &r_m,
          E_MeV_u[k],
          particle_no[k],
          material_no,
          rdd_model,
          rdd_parameters,
          er_model,
          &d_j_Gy);

      // (5) Add dose
      d_Gy += d_j_Gy;
    }

    // Fill dose into histogram
    if (d_Gy == 0.0){
      f0 += weight;
    }
    else{
      bin_no = floor((log10(d_Gy) - min_bin_Gy + 3.0*dd_bin_Gy/2.0) / dd_bin_Gy);
      if (bin_no > n_bins_f) bin_no = n_bins_f;
      f[bin_no - 1]             += weight / f_dd_Gy[bin_no - 1];
    }
    if(i%100 == 0){
      printf("Run %ld done.\n", i);
    }
  }

  // Normalize f
  double norm    = 0.0;
  double d_check = 0.0;
  for (i = 0; i < n_bins_f; i++){
    norm       += f_dd_Gy[i] * f[i];
  }
  norm += f0;
  for (i = 0; i < n_bins_f; i++){
    f[i]       /= norm;
    d_check    += f_d_Gy[i]*f_dd_Gy[i]*f[i];
  }

  if( write_output ){
    fprintf(output_file, "SPISS\n");
    fprintf(output_file, "number of runs: %ld\n",   n_runs);
    fprintf(output_file, "check D / Gy:   %4.3e\n", d_check);
    fprintf(output_file, "norm:           %4.3e\n", norm);
    fprintf(output_file, "f_n (%ld bins)\n", n_bins_f);
    fprintf(output_file, "f0: %4.2e\n", f0);
    for (i = 0; i < n_bins_f; i++){
      fprintf(output_file, "%ld; %4.2e; %4.2e; %4.2e\n", i+1, f_d_Gy[i], f_dd_Gy[i], f[i]);
    }
    fprintf(output_file, "f_1 (%ld bins)\n", n_bins_f1);
    for (i = 0; i < n_bins_f1; i++){
      fprintf(output_file, "%ld; %4.2e; %4.2e; %4.2e\n", i+1, f1_d_Gy[i], f1_dd_Gy[i], f1[i]);
    }
    fprintf(output_file, "\n");
    fprintf(output_file, "AmTrack SPISS run finished.\n");
    fprintf(output_file, "############################################################\n");
    fprintf(output_file, "############################################################\n");
    fclose(output_file);
  }

  free(accu_fluence);

  /* TODO memory might be not freed before !!!!
        free(f1_parameters);
        free(f1_d_Gy);
        free(f1_dd_Gy);
        free(f1);
        free(f_d_Gy);
//      free(f_dd_Gy);
        free(f);
   */

}

