/**
 * @file
 * @brief libamtrack main file holding the amorphous track routines for RE/RBE calculation
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

void AT_run_SPIFF_method(  const long  n,
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
    double        results[])
{

  long     n_bins_f1 = AT_n_bins_for_singe_impact_local_dose_distrib(  n,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      N2);

  double*  f1_parameters      =  (double*)calloc(AT_SC_F1_PARAMETERS_SINGLE_LENGTH * n, sizeof(double));

  AT_RDD_f1_parameters_mixed_field( n,
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

  AT_SC_get_f1(  n,
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

  double*  fluence_cm2    =  (double*)calloc(n, sizeof(double));

  long i;

  if(fluence_cm2_or_dose_Gy[0] < 0){
    double*  dose_Gy      =  (double*)calloc(n, sizeof(double));
    for (i = 0; i < n; i++){
      dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    // convert dose to fluence
    AT_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        fluence_cm2);
    free( dose_Gy );
  }else{
    for (i = 0; i < n; i++){
      fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
    }
  }

  const double u  =       AT_total_u(     n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      er_model);

  free( fluence_cm2 );

  long      n_bins_f;
  double    u_start;
  long      n_convolutions;

  AT_SC_get_f_array_size(  u,
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

  double*  f_d_Gy       =  (double*)calloc(n_bins_f, sizeof(double));
  double*  f_dd_Gy      =  (double*)calloc(n_bins_f, sizeof(double));
  double*  f            =  (double*)calloc(n_bins_f, sizeof(double));
  double*  fdd          =  (double*)calloc(n_bins_f, sizeof(double));
  double*  dfdd         =  (double*)calloc(n_bins_f, sizeof(double));
  double   f0           =  0.0;
  double   d_check      =  0.0;

  AT_SC_get_f_start(  n_bins_f1,
      N2,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      n_bins_f,
      f_d_Gy,
      f_dd_Gy,
      f);

  AT_SuccessiveConvolutions(  u,
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
      &d_check,
      write_output,
      shrink_tails,
      shrink_tails_under,
      adjust_N2);

  long     n_bins_f_used  = n_bins_f1;

  double*  S            =  (double*)calloc(n_bins_f_used, sizeof(double));
  double   S_HCP, S_gamma, efficiency;

  AT_get_response_distribution_from_dose_distribution(  n_bins_f_used,
      f_d_Gy,
      f,
      gamma_model,
      gamma_parameters,
      lethal_events_mode,
      S);

  S_HCP = AT_get_ion_response_from_response_distribution( n_bins_f_used,
		  f_dd_Gy,
		  f,
		  S);

  S_gamma = AT_get_gamma_response_for_average_dose( n_bins_f_used,
		  f_d_Gy,
		  f_dd_Gy,
		  f,
		  gamma_model,
		  gamma_parameters,
		  lethal_events_mode);

  efficiency = S_HCP / S_gamma;

//  AT_get_gamma_response(  n_bins_f_used,
//
//      f_d_Gy,
//      f_dd_Gy,
//
//      f,
//      f0,
//      gamma_model,
//      gamma_parameters,
//      lethal_events_mode,
//      // return
//
//      S,
//      &S_HCP,
//      &S_gamma,
//      &efficiency);

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
  double lower_Jensen_bound = f0*log_s0;
  for (i = 0; i < n_bins_f_used; i++){
	  double	log_s	= 0.0;
	  if (S[i] > 0){
		  log_s 	= log(S[i]);
	  }
	  lower_Jensen_bound 	+= f[i] * log_s * f_dd_Gy[i];
  }
  lower_Jensen_bound	=	exp(lower_Jensen_bound);

  /* Compute upper bound */
  double upper_Jensen_bound = f0*s0;
  for (i = 0; i < n_bins_f_used; i++){
	  upper_Jensen_bound 	+= f[i] * S[i] * f_dd_Gy[i];
  }

  results[0]      =  efficiency;        // 0 - 4: algo independent results
  results[1]      =  d_check;
  results[2]      =  S_HCP;
  results[3]      =  S_gamma;
  results[5]      =  u;                 // 5 - 9: algo specific: u
  results[6]      =  u_start;
  results[7]      =  n_convolutions;
  results[8]	  =  lower_Jensen_bound;
  results[9]	  =  upper_Jensen_bound;

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



void AT_run_IGK_method(  const long  n,
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
    double  results[])
{

  ////////////////////////////////////////////////////////////////////////////////////////////
  // 1. normalize fluence, get total fluence and dose

  // if fluence_cm2 < 0 the user gave doses in Gy rather than fluences, so in that case convert them first
  // only the first entry will be check
  long   i;
  double*  fluence_cm2    =  (double*)calloc(n, sizeof(double));
  double*  dose_Gy        =  (double*)calloc(n, sizeof(double));

  if(fluence_cm2_or_dose_Gy[0] < 0){
    for (i = 0; i < n; i++){
      dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    AT_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        fluence_cm2);
  }else{
    for (i = 0; i < n; i++){
      fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
    }
    AT_dose_Gy_from_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        fluence_cm2,
        material_no,
        dose_Gy);
  }

  double total_dose_Gy = 0.0;
  double total_fluence_cm2    = 0.0;

  for (i = 0; i < n; i++){
    total_dose_Gy      +=  dose_Gy[i];
    total_fluence_cm2  +=  fluence_cm2[i];
  }

  free( dose_Gy );

  double u_single;

  double*  norm_fluence          =  (double*)calloc(n, sizeof(double));
  double*  dose_contribution_Gy  =  (double*)calloc(n, sizeof(double));

  /* TODO: Replace by explicit functions */
  for (i = 0; i < n; i++){
        double LET_MeV_cm2_g = AT_LET_MeV_cm2_g_single(E_MeV_u[i], particle_no[i], material_no);
        double single_impact_fluence_cm2 = AT_single_impact_fluence_cm2_single(E_MeV_u[i], material_no, er_model);
    norm_fluence[i]          =  fluence_cm2[i] / total_fluence_cm2;
    u_single                 =  fluence_cm2[i] / single_impact_fluence_cm2;
    double single_impact_dose_Gy = AT_single_impact_dose_Gy_single(LET_MeV_cm2_g, single_impact_fluence_cm2);
    dose_contribution_Gy[i]  =  u_single * single_impact_dose_Gy;
  }

  free( fluence_cm2 );

  // Get accumulated normalized fluence for later sampling of particle type
  double*  accu_fluence          =  (double*)calloc(n, sizeof(double));
  accu_fluence[0]            =  norm_fluence[0];
  if(n > 1){
    for (i = 1; i < n; i++){
      accu_fluence[i] += accu_fluence[i-1] + norm_fluence[i];
    }
  }
  free(accu_fluence);
  // TODO do we really need accu_fluence ? it is not needed anywhere else

  //TODO rename KatseMitGlatse to something more reasonable
  FILE*    output_file = NULL;
  if( write_output ){
    output_file    =  fopen("KatseMitGlatse.log","w");
    if (output_file == NULL) return;                      // File error

    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "This is SGP efficiency Katz, version(2009/10/08).\n");
    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "\n\n\n");
  }

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

  long   n_tmp = 1;

  // Browse gamma parameters
  long   n_components         = 0;
  long   n_gamma_parameters   = 0;
  while  (gamma_parameters[n_gamma_parameters] != 0){
    n_gamma_parameters  += 4;
    n_components        += 1;
  }

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

  double   S_HCP                = 0.0;
  double   S_gamma              = 0.0;

  double   sI_m2                = 0.0;
  double   sI_cm2               = 0.0;
  double   P_I                  = 0.0;
  double   P_g                  = 0.0;
  double   gamma_contribution   = 0.0;
  double   cross_section_ratio  = 0.0;

  for(i = 0; i < n_components; i++){
    long j;
    for (j = 1; j < 4; j++){
      params->gamma_parameters[j] = gamma_parameters[i*4 + j];
    }
    // First: get activation cross section for ion mode
    gsl_set_error_handler_off();

    double error;
    gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (10000);
    gsl_function F;
    F.function           = &AT_sI_int;
    F.params             = (void*)params;
    double   lower_lim_m  = 0.0;
    if(rdd_model == RDD_KatzPoint){
      lower_lim_m = rdd_parameters[0];
    }
    double   upper_lim_m = AT_max_electron_range_m( *E_MeV_u, (int)material_no, (int)er_model);
    // TODO energy is an array of size n, why do we calculate upper_lim_m only from one energy value ?

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

    sI_m2  *= 2.0 * M_PI;
    sI_cm2  = sI_m2 * 10000.0;

    gsl_integration_workspace_free (w1);

    // TODO: INTERCEPT Katz point RDD for m / c detectors here!

    // Get saturation cross-section
    double   s0_m2 = 0.0;
    double   a0_m  = 0.0;
    if(rdd_model == RDD_KatzExtTarget){
      a0_m = rdd_parameters[1];
    }
    if(rdd_model == RDD_Geiss ||
        rdd_model == RDD_KatzSite){
      a0_m  =  rdd_parameters[0];
    }
    s0_m2   = saturation_cross_section_factor * M_PI * gsl_pow_2(a0_m);

    // Ion-kill probability
    double   fluence_cm2  = norm_fluence[0] * total_fluence_cm2;    // norm. fluence for particle i * total_fluence
    double   D_Gy         = dose_contribution_Gy[0];                // dose by particle i

    cross_section_ratio  = sI_m2 / s0_m2;
    double S_HCP_component;

    if( (cross_section_ratio < 1) & (cross_section_ratio >= 0) & (params->gamma_parameters[2] > 1)){
      P_I                = exp(-1.0 * sI_cm2 * fluence_cm2); // prob of being activated by ion kill mode
      gamma_contribution = 1.0 - cross_section_ratio;
      double   gamma_D_Gy = gamma_contribution * D_Gy;
      AT_gamma_response(  n_tmp,
          &gamma_D_Gy,
          gamma_model,
          params->gamma_parameters,
          false,
          // return
          &P_g);
      P_g             = 1.0 - P_g;                                  // prob of being activated by gamma kill mode
      S_HCP_component = gamma_parameters[i*4] * (1.0 - P_I * P_g);  // activation prob, weighted by S0 for ith component
    }else{
      P_I             = 1.0 - exp(-1.0 * sI_cm2 * fluence_cm2);     // prob of being activated by ion kill mode
      S_HCP_component = gamma_parameters[i*4] * P_I;                // activation prob, weighted by S0 for ith component
    }

    S_HCP += S_HCP_component;

  }

  AT_gamma_response(  n_tmp,
      &total_dose_Gy,
      gamma_model,
      gamma_parameters,
      false,
      // return
      &S_gamma);

  results[0]              =       S_HCP / S_gamma;
  results[1]              =       0.0;
  results[2]              =       S_HCP;
  results[3]              =       S_gamma;
  results[4]              =       0.0;
  results[5]              =       sI_cm2;
  results[6]              =       gamma_contribution * dose_contribution_Gy[0];
  results[7]              =       P_I;
  results[8]              =       P_g;
  results[9]              =       0.0;


  free(norm_fluence);
  free(dose_contribution_Gy);

  free(params);
  fclose(output_file);
}


void AT_run_SPISS_method(  const long  n,
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

  // The histogram initialization and handling has been adapted to SPIFF
  // although some features are not used here
  long    n_bins_f1 = AT_n_bins_for_singe_impact_local_dose_distrib(          n,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      N2);

  double* f1_parameters        = (double*)calloc(AT_SC_F1_PARAMETERS_SINGLE_LENGTH * n, sizeof(double));

  AT_RDD_f1_parameters_mixed_field( n,
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

  AT_SC_get_f1(           n,
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

  double*  fluence_cm2    =  (double*)calloc(n, sizeof(double));
  double*  dose_Gy        =  (double*)calloc(n, sizeof(double));

  long i;
  if(fluence_cm2_or_dose_Gy[0] < 0){
    for (i = 0; i < n; i++){
      dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    AT_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        fluence_cm2);
  }else{
    for (i = 0; i < n; i++){
      fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
    }
    AT_dose_Gy_from_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        fluence_cm2,
        material_no,
        dose_Gy);
  }
  double*  norm_fluence                                 =  (double*)calloc(n, sizeof(double));

  // Normalize fluence vector
  AT_normalize(    n,
                fluence_cm2,
                norm_fluence);

  const double u  =       AT_total_u(     n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      er_model);

  free( fluence_cm2 );
  free( dose_Gy );

  double*  accu_fluence                                 =  (double*)calloc(n, sizeof(double));

  // Get accumulated normalized fluence for later sampling of particle type
  accu_fluence[0]   =   norm_fluence[0];

  if(n > 1){
    for (i = 1; i < n; i++){
      accu_fluence[i] +=  accu_fluence[i-1] + norm_fluence[i];
    }
  }

  free(norm_fluence);

  long     n_bins_f;
  double   u_start;
  long     n_convolutions;

  AT_SC_get_f_array_size(   u,
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

  AT_SC_get_f_start(  n_bins_f1,
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
      for (k = 0; k < n; k++){
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

void AT_GSM_calculate_multiple_dose_histograms( const long  number_of_field_components,
    const double   	E_MeV_u[],
    const double   	fluence_cm2[],
    const long     	particle_no[],
    const long     	material_no,
    const long     	rdd_model,
    const double   	rdd_parameter[],
    const long     	er_model,
    const long     	nX,
    const double   	pixel_size_m,
    const long		N_runs,
    const long		N_repetitions,
    const long     	number_of_bins,
    const double   	dose_bin_centers_Gy[],
    double    		dose_bin_width_Gy[],
    double *       	mean_d_check_Gy,
    double *       	sd_d_check_Gy,
    double *       	mean_zero_dose_fraction,
    double *       	sd_zero_dose_fraction,
    double        	mean_dose_frequency_Gy[],
    double        	sd_dose_frequency_Gy[]){

	long i,j,k;

	AT_histoOld_get_bin_widths(	number_of_bins,
								dose_bin_centers_Gy,
								dose_bin_width_Gy);

	*mean_d_check_Gy	= 0.0;
	*sd_d_check_Gy	= 0.0;

	for (i = 0; i < number_of_bins; i++){
		mean_dose_frequency_Gy[i]	= 0.0;
		sd_dose_frequency_Gy[i]		= 0.0;
	}

	double	zero_dose_fraction, zero_dose_fraction_run;
	double*	dose_frequency_Gy		= (double*)calloc(number_of_bins, sizeof(double));
	double*	dose_frequency_Gy_run	= (double*)calloc(number_of_bins, sizeof(double));

	/* Create and initialize random number generator */
	gsl_rng * rng  								= gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rng, 137);
	unsigned long  random_number_generator_seed	= gsl_rng_get(rng);


	for (i = 0; i < N_repetitions; i++){
		zero_dose_fraction		= 0.0;
		for (j = 0; j < number_of_bins; j++){
			dose_frequency_Gy[j]		= 0.0;
		}

		for (j = 0; j < N_runs; j++){
			zero_dose_fraction_run		= 0.0;
			for (k = 0; k < number_of_bins; k++){
				dose_frequency_Gy_run[k]		= 0.0;
			}

			AT_GSM_calculate_dose_histogram( number_of_field_components,
					E_MeV_u,
					fluence_cm2,
					particle_no,
					material_no,
					rdd_model,
					rdd_parameter,
					er_model,
					nX,
					pixel_size_m,
					number_of_bins,
					dose_bin_centers_Gy,
					&random_number_generator_seed,
					&zero_dose_fraction_run,
					dose_frequency_Gy_run);

			zero_dose_fraction		+= zero_dose_fraction_run;

			for (k = 0; k < number_of_bins; k++){
				dose_frequency_Gy[k] += dose_frequency_Gy_run[k];
			}
		}

		zero_dose_fraction			/= (double)N_runs;

		for(j = 0 ; j < number_of_bins; j++){
			dose_frequency_Gy[j]     		/=  (double)N_runs;
		}


		/* compute <d> */
		float cur_d_check_Gy	=	0.0;
		for (j = 0; j < number_of_bins; j++){
			cur_d_check_Gy		+=	dose_bin_centers_Gy[j] * dose_frequency_Gy[j]; // * dose_bin_width_Gy[j];
		}

		*mean_d_check_Gy		+=	cur_d_check_Gy;
		*sd_d_check_Gy			+=	cur_d_check_Gy * cur_d_check_Gy;
	}

	/* Effective calculation of running mean and stdev of x in n runs: */
	/* (1) add x and x^2 in every run (--> sum_x and sum_x2)           */
	/* (2) mean  = sum_x / n                                           */
	/* (3) stdev = sqrt((sum_x2 - sum_x * sum_x)/(n-1))                */

	*mean_d_check_Gy		/= 	(double)N_repetitions;

	if(N_repetitions > 1){
		*sd_d_check_Gy		= sqrt(*sd_d_check_Gy/((double)N_repetitions) - (*mean_d_check_Gy)*(*mean_d_check_Gy));
	}else{
		*sd_d_check_Gy		= 0.0f;
	}

	free(dose_frequency_Gy);
	free(dose_frequency_Gy_run);
}
