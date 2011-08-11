/**
 * @brief Algorithms for ATMs based on Compound Poisson Processes (CPP)
 */

/*
 *    AT_Algorithms.c
 *    ===============
 *
 *    Created on: 18.11.2010
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

#include "AT_Algorithms_CPP.h"

void AT_run_CPPSC_method(  const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    stopping_power_source_no,
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

  /* Clear results */
  *relative_efficiency             = 0.0;
  *d_check                         = 0.0;
  *S_HCP                           = 0.0;
  *S_gamma                         = 0.0;
  *mean_number_of_tracks_contrib   = 0.0;
  *start_number_of_tracks_contrib  = 0.0;
  *n_convolutions                  = 0.0;
  *lower_Jensen_bound              = 0.0;
  *upper_Jensen_bound              = 0.0;

  /* Get array size for single impact dose
   * distribution for later memory allocation */
  long     n_bins_f1 = AT_n_bins_for_single_impact_local_dose_distrib(  number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      N2,
      stopping_power_source_no);

  /* Get f1 parameters - containing the most
   * relevant information on the tracks of
   * the mixed field, such as min/max local
   * dose, track radius etc. */
  double*  f1_parameters      =  (double*)calloc(AT_SC_F1_PARAMETERS_SINGLE_LENGTH * number_of_field_components, sizeof(double));
  AT_RDD_f1_parameters_mixed_field( number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      stopping_power_source_no,
      f1_parameters
  );

  /* Get local dose dictribution for the
   * impact of a single particle */
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
      stopping_power_source_no,
      f1_d_Gy,
      f1_dd_Gy,
      f1);

  /* Depending on user's input, convert
   * dose to fluence or vice versa */
  double*  fluence_cm2    =  (double*)calloc(number_of_field_components, sizeof(double));
  if(fluence_cm2_or_dose_Gy[0] < 0){
    double*  dose_Gy      =  (double*)calloc(number_of_field_components, sizeof(double));
    for (i = 0; i < number_of_field_components; i++){
      dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    AT_fluence_cm2_from_dose_Gy(  number_of_field_components,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        stopping_power_source_no,
        fluence_cm2);
    free( dose_Gy );
  }else{
    for (i = 0; i < number_of_field_components; i++){
      fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
    }
  }

  /* Compute the mean number of tracks that
   * deposit dose in a representative point
   * of the detector/cell */
  *mean_number_of_tracks_contrib  =       AT_mean_number_of_tracks_contrib(     number_of_field_components,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      er_model,
      stopping_power_source_no);

  free( fluence_cm2 );

  /* Get array size for low fluence local dose
   * distribution for later memory allocation */
  long      n_bins_f;
  AT_n_bins_for_low_fluence_local_dose_distribution(  *mean_number_of_tracks_contrib,
      fluence_factor,
      N2,
      n_bins_f1,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      &n_bins_f,
      start_number_of_tracks_contrib,
      n_convolutions);

  /* Get low fluence local dose distribution */
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

  /* Convolute this low fluence distribution
   * n_convolution times with itself to
   * get to the desired dose/fluence */
  AT_SuccessiveConvolutions(  *mean_number_of_tracks_contrib,
      n_bins_f,
      &N2,
      &n_bins_f1,
      f_d_Gy,
      f_dd_Gy,
      f,
      &f0,
      fdd,
      dfdd,
      d_check,
      write_output,
      shrink_tails,
      shrink_tails_under,
      adjust_N2);

  /* For the resulting local dose distribution
   * compute the response applying the gamma
   * response function first to each bin */
  long     n_bins_f_used  	= n_bins_f1;
  double*  S            	= (double*)calloc(n_bins_f_used, sizeof(double));
  AT_get_response_distribution_from_dose_distribution(  n_bins_f_used,
      f_d_Gy,
      f,
      gamma_model,
      gamma_parameters,
      lethal_events_mode,
      S);

  /* Then get the particle response by getting
   * the expected response value */
  *S_HCP = AT_get_ion_response_from_response_distribution( n_bins_f_used,
		  f_dd_Gy,
		  f,
		  S);

  /* Get the gamma response from the
   * expected dose value */
  *S_gamma = AT_get_gamma_response_for_average_dose( n_bins_f_used,
		  f_d_Gy,
		  f_dd_Gy,
		  f,
		  gamma_model,
		  gamma_parameters,
		  lethal_events_mode);

  /* Relative efficiency */
  *relative_efficiency = *S_HCP / *S_gamma;


  /* For Jensen's bounds:
   * get zero-dose reponse and its log */
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

  /* Free allocated memory and return*/
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

void AT_run_CPPSS_method(  const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    stopping_power_source_no,
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
  long    n_bins_f1 = AT_n_bins_for_single_impact_local_dose_distrib(          number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      N2,
      stopping_power_source_no);

  double* f1_parameters        = (double*)calloc(AT_SC_F1_PARAMETERS_SINGLE_LENGTH * number_of_field_components, sizeof(double));

  AT_RDD_f1_parameters_mixed_field( number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameters,
      er_model,
      stopping_power_source_no,
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
      stopping_power_source_no,
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
    AT_fluence_cm2_from_dose_Gy(  number_of_field_components,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        stopping_power_source_no,
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
        stopping_power_source_no,
        dose_Gy);
  }
  double*  norm_fluence                                 =  (double*)calloc(number_of_field_components, sizeof(double));

  // Normalize fluence vector
  AT_normalize(    number_of_field_components,
                fluence_cm2,
                norm_fluence);

  const double u  =       AT_mean_number_of_tracks_contrib(     number_of_field_components,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      er_model,
      stopping_power_source_no);

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
          stopping_power_source_no,
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
}
