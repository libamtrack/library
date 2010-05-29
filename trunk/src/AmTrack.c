/**
 * @file
 * @brief libamtrack main file holding the amorphous track routines for RE/RBE calculation
 */

/*
 *    AmTrack.c
 *    =========
 *
 *    Created on: 28.07.2009
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

#include "AmTrack.h"

void AT_PrintName(void){
  printf("This is libamtrack.\n");
}


int AT_GetNumber(void){ // TODO to be removed
  return 137;
}

void AT_run_SPIFF_method(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
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
  FILE* output_file = NULL;
  if(write_output){
    output_file    =  fopen("SuccessiveConvolutions.log","w");
    if (output_file == NULL) return;                      // File error

    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "This is LGC2.2 core - successive convolution mode (2008/08/12).\n");
  }

  long     n_bins_f1 = AT_SC_get_f1_array_size(  n,
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
      fluence_cm2,
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

  double*  fluence_cm2_local    =  (double*)calloc(n, sizeof(double));

  long i;

  if(fluence_cm2[0] < 0){
    double*  dose_Gy_local      =  (double*)calloc(n, sizeof(double));
    for (i = 0; i < n; i++){
      dose_Gy_local[i] = -1.0 * fluence_cm2[i];
    }
    // convert dose to fluence
    AT_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        dose_Gy_local,
        material_no,
        fluence_cm2_local);
    free( dose_Gy_local );
  }else{
    for (i = 0; i < n; i++){
      fluence_cm2_local[i] = fluence_cm2[i];
    }
  }

  const double u  =       AT_total_u(     n,
      E_MeV_u,
      particle_no,
      fluence_cm2_local,
      material_no,
      er_model);

  free( fluence_cm2_local );

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

  AT_get_gamma_response(  n_bins_f_used,

      f_d_Gy,
      f_dd_Gy,

      f,
      f0,
      gamma_model,
      gamma_parameters,
      lethal_events_mode,
      // return

      S,
      &S_HCP,
      &S_gamma,
      &efficiency);

  results[0]      =  efficiency;        // 0 - 4: algo independent results
  results[1]      =  d_check;
  results[2]      =  S_HCP;
  results[3]      =  S_gamma;
  results[5]      =  u;                 // 5 - 9: algo specific: u
  results[6]      =  u_start;
  results[7]      =  n_convolutions;

  //////////////////////////////////////////
  // Output results
  //////////////////////////////////////////
  if(write_output){
    fprintf(output_file, "\n\nResults:\n");
    fprintf(output_file, "\ndose check / Gy:         %4.3e Gy", results[1]);
    fprintf(output_file, "\nmean impact number u:    %4.3e Gy", results[5]);
    fprintf(output_file, "\nstart impact number:     %4.3e Gy", results[6]);
    fprintf(output_file, "\nnumber of convolutions:  %ld",      (long)(results[7]));
    fprintf(output_file, "\n\nf0: %4.3e\n", f0);
    for (i = 0; i < n_bins_f_used;i++){
      fprintf(output_file, "%ld; %4.2e; %4.2e; %4.2e\n",     i+1, f_d_Gy[i], f_dd_Gy[i], f[i]);
    }
    fprintf(output_file, "\nDone.\n###############################################\n");
    fprintf(output_file, "###############################################\n");
    fclose(output_file);
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


void AT_GSM_shoot_particles_on_grid(  const long  number_of_field_components,
		const double   fluence_cm2[],
		const double   sample_grid_size_m,
		long           number_of_particles_in_field_component[],
		double*        x_position[],
		double*        y_position[]){

	/* Create and initialize random number generator */
	gsl_rng * rng  =  gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rng, 12345678);

	/* Calculate total fluence */
	double total_fluence_cm2     =  AT_sum(  number_of_field_components, fluence_cm2);

	/* Normalize fluence vector */
	double* norm_fluence         =  (double*)calloc(number_of_field_components, sizeof(double));
	AT_normalize(  number_of_field_components,  fluence_cm2, norm_fluence);

	/* Mean number of particles from one component in the sample grid area */
	double mean_number_particles_in_field_component;

	long i;
	double sample_grid_area_cm2 = gsl_pow_2(100.0 * sample_grid_size_m);
	for (i = 0; i < number_of_field_components; i++){
		mean_number_particles_in_field_component   =  sample_grid_area_cm2 * total_fluence_cm2 * norm_fluence[i];

		/* Actual number of particles from one component in the sample grid area, taken from Poissonian distribution */
		number_of_particles_in_field_component[i]  =  (long)gsl_ran_poisson(rng, mean_number_particles_in_field_component);

		/* Now we know how many particles from one component we have, so position tables can be allocated */
		x_position[i] =  (double*)calloc(number_of_particles_in_field_component[i], sizeof(double));
		y_position[i] =  (double*)calloc(number_of_particles_in_field_component[i], sizeof(double));
	}
	free(norm_fluence);

	/* Sample particle positions for every component from uniform distribution */
	long j;
	for (i = 0; i < number_of_field_components; i++){
		for( j = 0 ; j < number_of_particles_in_field_component[i] ; j ++){
			x_position[i][j] = gsl_rng_uniform_pos(rng) * sample_grid_size_m;
			y_position[i][j] = gsl_rng_uniform_pos(rng) * sample_grid_size_m;
		}
	}
	gsl_rng_free(rng);
}



void AT_run_GSM_method(  const long  n,
    const double   E_MeV_u[],
    const long     particle_no[],
    const double   fluence_cm2[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameters[],
    const long     er_model,
    const long     gamma_model,
    const double   gamma_parameters[],
    const long     N_runs,
    const long     N2,
    const double   fluence_factor,
    const bool     write_output,
    const long     nX,
    const double   voxel_size_m,
    const bool     lethal_events_mode,
    double         results[])
{

	/** in case of debugging
	printf("n = %ld\n", n);
	printf("E_MeV_u[] = %g, %g, %g\n", E_MeV_u[0], E_MeV_u[1], E_MeV_u[2]);
	printf("particle_no[] = %ld, %ld, %ld\n", particle_no[0], particle_no[1], particle_no[2]);
	printf("fluence_cm2[] = %g, %g, %g\n", fluence_cm2[0], fluence_cm2[1], fluence_cm2[2]);
	printf("material_no = %ld\n", material_no);
	printf("rdd_model = %ld\n", rdd_model);
	printf("rdd_parameters[] = %g\n", rdd_parameters[0]);
	printf("er_model = %ld\n", er_model);
	printf("gamma_model = %ld\n", gamma_model);
	printf("gamma_parameters[] = %g, %g, %g, %g, %g\n", gamma_parameters[0], gamma_parameters[1], gamma_parameters[2], gamma_parameters[3], gamma_parameters[4]);
	printf("N_runs = %ld\n", N_runs);
	printf("N2 = %ld\n", N2);
	printf("fluence_factor = %g\n", fluence_factor);
	printf("write_output = %d\n", write_output);
	printf("nX = %ld\n", nX);
	printf("voxel_size_m = %g\n", voxel_size_m);
	printf("lethal_events_mode = %d\n", lethal_events_mode);
*/

  long    i, j, k, m;
  long    n_grid                =  gsl_pow_2(nX);
  double  calc_grid_size_m      =  voxel_size_m * nX;
  double  calc_grid_area_cm2    =  gsl_pow_2(calc_grid_size_m) * 10000;

  FILE*     output_file = NULL;
  if( write_output ){
    output_file    =  fopen("GridSummation.log","w");
    if (output_file == NULL) return;                      // File error

    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "This is grid summation algorithm, version(2009/06/30).\n");
    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "\n\n\n");
    fprintf(output_file, "calc grid:      %ld*%ld = %ld pixels\n", nX, nX, n_grid);
    fprintf(output_file, "calc grid size/m:     %e\n", calc_grid_size_m);
    fprintf(output_file, "calc grid area/cm2:  %e\n", calc_grid_area_cm2);
  }

  // Alloc checkerboard arrays
  double*  grid_d_Gy  =  (double*)calloc(n_grid, sizeof(double));
  double*  grid_S     =  (double*)calloc(n_grid, sizeof(double));

  // Clear results
  for (i = 0; i < 10; i++){
    results[i] = 0.0;
  }

  // Get f1, f parameters
  double  max_r_max_m          =  0.0;
//  long    n_bins_f1            =  AT_SC_get_f1_array_size( n,
//      E_MeV_u,
//      particle_no,
//      material_no,
//      rdd_model,
//      rdd_parameters,
//      er_model,
//      N2);

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

  for (i = 0; i < n; i++){
    max_r_max_m    =   GSL_MAX(max_r_max_m, f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 2]);
  }

  if( write_output ){
    fprintf(output_file, "f1 parameters for %ld particles\n", n);
    for (i = 0; i < n; i++){
      fprintf(output_file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",   f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 0],
          f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 1],
          f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 2],
          f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3],
          f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4],
          f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 5],
          f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 6],
          f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 7],
          f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 8]);
    }
    fprintf(output_file, "\nOverall r.max/m = %e\n\n",   max_r_max_m);
  }

//  double* f1_d_Gy              =  (double*)calloc(n_bins_f1, sizeof(double));
//  double* f1_dd_Gy             =  (double*)calloc(n_bins_f1, sizeof(double));
//  double* f1                   =  (double*)calloc(n_bins_f1, sizeof(double));
//
//  AT_SC_get_f1(  n,
//      E_MeV_u,
//      particle_no,
//      fluence_cm2,
//      material_no,
//      rdd_model,
//      rdd_parameters,
//      er_model,
//      N2,
//      n_bins_f1,
//      f1_parameters,
//      f1_d_Gy,
//      f1_dd_Gy,
//      f1);

  double*  fluence_cm2_local    =  (double*)calloc(n, sizeof(double));

  if(fluence_cm2[0] < 0){
    double*  dose_Gy_local        =  (double*)calloc(n, sizeof(double));
    for (i = 0; i < n; i++){
      dose_Gy_local[i] = -1.0 * fluence_cm2[i];
    }
    // dose to fluence
    AT_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        dose_Gy_local,
        material_no,
        fluence_cm2_local);
    free( dose_Gy_local );
  }else{
    for (i = 0; i < n; i++){
      fluence_cm2_local[i] = fluence_cm2[i];
    }
  }

  double* norm_fluence         =  (double*)calloc(n, sizeof(double));

  // Normalize fluence vector
  AT_normalize(    n,
                fluence_cm2_local,
                norm_fluence);

//  const double u  =       AT_total_u(     n,
//      E_MeV_u,
//      particle_no,
//      fluence_cm2_local,
//      material_no,
//      er_model);

  double total_fluence_cm2  =       AT_sum(   n,
      fluence_cm2_local);

  free( fluence_cm2_local );

  if( write_output ){
    for (i = 0; i < n; i++){
      fprintf(output_file, "Particle %ld: relative fluence: %4.2e\n", i+1, norm_fluence[i]);
    }
  }

//  long	  n_bins_f;
//  double  u_start;
//  long    n_convolutions;
//
//  AT_SC_get_f_array_size(   u,
//      fluence_factor,
//      N2,
//      n_bins_f1,
//      f1_d_Gy,
//      f1_dd_Gy,
//      f1,
//      // from here: return values
//      &n_bins_f,
//      &u_start,
//      &n_convolutions);
//
//  double*  f_d_Gy        =  (double*)calloc(n_bins_f, sizeof(double));
//  double*  f_dd_Gy       =  (double*)calloc(n_bins_f, sizeof(double));
//  double*  f             =  (double*)calloc(n_bins_f, sizeof(double));
//  double   f0            =  0.0;
//
//  AT_SC_get_f_start(  n_bins_f1,
//      N2,
//      f1_d_Gy,
//      f1_dd_Gy,
//      f1,
//      n_bins_f,
//      f_d_Gy,
//      f_dd_Gy,
//      f);

  // We are only interested in f_d_Gy and f_dd_Gy, so clear f
//  for (i = 0; i < n_bins_f; i++){
//    f[i] = 0;
//  }
//
//  double max_bin_Gy	= log10(f_d_Gy[n_bins_f-1]);
//  double min_bin_Gy	= log10(f_d_Gy[0]);
//  double dd_bin_Gy	= (max_bin_Gy - min_bin_Gy) / n_bins_f;

  // Largest r.max --> calculate size of sample area
  double sample_grid_size_m    = calc_grid_size_m + 2.01 * max_r_max_m;
  double sample_grid_area_cm2  = sample_grid_size_m * sample_grid_size_m * 10000;
  if( write_output ){
    fprintf(output_file, "sample grid size/m   = %e\n", sample_grid_size_m);
    fprintf(output_file, "sample grid area/cm2 = %e\n", sample_grid_area_cm2);
  }

  // mean and actual number of particles on sample_area
  double* mean_number_particles =  (double*)calloc(n, sizeof(double));
  long*  act_number_particles   =  (long*)calloc(n, sizeof(double));
  for (i = 0; i < n; i++){
    mean_number_particles[i] = sample_grid_area_cm2 * total_fluence_cm2 * norm_fluence[i];        // Area * Total_fluence (particle i)
  }

  free(norm_fluence);

  // create and initialize RNGs
  gsl_rng * rng1  =  gsl_rng_alloc(gsl_rng_taus);
  gsl_rng * rng2  =  gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng1, 12345678);
  gsl_rng_set(rng2, 2345678);

  // run loop
  long n_particles;
  for (m = 0; m < N_runs; m++){
    double run_results[10];

    // sample particles numbers
    n_particles = 0;
    for (i = 0; i < n; i++){
      act_number_particles[i] =  (long)gsl_ran_poisson(rng1, mean_number_particles[i]);
      n_particles             += act_number_particles[i];
    }

    if( write_output ){
      if(N_runs <= 20){
        fprintf(output_file, "\n\nRun %ld:\n", m + 1);
        fprintf(output_file, "Actual number of particles (mean)\n");
        for (i = 0; i < n; i++){
          fprintf(output_file, "particle %ld: %ld (%e)\n", i, act_number_particles[i], mean_number_particles[i]);
        }
        fprintf(output_file, "\nIn total: %ld\n", n_particles);
      }
    }

    // alloc particle array
    double*  x_pos          = (double*)calloc(n_particles, sizeof(double));
    double*  y_pos          = (double*)calloc(n_particles, sizeof(double));
    long*    particle_index = (long*)calloc(n_particles, sizeof(long));
    double*  r_m            = (double*)calloc(n_particles, sizeof(double));
    double*  r_max_m        = (double*)calloc(n_particles, sizeof(double));
    long     n_tmp          = 1;
    double   d_tmp_Gy       = 0.0;

    // fill in index / r_max_m
    j = 0;
    k = 0;
    for (i = 0; i < n_particles; i++){
      if(k >= act_number_particles[j]){
        k = 0;
        j++;
      }
      k++;
      particle_index[i] = j;
      r_max_m[i]        = f1_parameters[j*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 2];
    }

//    if( write_output ){
//      fprintf(output_file, "particle.no; particle.index; r.max.m\n");
//      for (i=0; i < n_particles; i++){
//        fprintf(output_file, "%d; %d; %e\n", i, particle_index[i], r_max_m[i]);
//      }
//    }

    // sample particle positions
    for (i = 0; i < n_particles; i++){
      x_pos[i] = gsl_rng_uniform_pos(rng2) * sample_grid_size_m;
      y_pos[i] = gsl_rng_uniform_pos(rng2) * sample_grid_size_m;
    }

    //TODO include this speedup trick also for n > 1

    if( n == 1){

      // calculate total number of contributing parts
      long no_contr_part = 0;
      double cur_x_pos = 0;
      double cur_y_pos = 0;
      for (j = 0; j < nX; j++){                                      // y
        cur_y_pos    = max_r_max_m + ((double)j + 0.5f)*voxel_size_m;
        for (i = 0; i < nX; i++){                              // x
          cur_x_pos  = max_r_max_m + ((double)i + 0.5f)*voxel_size_m;
          for (k = 0; k < n_particles; k++){      // particles
            r_m[k]   = sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
            if(r_m[k] <= r_max_m[k]){               // does particle contribute?
              no_contr_part++;
            }
          }
        }
      }

      // allocate memory for vector of distances between track cores and all grid points
      double*  distances = (double*)calloc(no_contr_part, sizeof(double));

      // calculate distances between track cores and all grid points and save them into allocated table
      long distances_index = 0;
      for (j = 0; j < nX; j++){                                      // y
        double cur_y_pos       = max_r_max_m + ((double)j + 0.5f)*voxel_size_m;
        for (i = 0; i < nX; i++){                              // x
          double cur_x_pos     = max_r_max_m + ((double)i + 0.5f)*voxel_size_m;
          grid_d_Gy[j * nX + i] = 0.0;
          for (k = 0; k < n_particles; k++){      // particles
            r_m[k]            = sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));

            if(r_m[k] <= r_max_m[k]){               // does particle contribute?
              distances[distances_index] = r_m[k];
              distances_index++;
            }
          }// particle loop
        } // x loop
      } // y loop

      // just a crosscheck, distances_index should go up to no_contr_part
      long dinstances_index_max = distances_index;

      // allocate memory for doses values at given distances stored in _distances_ table
      double*  doses = (double*)calloc(dinstances_index_max, sizeof(double));

      // calculate doses at given distances
      AT_D_RDD_Gy( dinstances_index_max,
          distances,
          E_MeV_u[particle_index[0]],
          particle_no[particle_index[0]],
          material_no,
          rdd_model,
          rdd_parameters,
          er_model,
          doses);

      free(distances);

      // fill grid using doses values stored in doses table
      distances_index = 0;
      for (j = 0; j < nX; j++){                                      // y
        double cur_y_pos   = max_r_max_m + ((double)j + 0.5f)*voxel_size_m;
        for (i = 0; i < nX; i++){                              // x
          double cur_x_pos = max_r_max_m + ((double)i + 0.5f)*voxel_size_m;
          for (k = 0; k < n_particles; k++){      // particles
            r_m[k]        = sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
            if(r_m[k] <= r_max_m[k]){               // does particle contribute?
              grid_d_Gy[j * nX + i] += doses[distances_index];
              distances_index++;
            } // particle contribution
          }// particle loop
        } // x loop
      } // y loop

      free(doses);

    } else {

      // allocate memory for vector of vectors distances between track cores and all grid points
      double**  distances = (double**)calloc(n, sizeof(double*));
      long* no_contr_part = (long*)calloc(n, sizeof(long));
      long* distances_index = (long*)calloc(n, sizeof(long));

      for (i = 0; i < n; i++){
    	  no_contr_part[i] = 0;
    	  distances_index[i] = 0;
      }

      // grid loop
      for (j = 0; j < nX; j++){          // y
        double cur_y_pos   =  max_r_max_m + ((double)j + 0.5f)*voxel_size_m;
        for (i = 0; i < nX; i++){        // x
          double cur_x_pos =  max_r_max_m + ((double)i + 0.5f)*voxel_size_m;
          grid_d_Gy[j * nX + i]=  0.0;
          for (k = 0; k < n_particles; k++){  // particles
            r_m[k]        =  sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
            if(r_m[k] <= r_max_m[k]){    // does particle contribute?
              no_contr_part[particle_index[k]]++;
            } // particle contribution
          }// particle loop
        } // x loop
      } // y loop

      for (i = 0; i < n; i++){
    	  distances[i] = (double*)calloc( no_contr_part[i], sizeof(double) );
      }

      // calculate distances between track cores and all grid points and save them into allocated table
      for (j = 0; j < nX; j++){                                      // y
        double cur_y_pos       = max_r_max_m + ((double)j + 0.5f)*voxel_size_m;
        for (i = 0; i < nX; i++){                              // x
          double cur_x_pos     = max_r_max_m + ((double)i + 0.5f)*voxel_size_m;
          grid_d_Gy[j * nX + i] = 0.0;
          for (k = 0; k < n_particles; k++){      // particles
            r_m[k]            = sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));

            if(r_m[k] <= r_max_m[k]){               // does particle contribute?
              distances[particle_index[k]][distances_index[particle_index[k]]] = r_m[k];
              distances_index[particle_index[k]]++;
            }
          }// particle loop
        } // x loop
      } // y loop

      double**  doses = (double**)calloc(n, sizeof(double*));
      for (i = 0; i < n; i++){
    	  doses[i] = (double*)calloc( no_contr_part[i], sizeof(double) );
          // calculate doses at given distances
          AT_D_RDD_Gy( distances_index[i],
              distances[i],
              E_MeV_u[i],
              particle_no[i],
              material_no,
              rdd_model,
              rdd_parameters,
              er_model,
              doses[i]);
      }
      for (i = 0; i < n; i++){
    	  free(distances[i]);
      }

      for (i = 0; i < n; i++){
    	  distances_index[i] = 0;
      }
      free(distances);
      free(no_contr_part);

      // grid loop
      for (j = 0; j < nX; j++){          // y
        double cur_y_pos   =  max_r_max_m + ((double)j + 0.5f)*voxel_size_m;
        for (i = 0; i < nX; i++){        // x
          double cur_x_pos =  max_r_max_m + ((double)i + 0.5f)*voxel_size_m;
          grid_d_Gy[j * nX + i]=  0.0;
          for (k = 0; k < n_particles; k++){  // particles
            r_m[k]        =  sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
            if(r_m[k] <= r_max_m[k]){    // does particle contribute?
                grid_d_Gy[j * nX + i] += doses[particle_index[k]][distances_index[particle_index[k]]];
                distances_index[particle_index[k]]++;
            } // particle contribution
          }// particle loop
        } // x loop
      } // y loop

      for (i = 0; i < n; i++){
    	  free(doses[i]);
      }

      free(doses);
      free(distances_index);

      // grid loop - OLD APPROACH
//      for (j = 0; j < nX; j++){          // y
//        double cur_y_pos   =  max_r_max_m + ((double)j + 0.5f)*voxel_size_m;
//        for (i = 0; i < nX; i++){        // x
//          double cur_x_pos =  max_r_max_m + ((double)i + 0.5f)*voxel_size_m;
//          grid_d_Gy[j * nX + i]=  0.0;
//          for (k = 0; k < n_particles; k++){  // particles
//            r_m[k]        =  sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
//            if(r_m[k] <= r_max_m[k]){    // does particle contribute?
//              AT_D_RDD_Gy(  n_tmp,
//                  &r_m[k],
//                  E_MeV_u[particle_index[k]],
//                  particle_no[particle_index[k]],
//                  material_no,
//                  rdd_model,
//                  rdd_parameters,
//                  er_model,
//                  &d_tmp_Gy);
//              grid_d_Gy[j * nX + i]  +=  d_tmp_Gy;
//            } // particle contribution
//          }// particle loop
//        } // x loop
//      } // y loop

    }

    free(r_m);

    // get gamma response for local dose

    double d_total_Gy = 0.0;
    double S_HCP      = 0.0;

    if( lethal_events_mode ){
      if (gamma_model != GR_LinQuad){
        printf("##############################################################################\n");
        printf("Sorry, no Grid Summation model with other than Linear Quadratic gamma response\n");
        printf("Please choose gamma_model = %d. Exiting now...\n", GR_LinQuad);
        printf("##############################################################################\n");
        return;
      }

      double alpha  =  gamma_parameters[0];
      double beta   =  gamma_parameters[1];
      double D0     =  gamma_parameters[2];

      // averaging over number of lethal events
      for (i = 0; i < n_grid; i++){
        d_total_Gy += grid_d_Gy[i];
        if( grid_d_Gy[i] < D0 ){
          grid_S[i] = alpha * grid_d_Gy[i] + beta * gsl_pow_2(grid_d_Gy[i]);
        } else {
          grid_S[i] = alpha * D0 - beta * gsl_pow_2(D0) + ( alpha + 2 * beta * D0) * (grid_d_Gy[i] - D0);
        }
        S_HCP += grid_S[i];
      }
    } else {
      // get gamma response for local dose
      AT_gamma_response( n_grid,
          grid_d_Gy,
          gamma_model,
          gamma_parameters,
          grid_S);

      // averaging over the dose
      for (i = 0; i < n_grid; i++){
        d_total_Gy += grid_d_Gy[i];
        S_HCP      += grid_S[i];
      }
    }

    S_HCP       /=  n_grid;
    d_total_Gy  /= n_grid;

    if( lethal_events_mode ){
      S_HCP  = exp( - S_HCP );
    }

    double S_gamma  = 0.0;
    AT_gamma_response(  n_tmp,
        &d_total_Gy,
        gamma_model,
        gamma_parameters,
        &S_gamma);

    double efficiency  = 0.0;
    if(S_gamma > 0){
      efficiency = S_HCP / S_gamma;
    }

//    long bin_no;
//    // Fill dose into histogram
//    for(i = 0; i < n_grid; i++){
//      if (grid_d_Gy[i] == 0.0){
//        f0 += 1.0f;
//      }
//      else{
//        bin_no = floor((log10(grid_d_Gy[i]) - min_bin_Gy + 3.0*dd_bin_Gy/2.0) / dd_bin_Gy);
//        if (bin_no > n_bins_f) bin_no = n_bins_f;
//        f[bin_no - 1]		+= 1.0 / f_dd_Gy[bin_no - 1];
//      }
//    }

    // write graph (first run)
//    bool write_graph = false;
//    if(write_graph && (m == 0)){
//      FILE*    graph_file;
//      graph_file    =  fopen("GridGraph.csv","w");
//      if (graph_file == NULL) return;    // File error
//
//      fprintf(graph_file, "x.m; y.m; d.Gy; S\n");
//
//      for (j = 0; j < nX; j++){
//        for (i = 0; i < nX; i++){
//          fprintf(graph_file, "%e; %e; %e; %e\n",  max_r_max_m + ((double)i + 0.5)*voxel_size_m,
//              max_r_max_m + ((double)j + 0.5)*voxel_size_m,
//              grid_d_Gy[j * nX + i],
//              grid_S[j * nX + i]);
//        }
//      }
//      fclose(graph_file);
//
//    }

    run_results[0]    = efficiency;
    run_results[1]    = d_total_Gy;
    run_results[2]    = S_HCP;
    run_results[3]    = S_gamma;
    run_results[4]    = n_particles;

    // copy to results
    results[0]      += run_results[0];
    results[1]      += run_results[1];
    results[2]      += run_results[2];
    results[3]      += run_results[3];
    results[4]      += run_results[4];

    results[5]      += run_results[0]*run_results[0];
    results[6]      += run_results[1]*run_results[1];
    results[7]      += run_results[2]*run_results[2];
    results[8]      += run_results[3]*run_results[3];
    results[9]      += n_particles * n_particles;

    if( write_output ){
      if (N_runs <= 20){
        fprintf(output_file, "\n\nRun %ld results\n", m + 1);
        fprintf(output_file, "efficiency     = %e\n", run_results[0]);
        fprintf(output_file, "d.check.Gy     = %e\n", run_results[1]);
        fprintf(output_file, "S (HCP)        = %e\n", run_results[2]);
        fprintf(output_file, "S (gamma)      = %e\n", run_results[3]);
        fprintf(output_file, "no. particles  = %e\n", run_results[4]);
      }
    }

    free(x_pos);
    free(y_pos);
    free(particle_index);
    free(r_max_m);
  }// end run loop

  free(f1_parameters);

  gsl_rng_free(rng1);
  gsl_rng_free(rng2);

  free(mean_number_particles);
  free(act_number_particles);

  free(grid_d_Gy);
  free(grid_S);

  results[0]  /= N_runs;
  results[1]  /= N_runs;
  results[2]  /= N_runs;
  results[3]  /= N_runs;
  results[4]  /= N_runs;

  results[5]  /= N_runs;
  results[6]  /= N_runs;
  results[7]  /= N_runs;
  results[8]  /= N_runs;
  results[9]  /= N_runs;

  results[5]  -= results[0]*results[0];
  results[6]  -= results[1]*results[1];
  results[7]  -= results[2]*results[2];
  results[8]  -= results[3]*results[3];
  results[9]  -= results[4]*results[4];

  results[5]  = GSL_MAX(0, results[5]);
  results[6]  = GSL_MAX(0, results[6]);
  results[7]  = GSL_MAX(0, results[7]);
  results[8]  = GSL_MAX(0, results[8]);
  results[9]  = GSL_MAX(0, results[9]);

  if( N_runs > 1 ){
	  results[5]  = sqrt(results[5] / (N_runs - 1.));
	  results[6]  = sqrt(results[6] / (N_runs - 1.));
	  results[7]  = sqrt(results[7] / (N_runs - 1.));
	  results[8]  = sqrt(results[8] / (N_runs - 1.));
	  results[9]  = sqrt(results[9] / (N_runs - 1.));
  }

  if( write_output ){
    fprintf(output_file, "\n###############################################\nResults\n");
    fprintf(output_file, "efficiency     = %e +/- %e\n", results[0], results[5]);
    fprintf(output_file, "d.check.Gy     = %e +/- %e\n", results[1], results[6]);
    fprintf(output_file, "S (HCP)        = %e +/- %e\n", results[2], results[7]);
    fprintf(output_file, "S (gamma)      = %e +/- %e\n", results[3], results[8]);
    fprintf(output_file, "no. particles  = %e +/- %e\n", results[4], results[9]);
    fprintf(output_file, "###############################################\n");
  }

//  // Normalize f
//  double norm 		= 0.0;
//  double d_check	= 0.0;
//  for (i = 0; i < n_bins_f; i++){
//    norm	+= f_dd_Gy[i] * f[i];
//  }
//  norm += f0;
//  for (i = 0; i < n_bins_f; i++){
//    f[i]	/= norm;
//    d_check += f_d_Gy[i]*f_dd_Gy[i]*f[i];
//  }
//
//  if( write_output ){
//    fprintf(output_file, "Grid summation dose distribution\n");
//    fprintf(output_file, "check D / Gy:   %4.3e\n", d_check);
//    fprintf(output_file, "norm:           %4.3e\n", norm);
//    fprintf(output_file, "f_n (%ld bins)\n", n_bins_f);
//    fprintf(output_file, "f0: %4.2e\n", f0);
//    for (i = 0; i < n_bins_f; i++){
//      fprintf(output_file, "%ld; %4.2e; %4.2e; %4.2e\n", i+1, f_d_Gy[i], f_dd_Gy[i], f[i]);
//    }
//    fprintf(output_file, "f_1 (%ld bins)\n", n_bins_f1);
//    for (i = 0; i < n_bins_f1; i++){
//      fprintf(output_file, "%ld; %4.2e; %4.2e; %4.2e\n", i+1, f1_d_Gy[i], f1_dd_Gy[i], f1[i]);
//    }
//    fclose(output_file);
//  }
//
//  free(f_d_Gy);
//  free(f_dd_Gy);
//  free(f);
//
//  free(f1_d_Gy);
//  free(f1_dd_Gy);
//  free(f1);

}


void AT_run_IGK_method(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
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
  double*  fluence_cm2_local    =  (double*)calloc(n, sizeof(double));
  double*  dose_Gy_local        =  (double*)calloc(n, sizeof(double));

  if(fluence_cm2[0] < 0){
    for (i = 0; i < n; i++){
      dose_Gy_local[i] = -1.0 * fluence_cm2[i];
    }
    AT_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        dose_Gy_local,
        material_no,
        fluence_cm2_local);
  }else{
    for (i = 0; i < n; i++){
      fluence_cm2_local[i] = fluence_cm2[i];
    }
    AT_D_Gy(  n,
        E_MeV_u,
        particle_no,
        fluence_cm2_local,
        material_no,
        dose_Gy_local);
  }

  double total_dose_Gy = 0.0;
  double total_fluence_cm2    = 0.0;

  for (i = 0; i < n; i++){
    total_dose_Gy      +=  dose_Gy_local[i];
    total_fluence_cm2  +=  fluence_cm2_local[i];
  }

  free( dose_Gy_local );

  double u_single;

  double*  norm_fluence          =  (double*)calloc(n, sizeof(double));
  double*  dose_contribution_Gy  =  (double*)calloc(n, sizeof(double));

  for (i = 0; i < n; i++){
	double LET_MeV_cm2_g = AT_LET_MeV_cm2_g_single(E_MeV_u[i], particle_no[i], material_no);
	double single_impact_fluence_cm2 = AT_single_impact_fluence_cm2_single(E_MeV_u[i], material_no, er_model);
    norm_fluence[i]          =  fluence_cm2_local[i] / total_fluence_cm2;
    u_single                 =  fluence_cm2_local[i] / single_impact_fluence_cm2;
    double single_impact_dose_Gy = AT_single_impact_dose_Gy_single(LET_MeV_cm2_g, single_impact_fluence_cm2);
    dose_contribution_Gy[i]  =  u_single * single_impact_dose_Gy;
  }

  free( fluence_cm2_local );

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
    const double  fluence_cm2[],
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
  long    n_bins_f1 = AT_SC_get_f1_array_size(  	n,
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

  double*  f1_d_Gy          				=  (double*)calloc(n_bins_f1, sizeof(double));
  double*  f1_dd_Gy        				=  (double*)calloc(n_bins_f1, sizeof(double));
  double*  f1            				=  (double*)calloc(n_bins_f1, sizeof(double));

  AT_SC_get_f1(  	  n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
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

  double*  fluence_cm2_local    =  (double*)calloc(n, sizeof(double));
  double*  dose_Gy_local        =  (double*)calloc(n, sizeof(double));

  long i;
  if(fluence_cm2[0] < 0){
    for (i = 0; i < n; i++){
      dose_Gy_local[i] = -1.0 * fluence_cm2[i];
    }
    AT_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        dose_Gy_local,
        material_no,
        fluence_cm2_local);
  }else{
    for (i = 0; i < n; i++){
      fluence_cm2_local[i] = fluence_cm2[i];
    }
    AT_D_Gy(  n,
        E_MeV_u,
        particle_no,
        fluence_cm2_local,
        material_no,
        dose_Gy_local);
  }
  double*  norm_fluence                                 =  (double*)calloc(n, sizeof(double));

  // Normalize fluence vector
  AT_normalize(    n,
                fluence_cm2_local,
                norm_fluence);

  const double u  =       AT_total_u(     n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      er_model);

  free( fluence_cm2_local );
  free( dose_Gy_local );

  double*  accu_fluence                                 =  (double*)calloc(n, sizeof(double));

  // Get accumulated normalized fluence for later sampling of particle type
  accu_fluence[0]   =	norm_fluence[0];

  if(n > 1){
    for (i = 1; i < n; i++){
      accu_fluence[i] +=  accu_fluence[i-1] + norm_fluence[i];
    }
  }

  free(norm_fluence);

  long	   n_bins_f;
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

  double*  f_d_Gy         	=  (double*)calloc(n_bins_f, sizeof(double));
  double*  f_dd_Gy        	=  (double*)calloc(n_bins_f, sizeof(double));
  double*  f            	=  (double*)calloc(n_bins_f, sizeof(double));
  double   f0            	=  0.0;

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

  long	act_number_particles;
  double  d_Gy;
  double  d_j_Gy;
  double  weight;
  double  r_m;
  long	  n_tmp = 1;
  double  F;
  long	  bin_no;
  double  max_bin_Gy	= log10(f_d_Gy[n_bins_f-1]);
  double  min_bin_Gy	= log10(f_d_Gy[0]);
  double  dd_bin_Gy	= (max_bin_Gy - min_bin_Gy) / (double)n_bins_f;

  // Do n_runs runs
  for (i = 0; i < n_runs; i++){
    // Get actual number particles for this run from Poisson generator
    act_number_particles	  =   (long)gsl_ran_poisson(rng1, u);
    // Reset local dose for run i
    d_Gy		= 0.0;
    // Reset weight for importance sampling
    weight 		= 1.0;
    // Add n individual doses according to their distribution
    long j;
    for (j = 0; j < act_number_particles; j++){
      // (1) draw random number 0..1 and sample particle type
      F 		= gsl_rng_uniform (rng1);
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
      AT_D_RDD_Gy(  	n_tmp,
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
      f[bin_no - 1]		+= weight / f_dd_Gy[bin_no - 1];
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
