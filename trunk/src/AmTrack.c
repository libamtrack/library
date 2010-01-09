/**
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

void AT_interparticleDistance_m( const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  fluence_cm2,
    float*  results_m
){
  long i;
  float fluence;
  for( i = 0 ; i < *n ; i++ ){
    if( fluence_cm2[i] > 0 ){
      results_m[i] = 2.0f / sqrt(M_PI*1e4*fluence_cm2[i]);
    } else {
      fluence = (-fluence_cm2[i]) / (LET_MeV_cm2_g[i] * MeV_g_to_J_kg);
      results_m[i] = 2.0f / sqrt(M_PI*1e4*fluence);
    }
  }
}

void AT_inv_interparticleDistance_Gy( const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  distance_m,
    float*  results_Gy
){
  long i;
  float fluence;
  for( i = 0 ; i < *n ; i++ ){
    fluence = (2.0f/distance_m[i])*(2.0f/distance_m[i])*M_1_PI*1e-4;
    results_Gy[i] = fluence * (LET_MeV_cm2_g[i] * MeV_g_to_J_kg);
  }
}

void AT_inv_interparticleDistance_cm2( const long*   n,
    const float*  distance_m,
    float*  results_cm2
){
  long i;
  for( i = 0 ; i < *n ; i++ ){
    results_cm2[i] = (2.0f/distance_m[i])*(2.0f/distance_m[i])*M_1_PI*1e-4;
  }
}


void AT_SPIFF(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const float*  fluence_cm2,
    const long*  material_no,
    const long*  RDD_model,
    const float*  RDD_parameters,
    const long*  ER_model,
    const float*  ER_parameters,
    const long*  gamma_model,
    const float*  gamma_parameters,
    const long*  N2,
    const float*  fluence_factor,
    const bool*  write_output,
    const bool*  shrink_tails,
    const float*  shrink_tails_under,
    const bool*  adjust_N2,
    const bool*   lethal_events_mode,
    float*  results)
{

  long  n_bins_f1;
  float*  f1_parameters      =  (float*)calloc(9 * (*n), sizeof(float));

  AT_SC_get_f1_array_size(  n,

      E_MeV_u,
      particle_no,
      material_no,
      RDD_model,
      RDD_parameters,
      ER_model,
      ER_parameters,

      N2,
      &n_bins_f1,
      f1_parameters);

  float*  f_parameters      =  (float*)calloc(7, sizeof(float));

  float*  norm_fluence      =  (float*)calloc(*n, sizeof(float));
  float*  dose_contribution_Gy  =  (float*)calloc(*n, sizeof(float));

  float*  f1_d_Gy          =  (float*)calloc(n_bins_f1, sizeof(float));
  float*  f1_dd_Gy        =  (float*)calloc(n_bins_f1, sizeof(float));
  float*  f1            =  (float*)calloc(n_bins_f1, sizeof(float));

  AT_SC_get_f1(  n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      /* detector parameters */
      material_no,
      RDD_model,
      RDD_parameters,
      /* electron range model */
      ER_model,
      ER_parameters,
      /* algorithm parameters*/
      N2,
      &n_bins_f1,
      /* f1 parameters*/
      f1_parameters,
      // from here: return values
      norm_fluence,
      dose_contribution_Gy,
      f_parameters,
      f1_d_Gy,
      f1_dd_Gy,

      f1);

  long      n_bins_f;
  float      u_start;
  long      n_convolutions;


  AT_SC_get_f_array_size(  &f_parameters[0],      // = u
      fluence_factor,
      N2,
      &n_bins_f1,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      // from here: return values
      &n_bins_f,
      &u_start,
      &n_convolutions);

  float*  f_d_Gy          =  (float*)calloc(n_bins_f, sizeof(float));
  float*  f_dd_Gy          =  (float*)calloc(n_bins_f, sizeof(float));
  float*  f            =  (float*)calloc(n_bins_f, sizeof(float));
  float*  fdd            =  (float*)calloc(n_bins_f, sizeof(float));
  float*  dfdd          =  (float*)calloc(n_bins_f, sizeof(float));
  float  f0            =  0.0f;
  float  d_check          =  0.0f;

  AT_SC_get_f_start(  &u_start,
      &n_bins_f1,
      N2,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      &n_bins_f,
      // from here: return values
      f_d_Gy,
      f_dd_Gy,
      f);

  AT_SuccessiveConvolutions(  &f_parameters[0],    // u
      &n_bins_f,
      N2,
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
      &write_output,
      &shrink_tails,
      shrink_tails_under,
      &adjust_N2);

  long      n_bins_f_used  = n_bins_f1;

  float*  S            =  (float*)calloc(n_bins_f_used, sizeof(float));
  float  S_HCP, S_gamma, efficiency;

  AT_get_gamma_response(  &n_bins_f_used,

      f_d_Gy,
      f_dd_Gy,


      f,
      &f0,
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
  results[5]      =  f_parameters[0];    // 5 - 9: algo specific: u
  results[6]      =  u_start;
  results[7]      =  n_convolutions;

  free(f1_parameters);
  free(f_parameters);
  free(norm_fluence);
  free(dose_contribution_Gy);
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

void AT_GSM(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const float*  fluence_cm2,
    const long*  material_no,
    const long*  RDD_model,
    const float*  RDD_parameters,
    const long*  ER_model,
    const float*  ER_parameters,
    const long*  gamma_model,
    const float*  gamma_parameters,
    const long*  N_runs,
    const long*   N2,
    const float*  fluence_factor,
    const bool*   write_output,
    const long*   nX,
    const float*  grid_size_m,
    const bool*   lethal_events_mode,
    float*  results)
{
  FILE*    output_file;
  struct tm  *start_tm, *end_tm;
  time_t     start_t, end_t;
  start_t    = time(NULL);

  long     i, j, k, m;
  long    n_grid          = (*nX) * (*nX);
  float    calc_grid_size_m    = (*grid_size_m) * (*nX);
  float    calc_grid_area_cm2    = calc_grid_size_m * calc_grid_size_m * 10000;

  output_file    =  fopen("GridSummation.log","w");
  if (output_file == NULL) return;                      // File error

  fprintf(output_file, "##############################################################\n");
  fprintf(output_file, "##############################################################\n");
  fprintf(output_file, "This is grid summation algorithm, version(2009/06/30).\n");
  fprintf(output_file, "##############################################################\n");
  fprintf(output_file, "\n\n\n");
  start_tm     = localtime(&start_t);
  fprintf(output_file, "Start time and date: %s\n", asctime(start_tm));
  fprintf(output_file, "calc grid:      %ld*%ld = %ld pixels\n", *nX, *nX, n_grid);
  fprintf(output_file, "calc grid size/m:     %e\n", calc_grid_size_m);
  fprintf(output_file, "calc grid area/cm2:  %e\n", calc_grid_area_cm2);

  // Alloc checkerboard arrays
  float*  grid_d_Gy = (float*)calloc(n_grid, sizeof(float));
  float*  grid_S    = (float*)calloc(n_grid, sizeof(float));

  // Clear results
  for (i = 0; i < 10; i++){
    results[i] = 0.0;
  }

  // Get f1, f parameters
  float* f1_parameters        = (float*)calloc((*n)*9, sizeof(float));
  float  max_r_max_m          = 0.0f;
  long   n_bins_f1			= 0;


  AT_SC_get_f1_array_size(  	n,
      E_MeV_u,
      particle_no,
      material_no,
      RDD_model,
      RDD_parameters,
      ER_model,
      ER_parameters,
      N2,
      &n_bins_f1,
      f1_parameters);

  fprintf(output_file, "f1 parameters for %ld particles\n", *n);
  for (i = 0; i < *n; i++){
    fprintf(output_file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",   f1_parameters[i*9 + 0],
        f1_parameters[i*9 + 1],
        f1_parameters[i*9 + 2],
        f1_parameters[i*9 + 3],
        f1_parameters[i*9 + 4],
        f1_parameters[i*9 + 5],
        f1_parameters[i*9 + 6],
        f1_parameters[i*9 + 7],
        f1_parameters[i*9 + 8]);
    max_r_max_m        =   fmaxf(max_r_max_m, f1_parameters[i*9 + 2]);
  }

  fprintf(output_file, "\nOverall r.max/m = %e\n\n",   max_r_max_m);

  float* f_parameters         = (float*)calloc(7, sizeof(float));
  float* f1_d_Gy          	=  (float*)calloc(n_bins_f1, sizeof(float));
  float* f1_dd_Gy        	    =  (float*)calloc(n_bins_f1, sizeof(float));
  float* f1            	  	=  (float*)calloc(n_bins_f1, sizeof(float));
  float* norm_fluence         = (float*)calloc(*n, sizeof(float));
  float* dose_contribution_Gy = (float*)calloc(*n, sizeof(float));

  AT_SC_get_f1(  n,            // for f parameters only
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      RDD_model,
      RDD_parameters,
      ER_model,
      ER_parameters,

      N2,
      &n_bins_f1,
      f1_parameters,
      norm_fluence,
      dose_contribution_Gy,
      f_parameters,
      f1_d_Gy,
      f1_dd_Gy,
      f1);

  long	n_bins_f;
  float   u_start;
  long    n_convolutions;
  float	u				= f_parameters[0];

  AT_SC_get_f_array_size(   &u,
      fluence_factor,
      N2,
      &n_bins_f1,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      // from here: return values
      &n_bins_f,
      &u_start,
      &n_convolutions);

  float*  f_d_Gy         	=  (float*)calloc(n_bins_f, sizeof(float));
  float*  f_dd_Gy        	=  (float*)calloc(n_bins_f, sizeof(float));
  float*  f            	=  (float*)calloc(n_bins_f, sizeof(float));
//  float*  fdd            	=  (float*)calloc(n_bins_f, sizeof(float)); // UNUSED
//  float*  dfdd          	=  (float*)calloc(n_bins_f, sizeof(float)); // UNUSED
  float   f0            	=  0.0f;

  AT_SC_get_f_start(  &u_start,
      &n_bins_f1,
      N2,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      &n_bins_f,
      // from here: return values
      f_d_Gy,
      f_dd_Gy,
      f);

  // We are only interested in f_d_Gy and f_dd_Gy, so clear f
  for (i = 0; i < n_bins_f; i++){
    f[i] = 0;
  }

  float	max_bin_Gy	= log10(f_d_Gy[n_bins_f-1]);
  float	min_bin_Gy	= log10(f_d_Gy[0]);
  float	dd_bin_Gy	= (max_bin_Gy - min_bin_Gy) / n_bins_f;

  fprintf(output_file, "f parameters\n");
  fprintf(output_file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n",   f_parameters[0],
      f_parameters[1],
      f_parameters[2],
      f_parameters[3],
      f_parameters[4],
      f_parameters[5],
      f_parameters[6]);

  // Largest r.max --> calculate size of sample area
  float sample_grid_size_m  = calc_grid_size_m + 2.01f * max_r_max_m;
  float sample_grid_area_cm2  = sample_grid_size_m * sample_grid_size_m * 10000;
  fprintf(output_file, "sample grid size/m   = %e\n", sample_grid_size_m);
  fprintf(output_file, "sample grid area/cm2 = %e\n", sample_grid_area_cm2);

  // mean and actual number of particles on sample_area
  float* mean_number_particles = (float*)calloc(*n, sizeof(float));
  long*  act_number_particles  =  (long*)calloc(*n, sizeof(float));
  for (i = 0; i < *n; i++){
    mean_number_particles[i] = sample_grid_area_cm2 * f_parameters[1] * norm_fluence[i];        // Area * Total_fluence (particle i)
  }

  // create and initialize RNGs
  gsl_rng * rng1 = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng * rng2 = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng1, 12345678);
  gsl_rng_set(rng2, 2345678);


  // run loop
  long n_particles;
  for (m = 0; m < *N_runs; m++){
    float* run_results = (float*)calloc(10, sizeof(float));

    // sample particles numbers
    n_particles = 0;
    for (i = 0; i < *n; i++){
      act_number_particles[i] =  (long)gsl_ran_poisson(rng1, mean_number_particles[i]);
      n_particles             += act_number_particles[i];
    }

    if(*N_runs <= 20){
      fprintf(output_file, "\n\nRun %ld:\n", m + 1);
      fprintf(output_file, "Actual number of particles (mean)\n");
      for (i = 0; i < *n; i++){
        fprintf(output_file, "particle %ld: %ld (%e)\n", i, act_number_particles[i], mean_number_particles[i]);
      }
      fprintf(output_file, "\nIn total: %ld\n", n_particles);
    }

    // alloc particle array
    float*  x_pos          = (float*)calloc(n_particles, sizeof(float));
    float*  y_pos          = (float*)calloc(n_particles, sizeof(float));
    long*   particle_index =  (long*)calloc(n_particles, sizeof(long));
    float*  r_m            = (float*)calloc(n_particles, sizeof(float));
    float*  r_max_m        = (float*)calloc(n_particles, sizeof(float));
    long    n_tmp          = 1;
    float   d_tmp_Gy       = 0.0;

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
      r_max_m[i]        = f1_parameters[j*9 + 2];
    }

    //    fprintf(output_file, "particle.no; particle.index; r.max.m\n");
    //    for (i=0; i < n_particles; i++){
    //      fprintf(output_file, "%d; %d; %e\n", i, particle_index[i], r_max_m[i]);
    //    }

    // sample particle positions
    for (i = 0; i < n_particles; i++){
      x_pos[i] = (float)gsl_rng_uniform_pos(rng2) * sample_grid_size_m;
      y_pos[i] = (float)gsl_rng_uniform_pos(rng2) * sample_grid_size_m;
    }

    // TODO include this speedup trick also for n > 1

    if( n_particles == 1){

      // calculate total number of contributing parts
      long no_contr_part = 0;
      float cur_x_pos = 0;
      float cur_y_pos = 0;
      for (j = 0; j < *nX; j++){                                      // y
        cur_y_pos    = max_r_max_m + ((float)j + 0.5f)*(*grid_size_m);
        for (i = 0; i < *nX; i++){                              // x
          cur_x_pos  = max_r_max_m + ((float)i + 0.5f)*(*grid_size_m);
          for (k = 0; k < n_particles; k++){      // particles
            r_m[k]   = sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
            if(r_m[k] <= r_max_m[k]){               // does particle contribute?
              no_contr_part++;
            }
          }
        }
      }

      //printf("\nContributing items: %ld\n", no_contr_part);
      // allocate memory for vector of distances between track cores and all grid points
      float*  distances = (float*)calloc(no_contr_part, sizeof(float));

      // calculate distances between track cores and all grid points and save them into allocated table
      long distances_index = 0;
      for (j = 0; j < *nX; j++){                                      // y
        float cur_y_pos       = max_r_max_m + ((float)j + 0.5f)*(*grid_size_m);
        for (i = 0; i < *nX; i++){                              // x
          float cur_x_pos     = max_r_max_m + ((float)i + 0.5f)*(*grid_size_m);
          grid_d_Gy[j * (*nX) + i] = 0.0f;
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
      float*  doses = (float*)calloc(dinstances_index_max, sizeof(float));

      // calculate doses at given distances
      AT_D_RDD_Gy( &dinstances_index_max,

          distances,
          &E_MeV_u[particle_index[0]],
          &particle_no[particle_index[0]],
          material_no,

          RDD_model,
          RDD_parameters,

          ER_model,
          ER_parameters,

          doses);

      free(distances);

      // fill grid using doses values stored in doses table
      distances_index = 0;
      for (j = 0; j < *nX; j++){                                      // y
        float cur_y_pos   = max_r_max_m + ((float)j + 0.5f)*(*grid_size_m);
        for (i = 0; i < *nX; i++){                              // x
          float cur_x_pos = max_r_max_m + ((float)i + 0.5f)*(*grid_size_m);
          for (k = 0; k < n_particles; k++){      // particles
            r_m[k]        = sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
            if(r_m[k] <= r_max_m[k]){               // does particle contribute?
              grid_d_Gy[j * (*nX) + i] += doses[distances_index];
              distances_index++;
            } // particle contribution
          }// particle loop
        } // x loop
      } // y loop

      free(doses);


    } else {
      // grid loop
      for (j = 0; j < *nX; j++){          // y
        float cur_y_pos   =  max_r_max_m + ((float)j + 0.5f)*(*grid_size_m);
        for (i = 0; i < *nX; i++){        // x
          float cur_x_pos =  max_r_max_m + ((float)i + 0.5f)*(*grid_size_m);
          grid_d_Gy[j * (*nX) + i]=  0.0f;
          for (k = 0; k < n_particles; k++){  // particles
            r_m[k]        =  sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) + (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
            if(r_m[k] <= r_max_m[k]){    // does particle contribute?
              AT_D_RDD_Gy(  &n_tmp,
                  &r_m[k],
                  &E_MeV_u[particle_index[k]],
                  &particle_no[particle_index[k]],
                  material_no,
                  RDD_model,
                  RDD_parameters,
                  ER_model,
                  ER_parameters,
                  &d_tmp_Gy);
              grid_d_Gy[j * (*nX) + i]  +=  d_tmp_Gy;
            } // particle contribution
          }// particle loop
        } // x loop
      } // y loop

    }

    // get gamma response for local dose

    float d_total_Gy = 0.0f;
    float S_HCP      = 0.0f;

    if( *lethal_events_mode ){
      if (*gamma_model != GR_LinQuad){
        printf("##############################################################################\n");
        printf("Sorry, no Grid Summation model with other than Linear Quadratic gamma response\n");
        printf("Please choose gamma_model = %d. Exiting now...\n", GR_LinQuad);
        printf("##############################################################################\n");
        return;
      }

      float alpha = gamma_parameters[0];
      float beta = gamma_parameters[1];
      float D0 = gamma_parameters[2];

      // averaging over number of lethal events
      for (i = 0; i < n_grid; i++){
        d_total_Gy += grid_d_Gy[i];
        if( grid_d_Gy[i] < D0 ){

          grid_S[i] = alpha * grid_d_Gy[i] + beta * grid_d_Gy[i] * grid_d_Gy[i];
        } else {
          grid_S[i] = alpha * D0 - beta * D0 * D0 + ( alpha + 2 * beta * D0) * (grid_d_Gy[i] - D0);
        }
        S_HCP += grid_S[i];
      }
    } else {
      // get gamma response for local dose
      AT_gamma_response( &n_grid,
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

    if( *lethal_events_mode ){
      S_HCP  = expf( - S_HCP );
    }

    float S_gamma  = 0.0f;
    AT_gamma_response(  &n_tmp,
        &d_total_Gy,
        gamma_model,
        gamma_parameters,
        &S_gamma);

    float efficiency  = 0.0f;
    if(S_gamma > 0){
      efficiency = S_HCP / S_gamma;
    }

    long bin_no;
    // Fill dose into histogram
    for(i = 0; i < n_grid; i++){
      if (grid_d_Gy[i] == 0.0f){
        f0 += 1.0f;
      }
      else{
        bin_no = floor((log10(grid_d_Gy[i]) - min_bin_Gy + 3*dd_bin_Gy/2) / dd_bin_Gy);
        if (bin_no > n_bins_f) bin_no = n_bins_f;
        f[bin_no - 1]		+= 1.0f / f_dd_Gy[bin_no - 1];
      }
    }

    // write graph (first run)
    bool write_graph = true;
    if(write_graph & (m == 0)){
      FILE*    graph_file;
      graph_file    =  fopen("GridGraph.csv","w");
      if (graph_file == NULL) return;    // File error

      fprintf(graph_file, "x.m; y.m; d.Gy; S\n");

      for (j = 0; j < *nX; j++){
        for (i = 0; i < *nX; i++){
          fprintf(graph_file, "%e; %e; %e; %e\n",  max_r_max_m + ((float)i + 0.5f)*(*grid_size_m),
              max_r_max_m + ((float)j + 0.5f)*(*grid_size_m),
              grid_d_Gy[j * (*nX) + i],
              grid_S[j * (*nX) + i]);
        }
      }
    }

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

    if (*N_runs <= 20){
      fprintf(output_file, "\n\nRun %ld results\n", m + 1);
      fprintf(output_file, "efficiency     = %e\n", run_results[0]);
      fprintf(output_file, "d.check.Gy     = %e\n", run_results[1]);
      fprintf(output_file, "S (HCP)       = %e\n", run_results[2]);
      fprintf(output_file, "S (gamma)     = %e\n", run_results[3]);
      fprintf(output_file, "no. particles    = %e\n", run_results[4]);
    }

    free(x_pos);
    free(y_pos);
    free(particle_index);
    free(r_max_m);
  }// end run loop

  results[0]  /= *N_runs;
  results[1]  /= *N_runs;
  results[2]  /= *N_runs;
  results[3]  /= *N_runs;
  results[4]  /= *N_runs;

  results[5]  /= *N_runs;
  results[6]  /= *N_runs;
  results[7]  /= *N_runs;
  results[8]  /= *N_runs;
  results[9]  /= *N_runs;

  results[5]  -= results[0]*results[0];
  results[6]  -= results[1]*results[1];
  results[7]  -= results[2]*results[2];
  results[8]  -= results[3]*results[3];
  results[9]  -= results[4]*results[4];

  results[5]  = fmaxf(0, results[5]);
  results[6]  = fmaxf(0, results[6]);
  results[7]  = fmaxf(0, results[7]);
  results[8]  = fmaxf(0, results[8]);
  results[9]  = fmaxf(0, results[9]);

  results[5]  = sqrt(results[5] / (*N_runs - 1));
  results[6]  = sqrt(results[6] / (*N_runs - 1));
  results[7]  = sqrt(results[7] / (*N_runs - 1));
  results[8]  = sqrt(results[8] / (*N_runs - 1));
  results[9]  = sqrt(results[9] / (*N_runs - 1));

  fprintf(output_file, "\n###############################################\nResults\n");
  fprintf(output_file, "efficiency     = %e +/- %e\n", results[0], results[5]);
  fprintf(output_file, "d.check.Gy     = %e +/- %e\n", results[1], results[6]);
  fprintf(output_file, "S (HCP)       = %e +/- %e\n", results[2], results[7]);
  fprintf(output_file, "S (gamma)     = %e +/- %e\n", results[3], results[8]);
  fprintf(output_file, "no. particles    = %e +/- %e\n", results[4], results[9]);
  fprintf(output_file, "###############################################\n");
  end_t  = time(NULL);
  end_tm = localtime(&end_t);
  float timespan_s  = difftime(end_t, start_t);
  fprintf(output_file, "\nEnd time and date: %s\n", asctime(end_tm));
  fprintf(output_file, "\nTime per run [s]:      %4.2e\n", timespan_s / *N_runs);
  fprintf(output_file, "\nTime per pixel [s]:    %4.2e\n", timespan_s / (*N_runs * n_grid));
  fprintf(output_file, "###############################################\n");
  fprintf(output_file, "###############################################\n");

  // Normalize f
  float norm 		= 0.0f;
  float d_check	= 0.0f;
  for (i = 0; i < n_bins_f; i++){
    norm	+= f_dd_Gy[i] * f[i];
  }
  norm += f0;
  for (i = 0; i < n_bins_f; i++){
    f[i]	/= norm;
    d_check += f_d_Gy[i]*f_dd_Gy[i]*f[i];
  }

  // Write f, f1
  fprintf(output_file, "Grid summation dose distribution\n");
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


  gsl_rng_free(rng1);
  gsl_rng_free(rng2);

  free(f1_parameters);
  free(f_parameters);
  free(norm_fluence);
  free(dose_contribution_Gy);

  free(mean_number_particles);
  free(act_number_particles);

  free(grid_d_Gy);
  free(grid_S);

  fclose(output_file);
}


void AT_IGK(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const float*  fluence_cm2,
    const long*  material_no,
    const long*  RDD_model,
    const float*  RDD_parameters,
    const long*  ER_model,
    const float*  ER_parameters,
    const long*  gamma_model,
    const float*  gamma_parameters,
    float*  results)
{
  FILE*    output_file;
  struct tm  *start_tm;
//  struct tm *end_tm; // NOT USED
  time_t start_t;
// time_t end_t; // NOT USED
  start_t    = time(NULL);

  output_file    =  fopen("KatseMitGlatse.log","w");
  if (output_file == NULL) return;                      // File error

  fprintf(output_file, "##############################################################\n");
  fprintf(output_file, "##############################################################\n");
  fprintf(output_file, "This is SGP efficiency Katz, version(2009/07/13).\n");
  fprintf(output_file, "##############################################################\n");
  fprintf(output_file, "\n\n\n");
  start_tm     = localtime(&start_t);
  fprintf(output_file, "Start time and date: %s\n", asctime(start_tm));

  if (*gamma_model != 1){
    fprintf(output_file, "##############################################################\n");
    fprintf(output_file, "Sorry, no Katz with other than the general hit-target model\n");
    fprintf(output_file, "Please choose gamma_model = 1. Exiting now...\n");
    fprintf(output_file, "##############################################################\n");
    return;
  }

  long   i;
  // Browse gamma parameters and pick one-hit components
  long   n_components     = 0;
  long  n_gamma_parameters   = 0;
  while  (gamma_parameters[n_gamma_parameters] != 0){
    n_gamma_parameters  += 4;
  }
  n_components        = n_gamma_parameters / 4;
  bool*  bOneHit        = (bool*)calloc(n_components, sizeof(bool));

  for(i = 0; i < n_components; i++){
    if(bOneHit[i]){

      ////////////////////////////////////////////////////////////////////////////////////////////
      /*
        gsl_set_error_handler_off();

      double ext_integral_Gy;
      double error;
      gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (10000);
      gsl_function F;
      F.function = &AT_RDD_Katz_ext_integrand_Gy;
//      float params[] = {*r_m,*a0_m,*alpha,*r_min_m,*r_max_m,*Katz_point_coeff_Gy};
      F.params = params;
//      int status = gsl_integration_qags (&F, int_lim_m, (*r_m)+(*a0_m), 1e-9, 1e-4, 10000, w1, &ext_integral_Gy, &error);
      if (status == GSL_EROUND || status == GSL_ESING){
    #ifdef _DEBUG
        indnt_init();
        fprintf(debf,"%s r=%g, integration from %g to %g , error no == %d\n",isp,*r_m,int_lim_m,(*r_m)+(*a0_m),status);
    #endif
        ext_integral_Gy = -1.0f;
      }
      gsl_integration_workspace_free (w1);
       */

      ////////////////////////////////////////////////////////////////////////////////////////////
    }
  }
}


/**
* SPISS routine:
* Compound Poisson approach for local dose distribution
* similar to SPIFF but using statistical sampling
* in order to evaluate Fn(t)
* See SPIFF pamphlet 27.2./29.2.2008
* This was routine "LGC_full_simulation_C" in LGC 2.1
*/
void AT_SPISS(	const long* n,
    const float* E_MeV_u,
    const long*  particle_no,
    const float* fluence_cm2,
    const long*  material_no,
    const long*  RDD_model,
    const float* RDD_parameters,
    const long*  ER_model,
    const float* ER_parameters,
    const long*  gamma_model,
    const float* gamma_parameters,
    const long*	 n_runs,
    const long*  N2,
    const float* fluence_factor,
    const int*   write_output,
    const long*  importance_sampling,
    float*  results)
{
  // Welcome screen
  printf("\n############################################################\n");
  printf("\n############################################################\n");
  printf("This is AmTrack - SPISS algorithm\n");
  printf("\n");

  FILE*    output_file;
  output_file    =  fopen("SPISS.log","w");
  if (output_file == NULL) return;                      // File error

  // The histogram initialization and handling has been adapted to SPIFF
  // although some features are not used here
  long  				n_bins_f1;
  float*  			f1_parameters      =  (float*)calloc(9 * (*n), sizeof(float));

  AT_SC_get_f1_array_size(  	n,
      E_MeV_u,
      particle_no,
      material_no,
      RDD_model,
      RDD_parameters,
      ER_model,
      ER_parameters,
      N2,
      &n_bins_f1,
      f1_parameters);

  float*  f_parameters      				=  (float*)calloc(7, sizeof(float));

  float*  norm_fluence      				=  (float*)calloc(*n, sizeof(float));
  float*  accu_fluence      				=  (float*)calloc(*n, sizeof(float));
  float*  dose_contribution_Gy  			=  (float*)calloc(*n, sizeof(float));

  float*  f1_d_Gy          				=  (float*)calloc(n_bins_f1, sizeof(float));
  float*  f1_dd_Gy        				=  (float*)calloc(n_bins_f1, sizeof(float));
  float*  f1            					=  (float*)calloc(n_bins_f1, sizeof(float));

  AT_SC_get_f1(  	  n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      RDD_model,
      RDD_parameters,
      ER_model,
      ER_parameters,
      N2,
      &n_bins_f1,
      f1_parameters,
      // from here: return values
      norm_fluence,
      dose_contribution_Gy,
      f_parameters,
      f1_d_Gy,
      f1_dd_Gy,
      f1);

  // Get accumulated normalized fluence for later sampling of particle type
  accu_fluence[0]							=	norm_fluence[0];
  long i;
  if(*n > 1){
    for (i = 1; i < *n; i++){
      accu_fluence[i]     	 				+=  accu_fluence[i-1] + norm_fluence[i];
    }
  }

  long	n_bins_f;
  float   u_start;
  long    n_convolutions;
  float	u				= f_parameters[0];

  AT_SC_get_f_array_size(   &u,
      fluence_factor,
      N2,
      &n_bins_f1,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      // from here: return values
      &n_bins_f,
      &u_start,
      &n_convolutions);

  float*  f_d_Gy         	=  (float*)calloc(n_bins_f, sizeof(float));
  float*  f_dd_Gy        	=  (float*)calloc(n_bins_f, sizeof(float));
  float*  f            	=  (float*)calloc(n_bins_f, sizeof(float));
//  float*  fdd            	=  (float*)calloc(n_bins_f, sizeof(float)); // NOT USED
//  float*  dfdd          	=  (float*)calloc(n_bins_f, sizeof(float)); // NOT USED
  float   f0            	=  0.0f;

  AT_SC_get_f_start(  &u_start,
      &n_bins_f1,
      N2,
      f1_d_Gy,
      f1_dd_Gy,
      f1,
      &n_bins_f,
      // from here: return values
      f_d_Gy,
      f_dd_Gy,
      f);

  // We are only interested in f_d_Gy and f_dd_Gy, so clear f
  for (i = 0; i < n_bins_f; i++){
    f[i] = 0;
  }

  if(*importance_sampling){
    printf("\n");
    printf("Importance sampling chosen. Biasing function G(r)=r^%ld\n", *importance_sampling);}
  else{
    printf("\n");
    printf("No importance sampling chosen.\n");}

  // init RNG
  gsl_rng * rng1   = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rng1, 12345678);

  long	act_number_particles;
  float	d_Gy;
  float	d_j_Gy;
  float	weight;
  float	r_m;
  long	n_tmp = 1;
  float	F;
  long	bin_no;
  float	max_bin_Gy	= log10(f_d_Gy[n_bins_f-1]);
  float	min_bin_Gy	= log10(f_d_Gy[0]);
  float	dd_bin_Gy	= (max_bin_Gy - min_bin_Gy) / n_bins_f;

  // Do n_runs runs
  for (i = 0; i < *n_runs; i++){
    // Get actual number particles for this run from Poisson generator
    act_number_particles	  =   (long)gsl_ran_poisson(rng1, u);
    // Reset local dose for run i
    d_Gy		= 0.0f;
    // Reset weight for importance sampling
    weight 		= 1.0f;
    // Add n individual doses according to their distribution
    long j;
    for (j = 0; j < act_number_particles; j++){
      // (1) draw random number 0..1 and sample particle type
      F 		= (float)gsl_rng_uniform (rng1);
      long k;
      for (k = 0; k < *n; k++){
        if (accu_fluence[k] >= F){
          break;
        }
      }

      // (2) draw again random number 0..1 for radius sampling
      F = (float)gsl_rng_uniform (rng1);

      // (3) Apply importance sampling / weighting
      if (*importance_sampling){
        weight	*=	(*importance_sampling) * pow(F, (*importance_sampling)-1.0f);
        F		 =	pow(F, (*importance_sampling));
      }

      // (4) get dose d_Gy[j](r_max * F)
      r_m				= f1_parameters[k*9 + 2] * sqrt(F);						// r_max for particle type k * 0..1
      AT_D_RDD_Gy(  	&n_tmp,
          &r_m,
          &E_MeV_u[k],
          &particle_no[k],
          material_no,
          RDD_model,
          RDD_parameters,
          ER_model,
          ER_parameters,
          &d_j_Gy);

      // (5) Add dose
      d_Gy		+=	d_j_Gy;
    }

    // Fill dose into histogram
    if (d_Gy == 0.0f){
      f0 += weight;
    }
    else{
      bin_no = floor((log10(d_Gy) - min_bin_Gy + 3*dd_bin_Gy/2) / dd_bin_Gy);
      if (bin_no > n_bins_f) bin_no = n_bins_f;
      f[bin_no - 1]		+= weight / f_dd_Gy[bin_no - 1];
    }
    if(i%100 == 0){
      printf("Run %ld done.\n", i);
    }
  }

  // Normalize f
  float norm 		= 0.0f;
  float d_check	= 0.0f;
  for (i = 0; i < n_bins_f; i++){
    norm	+= f_dd_Gy[i] * f[i];
  }
  norm += f0;
  for (i = 0; i < n_bins_f; i++){
    f[i]	/= norm;
    d_check += f_d_Gy[i]*f_dd_Gy[i]*f[i];
  }

  fprintf(output_file, "SPISS\n");
  fprintf(output_file, "number of runs: %ld\n",    *n_runs);
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
  /*
	free(f1_parameters);
	free(f_parameters);
	free(norm_fluence);
	free(dose_contribution_Gy);
	free(accu_fluence);
	free(f1_d_Gy);
	free(f1_dd_Gy);
	free(f1);
	free(f_d_Gy);
//	free(f_dd_Gy);
	free(f);
	free(fdd);
	free(dfdd);
   */
  fclose(output_file);

  printf("\n");
  printf("AmTrack SPISS run finished.\n");
  printf("############################################################\n");
  printf("############################################################\n");
  return;
}
