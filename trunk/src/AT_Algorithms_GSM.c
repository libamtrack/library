/**
 * @file
 * @brief TODO
 */

/*
 *    AT_Algorithms_GSM.c
 *    =========
 *
 *    Created on: 30.09.2010
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

#include "AT_Algorithms_GSM.h"

void AT_GSM_shoot_particles_on_grid(  const long  number_of_field_components,
                const double         fluence_cm2[],
                const double         sample_grid_size_m,
                unsigned long* 		 random_number_generator_seed,
                long                 number_of_particles_in_field_component[],
                double*              x_position[],
                double*              y_position[]){

        /* Create and initialize random number generator */
        gsl_rng * rng  =  gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(rng, *random_number_generator_seed);

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

        // get random integer as next seed
        *random_number_generator_seed = gsl_rng_get(rng);

        gsl_rng_free(rng);
}


void AT_GSM_calculate_dose_pattern( const long  number_of_field_components,
    const double   E_MeV_u[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     number_of_particles_in_field_component[],
    const double*  x_position[],
    const double*  y_position[],
    const long     nX,
    const double   pixel_size_m,
    double**       grid_D_Gy){

  long i,j,k,l;

  /* calculate maximum delta-electron range for all components */
  double* r_max_m   = (double*)calloc(number_of_field_components, sizeof(double));
  AT_max_electron_ranges_m( number_of_field_components, E_MeV_u, material_no, er_model, r_max_m );

  /* find maximum of maximal delta-electron ranges */
  double max_r_max_m = 0.0;
  for (i = 0; i < number_of_field_components; i++){
    max_r_max_m    =   GSL_MAX(max_r_max_m, r_max_m[i]);
  }

  /* allocate memory for vector of vectors distances between track cores and grid points */
  double**  distances   = (double**)calloc(number_of_field_components, sizeof(double*));

  /* table of number (one number per component of multi field) of visible pixels from all particles */
  long* number_of_visible_pixels   = (long*)calloc(number_of_field_components, sizeof(long));

  /* vector of indices for distances vectors */
  long* distances_index = (long*)calloc(number_of_field_components, sizeof(long));

  for (i = 0; i < number_of_field_components; i++){
    number_of_visible_pixels[i] = 0;
    distances_index[i] = 0;
  }

  /* calculate how many grid pixels are visible for all particles from given component */
  for (j = 0; j < nX; j++){          // y
    double cur_y_pos   =  max_r_max_m + ((double)j + 0.5)*pixel_size_m;
    for (i = 0; i < nX; i++){        // x
      double cur_x_pos =  max_r_max_m + ((double)i + 0.5)*pixel_size_m;
      for (k = 0; k < number_of_field_components; k++){  // field components
        for (l = 0; l < number_of_particles_in_field_component[k]; l++){  // particles
          double r_m  =  sqrt( gsl_pow_2(x_position[k][l] - cur_x_pos) + gsl_pow_2(y_position[k][l] - cur_y_pos));
          if(r_m <= r_max_m[k]){    // does particle contribute?
            number_of_visible_pixels[k]++;
          }
        } // particle contribution
      } // field components
    } // x loop
  } // y loop

  for (k = 0; k < number_of_field_components; k++){
    distances[k] = (double*)calloc( number_of_visible_pixels[k], sizeof(double) );
  }

  /* save distances between visible grid pixels and all particles from given component in linear vectors */
  for (j = 0; j < nX; j++){          // y
    double cur_y_pos   =  max_r_max_m + ((double)j + 0.5)*pixel_size_m;
    for (i = 0; i < nX; i++){        // x
      double cur_x_pos =  max_r_max_m + ((double)i + 0.5)*pixel_size_m;
      for (k = 0; k < number_of_field_components; k++){  // field components
        for (l = 0; l < number_of_particles_in_field_component[k]; l++){  // particles
          double r_m  =  sqrt( gsl_pow_2(x_position[k][l] - cur_x_pos) + gsl_pow_2(y_position[k][l] - cur_y_pos));
          if(r_m <= r_max_m[k]){    // does particle contribute?
            distances[k][distances_index[k]] = r_m;
            distances_index[k]++;
          }
        } // particle contribution
      } // field components
    } // x loop
  } // y loop

  /* for every component calculate doses delivered to points in given distance from particle track */
  double**  doses = (double**)calloc(number_of_field_components, sizeof(double*));
  for (i = 0; i < number_of_field_components; i++){
    doses[i] = (double*)calloc( number_of_visible_pixels[i], sizeof(double) );

    AT_D_RDD_Gy( number_of_visible_pixels[i],
        distances[i],
        E_MeV_u[i],
        particle_no[i],
        material_no,
        rdd_model,
        rdd_parameter,
        er_model,
        doses[i]);
  }

  for (i = 0; i < number_of_field_components; i++){
    free(distances[i]);
    distances_index[i] = 0;
  }

  free(distances);
  free(number_of_visible_pixels);

  /* sum up doses for every grid cell from every component */
  for (j = 0; j < nX; j++){          // y
    double cur_y_pos   =  max_r_max_m + ((double)j + 0.5)*pixel_size_m;
    for (i = 0; i < nX; i++){        // x
      double cur_x_pos =  max_r_max_m + ((double)i + 0.5)*pixel_size_m;
      grid_D_Gy[i][j]  =  0.0;

      for (k = 0; k < number_of_field_components; k++){  // field components
        for (l = 0; l < number_of_particles_in_field_component[k]; l++){  // particles
          double r_m  =  sqrt( gsl_pow_2(x_position[k][l] - cur_x_pos) + gsl_pow_2(y_position[k][l] - cur_y_pos));
          if(r_m <= r_max_m[k]){    // does particle contribute?
            grid_D_Gy[i][j] += doses[k][distances_index[k]];
            distances_index[k]++;
          }
        } // particle contribution
      } // field components
    } // x loop
  } // y loop

  for (i = 0; i < number_of_field_components; i++){
    free(doses[i]);
  }
  free(doses);
  free(distances_index);
  free(r_max_m);
}


void AT_GSM_calculate_histogram_from_grid( const long     nX,
    const double** grid,
    const long     number_of_bins,
    const double   bin_centers_Gy[],
    double *       zero_fraction,
    double         frequency[]){

  long i, j, bin_no;

  // Fill dose into histogram
  /* On histograms and binning in libamtrack, see documentation */
  for(i = 0; i < nX; i++){
    for(j = 0; j < nX; j++){
      if (grid[i][j] == 0.0){
        (*zero_fraction) += 1.0;
      }
      else{
        bin_no = AT_histo_bin_no(	number_of_bins,
									bin_centers_Gy,
									grid[i][j]);
        frequency[bin_no]  += 1.0;
      }
    }
  }
  for(i = 0; i < number_of_bins; i++){
    frequency[i]   /= gsl_pow_2(nX);
  }
  (*zero_fraction) /= gsl_pow_2(nX);
}

void AT_GSM_calculate_dose_histogram( const long  number_of_field_components,
    const double   E_MeV_u[],
    const double   fluence_cm2[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     nX,
    const double   pixel_size_m,
    const long     number_of_bins,
    const double   dose_bin_centers_Gy[],
    unsigned long* random_number_generator_seed,
    double *       zero_dose_fraction,
    double         dose_frequency_Gy[]){

  long    i;

  /* find maximum of maximal delta-electron ranges */
  double max_r_max_m = 0.0;
  for (i = 0; i < number_of_field_components; i++){
    max_r_max_m    =   GSL_MAX(max_r_max_m, AT_max_electron_range_m(E_MeV_u[i], material_no, er_model));
  }

  /* largest r.max --> calculate size of sample area */
  double sample_grid_size_m    = pixel_size_m * nX + 2.01 * max_r_max_m;

  long*  number_of_particles_in_field_component   =  (long*)calloc(number_of_field_components, sizeof(double));
  double** x_position = (double**)calloc(number_of_field_components, sizeof(double*));
  double** y_position = (double**)calloc(number_of_field_components, sizeof(double*));

  /* linearly allocated 2-D arrays, see http://c-faq.com/aryptr/dynmuldimary.html */
  double** grid_D_Gy = (double**)calloc(nX, sizeof(double*));
  grid_D_Gy[0] = (double*)calloc(nX * nX, sizeof(double));
  for(i = 1; i < nX; i++)
    grid_D_Gy[i] = grid_D_Gy[0] + i * nX;

  /* find random positions of particles on grid
   * allocate xy_position tables */
  AT_GSM_shoot_particles_on_grid( number_of_field_components,
                  fluence_cm2,
                  sample_grid_size_m,
                  random_number_generator_seed,
                  number_of_particles_in_field_component,
                  x_position,
                  y_position);

  /* calculate dose deposition pattern in grid cells */
  AT_GSM_calculate_dose_pattern( number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameter,
      er_model,
      number_of_particles_in_field_component,
      (const double**)x_position,
      (const double**)y_position,
      nX,
      pixel_size_m,
      grid_D_Gy);

  /* calculate dose frequency from dose pattern */
  AT_GSM_calculate_histogram_from_grid( nX,
      (const double**)grid_D_Gy,
      number_of_bins,
      dose_bin_centers_Gy,
      zero_dose_fraction,
      dose_frequency_Gy);

  /* free memory */
  for (i = 0; i < number_of_field_components; i++){
    free( x_position[i] );
    free( y_position[i] );
  }

  free( grid_D_Gy[0] );
  free( grid_D_Gy );

  free( number_of_particles_in_field_component );

  free( x_position );
  free( y_position );

}


void AT_GSM_calculate_local_response_grid( const long      nX,
    const long      gamma_model,
    const double    gamma_parameters[],
    const double**  grid_D_Gy,
    const bool      lethal_events_mode,
    double**        grid_response){

/* lethal_events_mode (dGSM) works with arbitrary gamma models */
/* It is however more efficient to use GR_LinQuad_Log if LQ is used */
//  if( lethal_events_mode && (gamma_model != GR_LinQuad) ){
//    printf("Sorry, no Grid Summation model (lethal events mode) with other than Linear Quadratic gamma response\n");
//    printf("Please choose gamma_model = %d. Exiting now...\n", GR_LinQuad);
//    return;
//  }

  // get gamma response for local dose
  long i;
  for (i = 0; i < nX; i++){
    AT_gamma_response( nX,
        grid_D_Gy[i],
        gamma_model,
        gamma_parameters,
        lethal_events_mode,
        grid_response[i]);
  } // i

}


void AT_run_GSM_method(  const long  n,
    const double   E_MeV_u[],
    const long     particle_no[],
    const double   fluence_cm2_or_dose_Gy[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameters[],
    const long     er_model,
    const long     gamma_model,
    const double   gamma_parameters[],
    const long     N_runs,
    const bool     write_output,
    const long     nX,
    const double   voxel_size_m,
    const bool     lethal_events_mode,
    double         results[10])
{

  /** in case of debugging
  printf("n = %ld\n", n);
  printf("E_MeV_u[] = %g, %g, %g\n", E_MeV_u[0], E_MeV_u[1], E_MeV_u[2]);
  printf("particle_no[] = %ld, %ld, %ld\n", particle_no[0], particle_no[1], particle_no[2]);
  printf("fluence_cm2[] = %g, %g, %g\n", fluence_cm2[0], fluence_cm2[1], fluence_cm2[2]);
  printf("AT_material_no = %ld\n", AT_material_no);
  printf("rdd_model = %ld\n", rdd_model);
  printf("rdd_parameters[] = %g\n", rdd_parameters[0]);
  printf("er_model = %ld\n", er_model);
  printf("gamma_model = %ld\n", gamma_model);
  printf("gamma_parameters[] = %g, %g, %g, %g, %g\n", gamma_parameters[0], gamma_parameters[1], gamma_parameters[2], gamma_parameters[3], gamma_parameters[4]);
  printf("N_runs = %ld\n", N_runs);
  printf("write_output = %d\n", write_output);
  printf("nX = %ld\n", nX);
  printf("voxel_size_m = %g\n", voxel_size_m);
  printf("lethal_events_mode = %d\n", lethal_events_mode);
  */

  long    i, j, k;

  /* Zero results */
  for (i = 0; i < 10; i++){
    results[i] = 0.0;
  }

  double*  fluence_cm2    =  (double*)calloc(n, sizeof(double));

  if(fluence_cm2_or_dose_Gy[0] < 0){
    double*  dose_Gy        =  (double*)calloc(n, sizeof(double));
    for (i = 0; i < n; i++){
      dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    // dose to fluence
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

  /* find maximum of maximal delta-electron ranges */
  double max_r_max_m = 0.0;
  for (i = 0; i < n; i++){
    max_r_max_m    =   GSL_MAX(max_r_max_m, AT_max_electron_range_m(E_MeV_u[i], material_no, er_model));
  }

  /* largest r.max --> calculate size of sample area */
  double sample_grid_size_m    = voxel_size_m * nX + 2.01 * max_r_max_m;

  long*  number_of_particles_in_field_component   =  (long*)calloc(n, sizeof(double));
  double** x_position = (double**)calloc(n, sizeof(double*));
  double** y_position = (double**)calloc(n, sizeof(double*));

  /* linearly allocated 2-D arrays, see http://c-faq.com/aryptr/dynmuldimary.html */
  double** grid_D_Gy = (double**)calloc(nX, sizeof(double*));
  grid_D_Gy[0] = (double*)calloc(nX * nX, sizeof(double));
  for(i = 1; i < nX; i++)
    grid_D_Gy[i] = grid_D_Gy[0] + i * nX;

  double** grid_response = (double**)calloc(nX, sizeof(double*));
  grid_response[0] = (double*)calloc(nX * nX, sizeof(double));
  for(i = 1; i < nX; i++)
    grid_response[i] = grid_response[0] + i * nX;

  gsl_rng * rng  =  gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, 137);
  unsigned long random_number_generator_seed = gsl_rng_get(rng);

  /* start loop over N_runs */
  for (k = 0; k < N_runs; k++){

    /* find random positions of particles on grid
     * allocate xy_position tables */

	AT_GSM_shoot_particles_on_grid( n,
                    fluence_cm2,
                    sample_grid_size_m,
                    &random_number_generator_seed,
                    number_of_particles_in_field_component,
                    x_position,
                    y_position);

    /* find total number of particles */
    long n_particles = 0;
    for (i = 0; i < n; i++){
      n_particles  += number_of_particles_in_field_component[i];
    }

    /* calculate dose deposition pattern in grid cells */
    AT_GSM_calculate_dose_pattern( n,
        E_MeV_u,
        particle_no,
        material_no,
        rdd_model,
        rdd_parameters,
        er_model,
        number_of_particles_in_field_component,
        (const double**)x_position,
        (const double**)y_position,
        nX,
        voxel_size_m,
        grid_D_Gy);

    for (i = 0; i < n; i++){
      free( x_position[i] );
      free( y_position[i] );
    }

    /* For debugging: write grid to disk */
    if(k == 0){
    	  FILE*    output_file = NULL;
    	  output_file    =  fopen("GSMdoseGrid.csv","w");
    	  if (output_file == NULL) return;                      // File error
   	      fprintf(output_file, "x.m;y.m;d.Gy\n");
   	    	  for (j = 0; j < nX; j++){          // y
   	    	    double cur_y_pos   =  max_r_max_m + ((double)j + 0.5)*voxel_size_m;
   	    	    for (i = 0; i < nX; i++){        // x
   	    	      double cur_x_pos =  max_r_max_m + ((double)i + 0.5)*voxel_size_m;
   	    	      fprintf(output_file,
   	    	    		  "%e;%e;%e\n",
   	    	    		  cur_x_pos,
   	    	    		  cur_y_pos,
   	    	    		  grid_D_Gy[i][j]);
   	    	    }
   	    	  }
   	       fclose(output_file);
    }

    /* calculate response pattern in grid cells, knowing dose pattern */
    AT_GSM_calculate_local_response_grid( nX,
        gamma_model,
        gamma_parameters,
        (const double**)grid_D_Gy,
        lethal_events_mode,
        grid_response);

    /* find average dose and average response */
    double total_dose_on_grid_ions_Gy = 0.0;
    double average_grid_response_ions = 0.0;

    for (i = 0; i < nX; i++){
      for (j = 0; j < nX; j++){
        total_dose_on_grid_ions_Gy  += grid_D_Gy[i][j];
        average_grid_response_ions  += grid_response[i][j];
      }
    }

    total_dose_on_grid_ions_Gy  /= gsl_pow_2(nX);
    average_grid_response_ions  /= gsl_pow_2(nX);

    if( lethal_events_mode ){
    	average_grid_response_ions  = exp(-1.0 * average_grid_response_ions);

    	if( gamma_model == GR_ExpSaturation ||
    			gamma_model == GR_GeneralTarget ||
    			gamma_model == GR_Radioluminescence ||
    			gamma_model == GR_Geiss){
    		average_grid_response_ions  = 1.0 - average_grid_response_ions;
    	}
    }

    /* calculate gamma response for total grid dose */
    double gamma_response  = 0.0;
    AT_gamma_response(  1,
        &total_dose_on_grid_ions_Gy,
        gamma_model,
        gamma_parameters,
        false,
        &gamma_response);

    /* finally calculate efficiency */
    double efficiency  = 0.0;
    if(gamma_response > 0){
      efficiency = average_grid_response_ions / gamma_response;
    }

    /* add to results */
    results[0]      += efficiency;
    results[1]      += total_dose_on_grid_ions_Gy;
    results[2]      += average_grid_response_ions;
    results[3]      += gamma_response;
    results[4]      += n_particles;

    results[5]      += gsl_pow_2(efficiency);
    results[6]      += gsl_pow_2(total_dose_on_grid_ions_Gy);
    results[7]      += gsl_pow_2(average_grid_response_ions);
    results[8]      += gsl_pow_2(gamma_response);
    results[9]      += gsl_pow_2(n_particles);

    /* Write intermediate results to file */

//    FILE*    output_file = NULL;
//    output_file    =  fopen("GSMdoseGrid.csv","w");
//    if (output_file == NULL) return;                      // File error
//    fprintf(output_file, "x.m;y.m;d.Gy\n");

  }// end run loop

  /* average over number of runs */
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

  results[5]  -= gsl_pow_2(results[0]);
  results[6]  -= gsl_pow_2(results[1]);
  results[7]  -= gsl_pow_2(results[2]);
  results[8]  -= gsl_pow_2(results[3]);
  results[9]  -= gsl_pow_2(results[4]);

  results[5]  = GSL_MAX(0., results[5]);
  results[6]  = GSL_MAX(0., results[6]);
  results[7]  = GSL_MAX(0., results[7]);
  results[8]  = GSL_MAX(0., results[8]);
  results[9]  = GSL_MAX(0., results[9]);

  if( N_runs > 1 ){
          results[5]  = sqrt(results[5] / (N_runs - 1.));
          results[6]  = sqrt(results[6] / (N_runs - 1.));
          results[7]  = sqrt(results[7] / (N_runs - 1.));
          results[8]  = sqrt(results[8] / (N_runs - 1.));
          results[9]  = sqrt(results[9] / (N_runs - 1.));
  }

  /* free memory */
  gsl_rng_free( rng );

  free( grid_D_Gy[0] );
  free( grid_D_Gy );

  free( grid_response[0] );
  free( grid_response );

  free( fluence_cm2 );
  free( number_of_particles_in_field_component );

  free( x_position );
  free( y_position );
}
