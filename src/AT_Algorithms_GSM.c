/**
 * @brief GSM algorithm
 */

/*
 *    AT_Algorithms_GSM.c
 *    =========
 *
 *    Created on: 30.09.2010
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

#include "AT_Algorithms_GSM.h"
#include <time.h>

void AT_GSM_sample_particle_positions(  const long  number_of_field_components,
                const double         fluence_cm2[],
                const double         sample_grid_size_m,
                unsigned long* 		 random_number_generator_seed,
                long                 number_of_particles_in_field_component[],
                double*              x_position[],
                double*              y_position[]){

        printf("start AT_GSM_sample_particle_positions\n");

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
        printf("number of field components: %ld\n", number_of_field_components);
        for (i = 0; i < number_of_field_components; i++){
                mean_number_particles_in_field_component   =  sample_grid_area_cm2 * total_fluence_cm2 * norm_fluence[i];

                /* Actual number of particles from one component in the sample grid area, taken from Poissonian distribution */
                number_of_particles_in_field_component[i]  =  (long)gsl_ran_poisson(rng, mean_number_particles_in_field_component);

                /* Now we know how many particles from one component we have, so position tables can be allocated */
                printf("number of particles in field component %ld: %ld\n", i, number_of_particles_in_field_component[i]);
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
                printf("x_position[%ld][0] = %f\n", i, x_position[i][0]);
        }

        // get random integer as next seed
        *random_number_generator_seed = gsl_rng_get(rng);

        gsl_rng_free(rng);
        printf("stop AT_GSM_sample_particle_positions\n");
}


void AT_GSM_dose_grid_from_particles_positions( const long  number_of_field_components,
    const double   E_MeV_u[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     stopping_power_source_no,
    const long     number_of_particles_in_field_component[],
    const double*  x_position[],
    const double*  y_position[],
    const long     nX,
    const double   pixel_size_m,
    double**       grid_D_Gy){

  long i,j,k,l;
  time_t begin, end;
  printf("start AT_GSM_dose_grid_from_particles_positions\n");

  begin = time(NULL);
  /* calculate maximum delta-electron range for all components */
  double* r_max_m   = (double*)calloc(number_of_field_components, sizeof(double));
  AT_max_electron_ranges_m( number_of_field_components, E_MeV_u, material_no, er_model, r_max_m );

  end = time(NULL);
  printf("step1 time %f [s]\n", difftime(end, begin));

  /* find maximum of maximal delta-electron ranges */
  double max_r_max_m = 0.0;
  for (i = 0; i < number_of_field_components; i++){
    max_r_max_m    =   GSL_MAX(max_r_max_m, r_max_m[i]);
  }

  end = time(NULL);
  printf("step2 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

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

  end = time(NULL);
  printf("step3 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

  /* calculate how many grid pixels are visible for all particles from given component */
  for (j = 0; j < nX; j++){          // y
    double cur_y_pos   =  max_r_max_m + ((double)j + 0.5)*pixel_size_m;
    for (i = 0; i < nX; i++){        // x
      double cur_x_pos =  max_r_max_m + ((double)i + 0.5)*pixel_size_m;
      for (k = 0; k < number_of_field_components; k++){  // field components
        for (l = 0; l < number_of_particles_in_field_component[k]; l++){  // particles
          double r_squared_m  =  gsl_pow_2(x_position[k][l] - cur_x_pos) + gsl_pow_2(y_position[k][l] - cur_y_pos);
          if(r_squared_m <= gsl_pow_2(r_max_m[k])){    // does particle contribute?, we compare r^2 with r_max^2 to avoid CPU costly sqrt
            number_of_visible_pixels[k]++;
          }
        } // particle contribution
      } // field components
    } // x loop
  } // y loop

  end = time(NULL);
  printf("step4 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

  for (k = 0; k < number_of_field_components; k++){
    distances[k] = (double*)calloc( number_of_visible_pixels[k], sizeof(double) );
  }

  end = time(NULL);
  printf("step5 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

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
  
  end = time(NULL);
  printf("step6 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

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
        stopping_power_source_no,
        doses[i]);
  }

  end = time(NULL);
  printf("step7 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

  for (i = 0; i < number_of_field_components; i++){
    free(distances[i]);
    distances_index[i] = 0;
  }

  free(distances);
  free(number_of_visible_pixels);

  end = time(NULL);
  printf("step8 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

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
  
  end = time(NULL);
  printf("step9 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

  for (i = 0; i < number_of_field_components; i++){
    free(doses[i]);
  }

  end = time(NULL);
  printf("step10 time %f [s]\n", difftime(end, begin));
  begin = time(NULL);

  free(doses);
  free(distances_index);
  free(r_max_m);
  printf("stop AT_GSM_dose_grid_from_particles_positions\n");

}


void AT_GSM_local_dose_distrib_from_dose_grid( const long     nX,
    const double** grid,
    const long     number_of_bins,
    const double   bin_centers_Gy[],
    double *       zero_fraction,
    double         frequency[]){

  long i, j, bin_no;

  printf("start AT_GSM_local_dose_distrib_from_dose_grid\n");

  printf("step 1\n");

  // Fill dose into histogram
  /* On histograms and binning in libamtrack, see documentation */
  for(i = 0; i < nX; i++){
    for(j = 0; j < nX; j++){
      if (grid[i][j] == 0.0){
        (*zero_fraction) += 1.0;
      }
      else{
        bin_no = AT_histoOld_bin_no(	number_of_bins,
									bin_centers_Gy,
									grid[i][j]);
        // user could provide bins not spanning whole dose range
        // we skip the values outside the range
        if((bin_no >= 0) && (bin_no < number_of_bins)){
          frequency[bin_no] += 1.0;
        }
      }
    }
  }

  printf("step 2\n");
  for(i = 0; i < number_of_bins; i++){
    frequency[i]   /= gsl_pow_2(nX);
  }

  printf("step 3\n");
  (*zero_fraction) /= gsl_pow_2(nX);
  printf("step 4\n");
}


void AT_GSM_response_grid_from_dose_grid( const long      nX,
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



void AT_run_GSM_method(  const long  number_of_field_components,
    const double   E_MeV_u[],
    const long     particle_no[],
    const double   fluence_cm2_or_dose_Gy[],
    const long     material_no,
    const long     stopping_power_source_no,
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
    double*        relative_efficiency,
    double*        d_check,
    double*        S_HCP,
    double*        S_gamma,
    double*        average_n_particles,
    double*        sd_relative_efficiency,
    double*        sd_d_check,
    double*        sd_S_HCP,
    double*        sd_S_gamma,
    double*        sd_n_particles)
{
  long    i, j, k;

  /* Zero results */
  *relative_efficiency    = 0.0;
  *d_check                = 0.0;
  *S_HCP                  = 0.0;
  *S_gamma                = 0.0;
  *average_n_particles    = 0.0;
  *sd_relative_efficiency = 0.0;
  *sd_d_check             = 0.0;
  *sd_S_HCP               = 0.0;
  *sd_S_gamma             = 0.0;
  *sd_n_particles         = 0.0;

  /* Convert dose to fluence or
   * vice versa depending on
   * user's input */
  double*  fluence_cm2    =  (double*)calloc(number_of_field_components, sizeof(double));
  if(fluence_cm2_or_dose_Gy[0] < 0){
    double*  dose_Gy        =  (double*)calloc(number_of_field_components, sizeof(double));
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

  /* Find maximum of track radius in mixed field */
  double max_r_max_m = 0.0;
  for (i = 0; i < number_of_field_components; i++){
    max_r_max_m    =   GSL_MAX(max_r_max_m, AT_max_electron_range_m(E_MeV_u[i], material_no, er_model));
  }

  /* From the largest radius, calculate size
   * of sample area (grid) */
  double sample_grid_size_m                       = voxel_size_m * nX + 2.01 * max_r_max_m;
  long*  number_of_particles_in_field_component   =  (long*)calloc(number_of_field_components, sizeof(double));
  double** x_position = (double**)calloc(number_of_field_components, sizeof(double*));
  double** y_position = (double**)calloc(number_of_field_components, sizeof(double*));

  /* Prepare linearly allocated 2-D arrays for
   * dose on grid and response on grid, see
   * http://c-faq.com/aryptr/dynmuldimary.html */
  double** grid_D_Gy = (double**)calloc(nX, sizeof(double*));
  grid_D_Gy[0] = (double*)calloc(nX * nX, sizeof(double));
  for(i = 1; i < nX; i++)
    grid_D_Gy[i] = grid_D_Gy[0] + i * nX;

  double** grid_response = (double**)calloc(nX, sizeof(double*));
  grid_response[0] = (double*)calloc(nX * nX, sizeof(double));
  for(i = 1; i < nX; i++)
    grid_response[i] = grid_response[0] + i * nX;

  /* Initialize random number generator */
  gsl_rng * rng  =  gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, 137);
  unsigned long random_number_generator_seed = gsl_rng_get(rng);

  /* Main loop
   * Run simulation of throwing particles on the
   * grid as many time as needed to achieve
   * the requested statistical uncertainty
   * (not yet implemented) or for as many times
   * as requested by user */
  for (k = 0; k < N_runs; k++){

    /* Find random positions of particles on sample grid
     * and store in tables */
	AT_GSM_sample_particle_positions( number_of_field_components,
                    fluence_cm2,
                    sample_grid_size_m,
                    &random_number_generator_seed,
                    number_of_particles_in_field_component,
                    x_position,
                    y_position);

    /* Get the total number of particles in this run*/
    long n_particles = 0;
    for (i = 0; i < number_of_field_components; i++){
      n_particles  += number_of_particles_in_field_component[i];
    }

    /* Calculate dose deposition pattern in grid cells */
    AT_GSM_dose_grid_from_particles_positions( number_of_field_components,
        E_MeV_u,
        particle_no,
        material_no,
        rdd_model,
        rdd_parameters,
        er_model,
        stopping_power_source_no,
        number_of_particles_in_field_component,
        (const double**)x_position,
        (const double**)y_position,
        nX,
        voxel_size_m,
        grid_D_Gy);

    /* After doing so, the position grid
     * can be freed again */
    for (i = 0; i < number_of_field_components; i++){
      free( x_position[i] );
      free( y_position[i] );
    }

    /* For debugging/illustration: write
     * grid of first run to disk */
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

    /* Calculate response pattern on grid
	 * from known dose pattern. In case
	 * of lethal event mode use survival */
    AT_GSM_response_grid_from_dose_grid( nX,
        gamma_model,
        gamma_parameters,
        (const double**)grid_D_Gy,
        lethal_events_mode,
        grid_response);

    /* Compute particle dose and response
     * by averaging all entries on grid */
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

    /* In case of lethal event mode
     * compute survival */
    if( lethal_events_mode ){
    	average_grid_response_ions  = exp(-1.0 * average_grid_response_ions);

    	if( gamma_model == GR_ExpSaturation ||
    			gamma_model == GR_GeneralTarget ||
    			gamma_model == GR_Radioluminescence ||
    			gamma_model == GR_Geiss){
    		average_grid_response_ions  = 1.0 - average_grid_response_ions;
    	}
    }

    /* Calculate gamma response for total grid dose */
    double gamma_response  = 0.0;
    AT_gamma_response(  1,
        &total_dose_on_grid_ions_Gy,
        gamma_model,
        gamma_parameters,
        false,
        &gamma_response);

    /* Calculate relative efficiency */
    double efficiency  = 0.0;
    if(gamma_response > 0){
      efficiency = average_grid_response_ions / gamma_response;
    }

    /* Add run results to overall
     * results (running mean and variance) */
    *relative_efficiency     += efficiency;
    *d_check                 += total_dose_on_grid_ions_Gy;
    *S_HCP                   += average_grid_response_ions;
    *S_gamma                 += gamma_response;
    *average_n_particles     += n_particles;

    *sd_relative_efficiency  += gsl_pow_2(efficiency);
    *sd_d_check              += gsl_pow_2(total_dose_on_grid_ions_Gy);
    *sd_S_HCP                += gsl_pow_2(average_grid_response_ions);
    *sd_S_gamma              += gsl_pow_2(gamma_response);
    *sd_n_particles          += gsl_pow_2(n_particles);
  }/* End of main loop */

  /* Normalize results to number of runs */
  *relative_efficiency     /= N_runs;
  *d_check                 /= N_runs;
  *S_HCP                   /= N_runs;
  *S_gamma                 /= N_runs;
  *average_n_particles     /= N_runs;

  *sd_relative_efficiency  /= N_runs;
  *sd_d_check              /= N_runs;
  *sd_S_HCP                /= N_runs;
  *sd_S_gamma              /= N_runs;
  *sd_n_particles          /= N_runs;

  /* Compute errors (standard deviations)
   * from running estimators */
  *sd_relative_efficiency  -= gsl_pow_2(*relative_efficiency);
  *sd_d_check              -= gsl_pow_2(*d_check);
  *sd_S_HCP                -= gsl_pow_2(*S_HCP);
  *sd_S_gamma              -= gsl_pow_2(*S_gamma);
  *sd_n_particles          -= gsl_pow_2(*average_n_particles);

  /* Only report if more then
   * one run */
  *sd_relative_efficiency  = GSL_MAX(0., *sd_relative_efficiency);
  *sd_d_check              = GSL_MAX(0., *sd_d_check);
  *sd_S_HCP                = GSL_MAX(0., *sd_S_HCP);
  *sd_S_gamma              = GSL_MAX(0., *sd_S_gamma);
  *sd_n_particles          = GSL_MAX(0., *sd_n_particles);
  if( N_runs > 1 ){
	  *sd_relative_efficiency  = sqrt(*sd_relative_efficiency / (N_runs - 1.));
	  *sd_d_check              = sqrt(*sd_d_check / (N_runs - 1.));
	  *sd_S_HCP                = sqrt(*sd_S_HCP / (N_runs - 1.));
	  *sd_S_gamma              = sqrt(*sd_S_gamma / (N_runs - 1.));
	  *sd_n_particles          = sqrt(*sd_n_particles / (N_runs - 1.));
  }

  /* Free allocated memory
   * and return */
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

void AT_GSM_multiple_local_dose_distrib( const long  number_of_field_components,
    const double   	E_MeV_u[],
    const double   	fluence_cm2[],
    const long     	particle_no[],
    const long     	material_no,
    const long     	rdd_model,
    const double   	rdd_parameter[],
    const long     	er_model,
    const long      stopping_power_source_no,
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

  // print all input arguments
  printf("number_of_field_components = %ld\n", number_of_field_components);
  printf("E_MeV_u = %g\n", E_MeV_u[0]);
  printf("fluence_cm2 = %g\n", fluence_cm2[0]);
  printf("particle_no = %ld\n", particle_no[0]);
  printf("material_no = %ld\n", material_no);
  printf("rdd_model = %ld\n", rdd_model);
  printf("rdd_parameter = %g\n", rdd_parameter[0]);
  printf("er_model = %ld\n", er_model);
  printf("stopping_power_source_no = %ld\n", stopping_power_source_no);
  printf("nX = %ld\n", nX);
  printf("pixel_size_m = %g\n", pixel_size_m);
  printf("N_runs = %ld\n", N_runs);
  printf("N_repetitions = %ld\n", N_repetitions);
  printf("number_of_bins = %ld\n", number_of_bins);
  printf("dose_bin_centers_Gy = %g\n", dose_bin_centers_Gy[0]);
  printf("dose_bin_width_Gy = %g\n", dose_bin_width_Gy[0]);
  printf("mean_d_check_Gy = %g\n", *mean_d_check_Gy);
  printf("sd_d_check_Gy = %g\n", *sd_d_check_Gy);
  printf("mean_zero_dose_fraction = %g\n", *mean_zero_dose_fraction);
  printf("sd_zero_dose_fraction = %g\n", *sd_zero_dose_fraction);
  printf("mean_dose_frequency_Gy = %g\n", mean_dose_frequency_Gy[0]);
  printf("sd_dose_frequency_Gy = %g\n", sd_dose_frequency_Gy[0]);


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

			AT_GSM_local_dose_distrib( number_of_field_components,
					E_MeV_u,
					fluence_cm2,
					particle_no,
					material_no,
					rdd_model,
					rdd_parameter,
					er_model,
					stopping_power_source_no,
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

void AT_GSM_local_dose_distrib( const long  number_of_field_components,
    const double   E_MeV_u[],
    const double   fluence_cm2[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     stopping_power_source_no,
    const long     nX,
    const double   pixel_size_m,
    const long     number_of_bins,
    const double   dose_bin_centers_Gy[],
    unsigned long* random_number_generator_seed,
    double *       zero_dose_fraction,
    double         dose_frequency_Gy[]){

  long    i;

  // print all input parameters
  printf("INPUT: number_of_field_components = %ld\n", number_of_field_components);
  printf("INPUT: E_MeV_u = %g\n", E_MeV_u[0]);
  printf("INPUT: fluence_cm2 = %g\n", fluence_cm2[0]);
  printf("INPUT: particle_no = %ld\n", particle_no[0]);
  printf("INPUT: material_no = %ld\n", material_no);
  printf("INPUT: rdd_model = %ld\n", rdd_model);
  printf("INPUT: rdd_parameter = %g\n", rdd_parameter[0]);
  printf("INPUT: er_model = %ld\n", er_model);
  printf("INPUT: stopping_power_source_no = %ld\n", stopping_power_source_no);
  printf("INPUT: nX = %ld\n", nX);
  printf("INPUT: pixel_size_m = %g\n", pixel_size_m);
  printf("INPUT: number_of_bins = %ld\n", number_of_bins);
  printf("INPUT: dose_bin_centers_Gy[0] = %g\n", dose_bin_centers_Gy[0]);
  printf("INPUT: dose_bin_centers_Gy[number_of_bins-1] = %g\n", dose_bin_centers_Gy[number_of_bins-1]);
  printf("INPUT: random_number_generator_seed = %ld\n", *random_number_generator_seed);
  printf("INPUT: zero_dose_fraction = %g\n", *zero_dose_fraction);
  printf("INPUT: dose_frequency_Gy = %g\n", dose_frequency_Gy[0]);


  
  printf("AT_GSM_local_dose_distrib\n");
  printf("input: number_of_field_components: %ld \n", number_of_field_components);

  /* find maximum of maximal delta-electron ranges */
  double max_r_max_m = 0.0;
  for (i = 0; i < number_of_field_components; i++){
    max_r_max_m    =   GSL_MAX(max_r_max_m, AT_max_electron_range_m(E_MeV_u[i], material_no, er_model));
    printf("Field %ld / %ld,  Rmax = %g [m]\n", i, number_of_field_components, max_r_max_m);
  }
  printf("info: max_r_max_m: %g \n", max_r_max_m);

  /* largest r.max --> calculate size of sample area */
  double sample_grid_size_m    = pixel_size_m * nX + 2.01 * max_r_max_m;
  printf("info: sample_grid_size_m: %g \n", sample_grid_size_m);

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
  AT_GSM_sample_particle_positions( number_of_field_components,
                  fluence_cm2,
                  sample_grid_size_m,
                  random_number_generator_seed,
                  number_of_particles_in_field_component,
                  x_position,
                  y_position);

  printf("x_position: %p \n", x_position);
  printf("x_position[0]: %p \n", x_position[0]);

  /* calculate dose deposition pattern in grid cells */
  AT_GSM_dose_grid_from_particles_positions( number_of_field_components,
      E_MeV_u,
      particle_no,
      material_no,
      rdd_model,
      rdd_parameter,
      er_model,
      stopping_power_source_no,
      number_of_particles_in_field_component,
      (const double**)x_position,
      (const double**)y_position,
      nX,
      pixel_size_m,
      grid_D_Gy);

  printf("x_position: %p \n", x_position);
  printf("x_position[0]: %p \n", x_position[0]);


  /* calculate dose frequency from dose pattern */
  AT_GSM_local_dose_distrib_from_dose_grid( nX,
      (const double**)grid_D_Gy,
      number_of_bins,
      dose_bin_centers_Gy,
      zero_dose_fraction,
      dose_frequency_Gy);

  printf("x_position: %p \n", x_position);
  printf("x_position[0]: %p \n", x_position[0]);
  printf("Free 1\n");

  /* free memory */
  for (i = 0; i < number_of_field_components; i++){
    printf("i = %ld \n", i);
    printf("freeing x_position \n");
    printf("x_position[i] = %p \n", x_position[i]);
    printf("x_position[i][0] = %g \n", x_position[i][0]);
    free( x_position[i] );
    printf("freeing y_position \n");
    free( y_position[i] );
  }

  printf("x_position: %p \n", x_position);
  printf("x_position[0]: %p \n", x_position[0]);

  printf("Free 2\n");
  free( grid_D_Gy[0] );

  printf("Free 3\n");
  free( grid_D_Gy );

  printf("Free 4\n");
  free( number_of_particles_in_field_component );
  printf("Free 5\n");
  printf("x_position: %p \n", x_position);
  free( x_position );
  printf("Free 6\n");
  free( y_position );
  printf("Free 7\n");
}

