/*
 * @file
 * @brief Dummy file to save dose pattern to output file.
 */

/*
 *    AT_test_GSM.c
 *    ===================
 *
 *    Created on: 2010-06-11
 *    Author: grzanka
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

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "AmTrack.h"
#include "AT_ElectronRange.h"

void AT_GSM_save_dose_pattern_to_file( const long  number_of_field_components,
    const double   E_MeV_u[],
    const double   fluence_cm2[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     nX,
    const double   pixel_size_m,
    const char *   filename){

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
			137,
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

	/* free memory */
	for (i = 0; i < number_of_field_components; i++){
		free( x_position[i] );
		free( y_position[i] );
	}

	FILE * outfile = fopen( filename, "w" );

	long j;
	for(i = 0; i < nX; i++)
		for(j = 0; j < nX; j++)
			fprintf(outfile,"%ld,%ld,%g\n",i,j,grid_D_Gy[i][j]);
	fclose(outfile);

	free( grid_D_Gy[0] );
	free( grid_D_Gy );

	free( number_of_particles_in_field_component );

	free( x_position );
	free( y_position );
}


int main( int argc, char * argv[] ){

	char filename1[200] = "grid1.dat";

	if( argc != 3){
		printf("Usage: %s nX pixel_size_m\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	printf("Saving dose pattern to file %s\n", filename1);

	long     number_of_field_components = 2;
	double   E_MeV_u[] = {1.,1.};
	double   fluence_cm2[] = {1e7,1e7};
	long     particle_no[] = {1001,6012};
	long     material_no = Water_Liquid;
	long     rdd_model   = RDD_Geiss;
	double   rdd_parameter[] = {5e-8};
	long     er_model    = ER_Waligorski;
	long     nX = atoi(argv[1]);
	double   pixel_size_m = atof(argv[2]);

	AT_GSM_save_dose_pattern_to_file( number_of_field_components,
	    E_MeV_u,
	    fluence_cm2,
	    particle_no,
	    material_no,
	    rdd_model,
	    rdd_parameter,
	    er_model,
	    nX,
	    pixel_size_m,
	    filename1);

	return EXIT_SUCCESS;
};
