/**
 * @brief AT_CPPSC wrapper
 */

/*
 *    AT_CPPSC.c
 *    ===================
 *
 *    Created on: 2010-10-11
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AT_Algorithms_CPP.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[10000];

	double E_MeV_u_single;

	long material_no = Water_Liquid;
	long particle_no_single = PARTICLE_PROTON_NUMBER;
	long rdd_model = RDD_Geiss;
	double rdd_parameters[RDD_MAX_NUMBER_OF_PARAMETERS];
	long er_model = ER_Geiss;
	long gamma_model = GR_GeneralTarget;
	double D0_Gy = 10.;
	double m = 1.;
	double c = 1.;
	double gamma_parameters[GR_MAX_NUMBER_OF_PARAMETERS];
	long source_no = PSTAR;

	double dose_Gy_single = 10.0;

	FILE *file;
	fflush(stdin);
	file = fopen(path, "a+");

	if (file == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), file) != 0) {
		if (strstr(Text, "E_MeV_u:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			E_MeV_u_single = atof(token);
		}
		if (strstr(Text, "dose_Gy:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			dose_Gy_single = atof(token);
		}
		if (strstr(Text, "material_no:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			material_no = atol(token);
		}
		if (strstr(Text, "particle_no:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			particle_no_single = atol(token);
		}
		if (strstr(Text, "rdd_model:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			rdd_model = atol(token);
		}
		if (strstr(Text, "er_model:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			er_model = atol(token);
		}
		if (strstr(Text, "D0_Gy:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			D0_Gy = atof(token);
		}
		if (strstr(Text, "c:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			c = atof(token);
		}
		if (strstr(Text, "m:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			m = atof(token);
		}
		if (strstr(Text, "source_no:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			source_no = atol(token);
		}
	}

	int i;

	long number_of_field_components = 1;
	double fluence_factor = 1.0;
	long N2 = 10;
	const bool    write_output = false;
	const bool    shrink_tails = true;
	const double  shrink_tails_under = 1e-30;
	const bool    adjust_N2 = true;
	const bool    lethal_events_mode = false;

	double       relative_efficiency;
    double       d_check;
	double       S_HCP;
	double       S_gamma;
	double       mean_number_of_tracks_contrib;
	double       start_number_of_tracks_contrib;
	long         n_convolutions;
	double		 lower_Jensen_bound;
	double       upper_Jensen_bound;

	int rdd_index = AT_RDD_index_from_RDD_number(rdd_model);
	for( i = 0; i < RDD_MAX_NUMBER_OF_PARAMETERS; i++){
		rdd_parameters[i] = AT_RDD_Data.parameter_default[rdd_index][i];
	}

	int gr_index = AT_Gamma_index_from_material_number(gamma_model);
	for( i = 0; i < GR_MAX_NUMBER_OF_PARAMETERS; i++){
		gamma_parameters[i] = AT_GR_Data.parameter_default[gr_index][i];
	}
	gamma_parameters[1] = D0_Gy;
	gamma_parameters[2] = c;
	gamma_parameters[3] = m;

	double fluence_cm2_or_dose_Gy_single = -dose_Gy_single;

	AT_run_CPPSC_method(  number_of_field_components,
	    &E_MeV_u_single,
	    &particle_no_single,
	    &fluence_cm2_or_dose_Gy_single,
	    material_no,
	    source_no,
	    rdd_model,
	    rdd_parameters,
	    er_model,
	    gamma_model,
	    gamma_parameters,
	    N2,
	    fluence_factor,
	    write_output,
	    shrink_tails,
	    shrink_tails_under,
	    adjust_N2,
	    lethal_events_mode,
	    &relative_efficiency,
	    &d_check,
	    &S_HCP,
	    &S_gamma,
	    &mean_number_of_tracks_contrib,
	    &start_number_of_tracks_contrib,
	    &n_convolutions,
	    &lower_Jensen_bound,
	    &upper_Jensen_bound);

	fprintf(file, "RE: %g\n", relative_efficiency);

	/* Get array size for single impact dose
	 * distribution for later memory allocation */
	long     n_bins_f1 = AT_n_bins_for_single_impact_local_dose_distrib(  number_of_field_components,
			&E_MeV_u_single,
			&particle_no_single,
			material_no,
			rdd_model,
			rdd_parameters,
			er_model,
			N2,
			source_no);

	/* Get f1 parameters - containing the most
	 * relevant information on the tracks of
	 * the mixed field, such as min/max local
	 * dose, track radius etc. */
	double*  f1_parameters      =  (double*)calloc(AT_SC_F1_PARAMETERS_SINGLE_LENGTH * number_of_field_components, sizeof(double));
	AT_RDD_f1_parameters_mixed_field( number_of_field_components,
			&E_MeV_u_single,
			&particle_no_single,
			material_no,
			rdd_model,
			rdd_parameters,
			er_model,
			source_no,
			f1_parameters
	);

	/* Get local dose distribution for the impact of a single particle */
	double*  f1_d_Gy       =  (double*)calloc(n_bins_f1, sizeof(double));
	double*  f1_dd_Gy      =  (double*)calloc(n_bins_f1, sizeof(double));
	double*  f1            =  (double*)calloc(n_bins_f1, sizeof(double));
	AT_single_impact_local_dose_distrib(  number_of_field_components,
			&E_MeV_u_single,
			&particle_no_single,
			&fluence_cm2_or_dose_Gy_single,
			material_no,
			rdd_model,
			rdd_parameters,
			er_model,
			N2,
			n_bins_f1,
			f1_parameters,
			source_no,
			f1_d_Gy,
			f1_dd_Gy,
			f1);

	fprintf(file , "f1_d_Gy:");
	for (i = 0; i < n_bins_f1; i++) {
		fprintf(file, " %g", f1_d_Gy[i]);
	}
	fprintf(file, "\n");
	fprintf(file , "f1:");
	for (i = 0; i < n_bins_f1; i++) {
		fprintf(file, " %g", f1_d_Gy[i]);
	}
	fprintf(file, "\n");


	/* Compute the mean number of tracks that
	 * deposit dose in a representative point
	 * of the detector/cell */

	double fluence_cm2 = AT_fluence_cm2_from_dose_Gy_single(
			E_MeV_u_single,
			particle_no_single,
			dose_Gy_single,
			material_no,
			source_no);
	mean_number_of_tracks_contrib  =       AT_mean_number_of_tracks_contrib(     number_of_field_components,
			&E_MeV_u_single,
			&particle_no_single,
			&fluence_cm2,
			material_no,
			er_model,
			source_no);

	/* Get array size for low fluence local dose
	 * distribution for later memory allocation */
	long      n_bins_f;
	AT_n_bins_for_low_fluence_local_dose_distribution(  mean_number_of_tracks_contrib,
			fluence_factor,
			N2,
			n_bins_f1,
			f1_d_Gy,
			f1_dd_Gy,
			f1,
			&n_bins_f,
			&start_number_of_tracks_contrib,
			&n_convolutions);

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
	AT_SuccessiveConvolutions(  mean_number_of_tracks_contrib,
			n_bins_f,
			&N2,
			&n_bins_f1,
			f_d_Gy,
			f_dd_Gy,
			f,
			&f0,
			fdd,
			dfdd,
			&d_check,
			write_output,
			shrink_tails,
			shrink_tails_under,
			adjust_N2);

	for (i = 0; i < n_bins_f; i++) {
		if( fdd[i] <= 0.0) f[i] = 0.0;
		if( f[i] <= 0.0) fdd[i] = 0.0;
	}


	fprintf(file , "fdd:");
	for (i = 0; i < n_bins_f; i++) {
		if( f_d_Gy[i] > 0. )
			fprintf(file, " %g", f_d_Gy[i]);
	}
	fprintf(file, "\n");
	fprintf(file , "f:");
	for (i = 0; i < n_bins_f1; i++) {
		if( f[i] > 0. )
			fprintf(file, " %g", f[i]);
	}
	fprintf(file, "\n");
	fclose(file);

	/* Free allocated memory and return*/
	free(f_d_Gy);
	free(f_dd_Gy);
	free(f);
	free(fdd);
	free(dfdd);

	free(f1_d_Gy);
	free(f1_dd_Gy);
	free(f1);

	free(f1_parameters);

	return EXIT_SUCCESS;
}
