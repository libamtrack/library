/**
 * @brief AT_Katz wrapper
 */

/*
 *    AT_Katz.c
 *    ===================
 *
 *    Created on: 2010-11-06
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

#include "AT_KatzModel.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[10000];

	double E_MeV_u = 0.0;

	long material_no = Water_Liquid;
	long particle_no_single = PARTICLE_PROTON_NUMBER;
	long rdd_model = RDD_Geiss;
	double rdd_parameters[RDD_MAX_NUMBER_OF_PARAMETERS];
	long er_model = ER_Geiss;
	long source_no = PSTAR;
	double D0_Gy = 0.;
	double m = 0.;
	double sigma0_m2 = 0.;
	double a0_m = 0.;

	double D_Gy_start = 0.;
	double D_Gy_end = 0.;
	long n_points = 0;

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
			E_MeV_u = atof(token);
		}
		if (strstr(Text, "D_Gy_start:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			D_Gy_start = atof(token);
		}
		if (strstr(Text, "D_Gy_end:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			D_Gy_end = atof(token);
		}
		if (strstr(Text, "n_points:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			n_points = atol(token);
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
		if (strstr(Text, "source_no:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			source_no = atol(token);
		}
		if (strstr(Text, "a0_m:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			a0_m = atof(token);
		}
		if (strstr(Text, "D0_Gy:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			D0_Gy = atof(token);
		}
		if (strstr(Text, "m:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			m = atof(token);
		}
		if (strstr(Text, "sigma0_m2:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			sigma0_m2 = atof(token);
		}
	}

	long i;

	int rdd_index = AT_RDD_index_from_RDD_number(rdd_model);
	for( i = 0; i < RDD_MAX_NUMBER_OF_PARAMETERS; i++){
		rdd_parameters[i] = AT_RDD_Data.parameter_default[rdd_index][i];
	}
	rdd_parameters[1] = a0_m;

	if( n_points < 1 ){
		fprintf(stderr, "Number of points should be greater than 0, but is equal to %ld\n", n_points);
		fclose(file);
		return EXIT_FAILURE;
	}

	double * D_Gy      = (double*)calloc(n_points, sizeof(double));
	double * survival  = (double*)calloc(n_points, sizeof(double));

	if( n_points > 1){
		for (i = 0; i < n_points; i++) {
			D_Gy[i] = D_Gy_start + (i/(double)(n_points-1)) * (D_Gy_end - D_Gy_start);
		}
	} else {
		D_Gy[0] = D_Gy_start;
	}

	double * fluence_cm2 = (double*)calloc(n_points, sizeof(double));

	for (i = 0; i < n_points; i++) {
		fluence_cm2[i] = AT_fluence_cm2_from_dose_Gy_single(E_MeV_u, particle_no_single, D_Gy[i], material_no, PSTAR );
	}

	int status = AT_KatzModel_single_field_survival_optimized_for_fluence_vector(
			n_points,
			fluence_cm2,
			E_MeV_u,
			particle_no_single,
			material_no,
			rdd_model,
			rdd_parameters,
			er_model,
			D0_Gy,
			m,
			sigma0_m2,
			source_no,
			survival);

	free(fluence_cm2);

	if( status != EXIT_SUCCESS ){
		fprintf(stderr,"Exit code of AT_KatzModel_single_field_survival_optimized_for_fluence_vector is %d\n", status);
		fclose(file);
		free(D_Gy);
		free(survival);
		return EXIT_FAILURE;
	}

	fprintf(file , "D_Gy:");
	for (i = 0; i < n_points; i++) {
		fprintf(file, " %g", D_Gy[i]);
	}
	fprintf(file, "\n");
	fprintf(file , "survival:");
	for (i = 0; i < n_points; i++) {
		fprintf(file, " %g", survival[i]);
	}
	fprintf(file, "\n");
	fclose(file);
	free(D_Gy);
	free(survival);

	return EXIT_SUCCESS;
}
