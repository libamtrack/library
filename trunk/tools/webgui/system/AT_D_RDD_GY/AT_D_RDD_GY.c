/**
 * @brief AT_CSDA_range wrapper
 */

/*
 *    AT_CSDA_range.c
 *    ===================
 *
 *    Created on: 2010-10-11
 *    Creator: christophkolb
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
#include "AT_RDD.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[6000];

	double r_m[5000];
	double RDD_GY[5000];
	double rdd_parameters[RDD_MAX_NUMBER_OF_PARAMETERS];
	long material_no;
	long particle_no_single;
	long rdd_model;
	long er_model;
	double E_MeV_u;

	double r_start_m;
	double r_stop_m;
	long n_points;
	long x_axis_type;

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "r_min")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			r_start_m = atof(token);
		}
		if (strstr(Text, "r_max")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			r_stop_m = atof(token);
		}
		if (strstr(Text, "n_points")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			n_points = atol(token);
		}
		if (strstr(Text, "x_axis_type")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			x_axis_type = atol(token);
		}
		if (strstr(Text, "material_no")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			material_no = atol(token);
		}
		if (strstr(Text, "particle_no")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			particle_no_single = atol(token);
		}
		if (strstr(Text, "E_MeV_u")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			E_MeV_u = atof(token);
		}
		if (strstr(Text, "rdd_model")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			rdd_model = atol(token);
		}
		if (strstr(Text, "er_model")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			er_model = atol(token);
		}
	}

	int i;

	if( x_axis_type == 2){
		for (i = 0; i < n_points; i++) {
			r_m[i] = r_start_m + (i/(double)(n_points-1)) * (r_stop_m - r_start_m);
		}
	} else if( x_axis_type == 1){
		for (i = 0; i < n_points; i++) {
			double logE = log(r_start_m) + (i/(double)(n_points-1)) * (log(r_stop_m) - log(r_start_m));
			r_m[i] = exp(logE);
		}
	} else {
		return EXIT_FAILURE;
	}

	int rdd_index = AT_RDD_index_from_RDD_number(rdd_model);
	for( i = 0; i < RDD_MAX_NUMBER_OF_PARAMETERS; i++){
		rdd_parameters[i] = AT_RDD_Data.parameter_default[rdd_index][i];
	}

	AT_D_RDD_Gy( n_points,
			r_m,
			E_MeV_u,
			particle_no_single,
			material_no,
			rdd_model,
			rdd_parameters,
			er_model,
			RDD_GY);

	fprintf(f , "r:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", r_m[i]);
	}
	fprintf(f, "\n");

	fprintf(f , "RDD:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", RDD_GY[i]);
	}
	fprintf(f, "\n");
	fclose(f);

	return EXIT_SUCCESS;
}
