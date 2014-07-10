/**
 * @brief AT_Gamma_Response wrapper
 */

/*
 *    AT_Gamma_Response.c
 *    ===================
 *
 *    Created on: 2011-04-27
 *    Creator: grzanka
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
#include "AT_GammaResponse.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[10000];

	double D_start = 0.1;
	double D_end = 10.;
	long n_points = 30;
	long gr_model = 1;
	long x_axis_type = 1;

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "D_start:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			D_start = atof(token);
		}
		if (strstr(Text, "D_end:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			D_end = atof(token);
		}
		if (strstr(Text, "n_points:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			n_points = atol(token);
		}
		if (strstr(Text, "x_axis_type:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			x_axis_type = atol(token);
		}
		if (strstr(Text, "gr_model:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			gr_model = atol(token);
		}
	}

	double * D    = (double*)calloc(n_points, sizeof(double));
	double * response = (double*)calloc(n_points, sizeof(double));

	long i;

	double gamma_parameters[GR_MAX_NUMBER_OF_PARAMETERS];
	int gr_index = AT_Gamma_index_from_material_number(gr_model);
	for( i = 0; i < GR_MAX_NUMBER_OF_PARAMETERS; i++){
		gamma_parameters[i] = AT_GR_Data.parameter_default[gr_index][i];
	}

	if( x_axis_type == 2){
		if( n_points > 1){
			for (i = 0; i < n_points; i++) {
				D[i] = D_start + (i/(double)(n_points-1)) * (D_end - D_start);
			}
		} else {
			D[0] = D_start;
		}
	} else if( x_axis_type == 1){
		if( D_start <= 0 ){
			fprintf(stderr, "D_start should be > 0\n");
			fclose(f);
			free(response);
			free(D);
			return EXIT_FAILURE;
		}
		if( n_points > 1 ){
			for (i = 0; i < n_points; i++) {
				double logD = log(D_start) + (i/(double)(n_points-1)) * (log(D_end) - log(D_start));
				D[i] = exp(logD);
			}
		} else {
			D[0] = D_start;
		}
	} else {
		fprintf(stderr, "X axis spacing type %ld not supported\n", x_axis_type);
		fclose(f);
		free(response);
		free(D);
		return EXIT_FAILURE;
	}

	AT_gamma_response(
			n_points,
			D,
			gr_model,
			gamma_parameters,
			false,
			response
	);

	fprintf(f , "D:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", D[i]);
	}
	fprintf(f, "\n");

	fprintf(f , "GR:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", response[i]);
	}
	fprintf(f, "\n");
	free(response);
	free(D);
	fclose(f);

	return EXIT_SUCCESS;
}
