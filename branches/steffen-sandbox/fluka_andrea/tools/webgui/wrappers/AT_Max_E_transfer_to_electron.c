/**
 * @brief AT_Max_E_transfer_to_electron wrapper
 */

/*
 *    AT_Max_E_transfer_to_electron.c
 *    ===================
 *
 *    Created on: 2011-06-29
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
#include "AT_PhysicsRoutines.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[10000];

	double E_start = 0.1;
	double E_end = 50.;
	long n_points = 30;
	long x_axis_type = 1;

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "E_start:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			E_start = atof(token);
		}
		if (strstr(Text, "E_end:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			E_end = atof(token);
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
	}

	double * E_MeV_u    = (double*)calloc(n_points, sizeof(double));
	double * electron_E_MeV = (double*)calloc(n_points, sizeof(double));

	long i;

	if( x_axis_type == 2){
		if( n_points > 1){
			for (i = 0; i < n_points; i++) {
				E_MeV_u[i] = E_start + (i/(double)(n_points-1)) * (E_end - E_start);
			}
		} else {
			E_MeV_u[0] = E_start;
		}
	} else if( x_axis_type == 1){
		if( n_points > 1 ){
			for (i = 0; i < n_points; i++) {
				double logE = log(E_start) + (i/(double)(n_points-1)) * (log(E_end) - log(E_start));
				E_MeV_u[i] = exp(logE);
			}
		} else {
			E_MeV_u[0] = E_start;
		}
	} else {
		fprintf(stderr, "X axis spacing type %ld not supported\n", x_axis_type);
		fclose(f);
		free(electron_E_MeV);
		free(E_MeV_u);
		return EXIT_FAILURE;
	}

	int status = AT_max_E_transfer_MeV( n_points,
			E_MeV_u,
			electron_E_MeV);


	if( status != EXIT_SUCCESS){
		fprintf(stderr, "Exit code from AT_Max_E_transfer_to_electron = %d\n", status);
		fclose(f);
		free(electron_E_MeV);
		free(E_MeV_u);
		return EXIT_FAILURE;
	}

	fprintf(f , "ionE:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", E_MeV_u[i]);
	}
	fprintf(f, "\n");

	fprintf(f , "electronE:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", electron_E_MeV[i]);
	}
	fprintf(f, "\n");
	free(electron_E_MeV);
	free(E_MeV_u);
	fclose(f);

	return EXIT_SUCCESS;
}
