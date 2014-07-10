/**
 * @brief AT_E_from_gamma wrapper
 */

/*
 *    AT_E_from_gamma.c
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

#include "AT_PhysicsRoutines.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[10000];

	double gamma_start = 0.;
	double gamma_end = 0.;
	long n_points = 0;
	long x_axis_type = 0;

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "gamma_start:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			gamma_start = atof(token);
		}
		if (strstr(Text, "gamma_end:")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			gamma_end = atof(token);
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

	if( n_points < 1 ){
		fprintf(stderr, "Number of points should be greater than 0, but is equal to %ld\n", n_points);
		fclose(f);
		return EXIT_FAILURE;
	}

	double * gamma = (double*)calloc(n_points, sizeof(double));
	double * E    = (double*)calloc(n_points, sizeof(double));

	long i;

	if( x_axis_type == 2){
		if( n_points > 1){
			for (i = 0; i < n_points; i++) {
				gamma[i] = gamma_start + (i/(double)(n_points-1)) * (gamma_end - gamma_start);
			}
		} else {
			gamma[0] = gamma_start;
		}
	} else if( x_axis_type == 1){
		if( n_points > 1 ){
			for (i = 0; i < n_points; i++) {
				double loggamma = log(gamma_start) + (i/(double)(n_points-1)) * (log(gamma_end) - log(gamma_start));
				gamma[i] = exp(loggamma);
			}
		} else {
			gamma[0] = gamma_start;
		}
	} else {
		fprintf(stderr, "X axis spacing type %ld not supported\n", x_axis_type);
		fclose(f);
		free(gamma);
		free(E);
		return EXIT_FAILURE;
	}

	int status = AT_E_from_gamma(n_points, gamma, E);

	if( status != EXIT_SUCCESS){
		fprintf(stderr, "Exit code from AT_E_from_gamma = %d\n", status);
		fclose(f);
		free(gamma);
		free(E);
		return EXIT_FAILURE;
	}

	fprintf(f , "gamma:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", gamma[i]);
	}
	fprintf(f, "\n");

	fprintf(f , "E:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", E[i]);
	}
	fprintf(f, "\n");

	fclose(f);
	free(gamma);
	free(E);

	return EXIT_SUCCESS;
}
