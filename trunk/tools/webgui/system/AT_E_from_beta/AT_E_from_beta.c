/**
 * @brief AT_E_from_beta wrapper
 */

/*
 *    AT_E_from_beta.c
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
#include "AT_Algorithms_IGK"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[600];

	long n = 0;
	double beta[500];

	double beta_start;
	double beta_end;
	long n_points;
	long x_axis_type;

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "beta_start")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			beta_start = atof(token);
		}
		if (strstr(Text, "beta_end")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			beta_end = atof(token);
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
	}

	double E[500];
	int i;
	if( x_axis_type == 2){
		for (i = 0; i < n_points; i++) {
			beta[i] = beta_start + (i/(double)(n_points-1)) * (beta_end - beta_start);
		}
	} else if( x_axis_type == 1){
		for (i = 0; i < n_points; i++) {
			double logbeta = log(beta_start) + (i/(double)(n_points-1)) * (log(beta_end) - log(beta_start));
			beta[i] = exp(logbeta);
		}
	} else {
		return EXIT_FAILURE;
	}

	AT_E_from_beta(n_points, beta, E);

	fprintf(f , "beta:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", beta[i]);
	}
	fprintf(f, "\n");

	fprintf(f , "E:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", E[i]);
	}
	fprintf(f, "\n");
	fclose(f);

	return EXIT_SUCCESS;
}
