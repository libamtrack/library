/**
 * @brief AT_LET_from_E wrapper
 */

/*
 *    AT_LET_from_E.c
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
#include "AT_DataLET.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[600];

	double E_MeV_u[500];
	double LET_keV_um[500];
	long material_no = Water_Liquid;
	long particle_no[500];
	long particle_no_single = PARTICLE_PROTON_NUMBER;

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
		if (strstr(Text, "E_start")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			E_start = atof(token);
		}
		if (strstr(Text, "E_end")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			E_end = atof(token);
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
	}

	int i;
	for( i = 0 ; i < n_points ; i++){
		particle_no[i] = particle_no_single;
	}

	if( x_axis_type == 2){
		for (i = 0; i < n_points; i++) {
			E_MeV_u[i] = E_start + (i/(double)(n_points-1)) * (E_end - E_start);
		}
	} else if( x_axis_type == 1){
		for (i = 0; i < n_points; i++) {
			double logE = log(E_start) + (i/(double)(n_points-1)) * (log(E_end) - log(E_start));
			E_MeV_u[i] = exp(logE);
		}
	} else {
		return EXIT_FAILURE;
	}

	AT_LET_keV_um(  n_points,
	    E_MeV_u,
	    particle_no,
	    material_no,
	    LET_keV_um);

	fprintf(f , "E:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", E_MeV_u[i]);
	}
	fprintf(f, "\n");

	fprintf(f , "LET:");
	for (i = 0; i < n_points; i++) {
		fprintf(f, " %g", LET_keV_um[i]);
	}
	fprintf(f, "\n");
	fclose(f);

	return EXIT_SUCCESS;
}
