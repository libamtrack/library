/**
 * @brief AT_fluence_from_dose wrapper
 */

/*
 *    AT_fluence_from_dose.c
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
	char Text[600];

	long n = 0;
	double D_Gy[500];
	double E_MeV_u[500];
	long particle_no[500];
	double E_MeV_u_single = 50.;
	long particle_no_single = PARTICLE_PROTON_NUMBER;
	long material_no = Water_Liquid;

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "E_MeV_u")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			E_MeV_u_single = atof(token);
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
		if (strstr(Text, "dose_input")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			char* nToken;
			nToken = strtok(token, " ");
			do {
				D_Gy[n] = atof(nToken);
				n++;
			} while ((nToken = strtok(NULL, " ")));
		}
	}

	int i;
	for( i = 0 ; i < n ; i++){
		E_MeV_u[i] = E_MeV_u_single;
		particle_no[i] = particle_no_single;
	}

	double fluence_cm2[500];

	AT_fluence_cm2_from_dose_Gy(  n,
	    E_MeV_u,
	    particle_no,
	    D_Gy,
	    material_no,
	    fluence_cm2);

	fprintf(f , "fluence:");
	for (i = 0; i < n; i++) {
		fprintf(f, " %g", fluence_cm2[i]);
	}
	fprintf(f, "\n");
	fclose(f);

	return EXIT_SUCCESS;
}
