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

	long n = 0;
	double E_MeV_u[500];
	double LET_keV_um[500];
	long material_no;
	long particle_no[500];
	long particle_no_single;

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "E_input")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			char* nToken;
			nToken = strtok(token, " ");
			do {
				E_MeV_u[n] = atof(nToken);
				n++;
			} while ((nToken = strtok(NULL, " ")));
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
	for( i = 0 ; i < n ; i++){
		particle_no[i] = particle_no_single;
	}

	AT_LET_keV_um(  n,
	    E_MeV_u,
	    particle_no,
	    material_no,
	    LET_keV_um);

	char str[] = { "LET_output:" };
	for (i = 0; i < n; i++) {
		char text[1024];
		sprintf(text, " %f", LET_keV_um[i]);
		strcat(str, text);
	}
	strcat(str, "\n");
	fputs(str, f);

	if (f) {
		fclose(f);
	}

	return EXIT_SUCCESS;
}
