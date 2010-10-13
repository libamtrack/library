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
#include "AmTrack.h"

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Wrong parameters");
		return EXIT_FAILURE;
	}
	char *path = argv[1];
	char Text[600];

	long n = 0;
	double beta[500];

	FILE *f;
	fflush(stdin);
	f = fopen(path, "a+");

	if (f == NULL) {
		return EXIT_FAILURE;
	}

	while (fgets(Text, sizeof(Text), f) != 0) {
		if (strstr(Text, "beta_input")) {
			strtok(Text, ":");
			char* token = strtok(NULL, ":");
			char* nToken;
			nToken = strtok(token, " ");
			do {
				beta[n] = atof(nToken);
				n++;
			} while ((nToken = strtok(NULL, " ")));
		}
	}

	double E[500];

	AT_E_from_beta(n, beta, E);
	char str[] = { "E_output:" };
	int i;
	for (i = 0; i < n; i++) {
		char text[1024];
		sprintf(text, " %f", E[i]);
		strcat(str, text);
	}
	strcat(str, "\n");
	fputs(str, f);

	if (f) {
		fclose(f);
	}

	return EXIT_SUCCESS;
}
