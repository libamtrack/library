/**
 * @brief A simple file to check the library and enable debugging
 */

/*
 *    AT_test.c
 *    ===================
 *
 *    Created on: 2009-06-08
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
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
	#include <sys/malloc.h>
#else
	#include <malloc.h>
#endif

#include "config.h"

// CAVE: Remove this by appropriate include file when using other functions from libamtrack
// As it can trigger very strange behaviour such as faulty return values etc.
// #include "AT.h"

int main(){
	 FILE *CSV;
  	 CSV = fopen("FLUKA_DEDX_WATER_76.8eV.txt", "r");

	 if(NULL == CSV) {
		 fprintf(stderr, "Error opening stopping power data file ...\n");
		 return EXIT_FAILURE;
	 }


	 char 	line[1000], test[1000];

	 int  	i = 0, ZZ;
	 float	EE, SS;

	 int	Z[10000];
	 float	E[10000];
	 float  S[10000];

	 while (fgets(line, sizeof line, CSV)) {
	     if (*line == EOF) break;
	     if (isalpha(*line)) continue; /* ignore comment line */
	     	 sscanf(line, "%s", test);
	     	 if ((sscanf(line, "%d %e %e", &ZZ, &EE, &SS) != 3)&&(strlen(test) != 0)&&(strlen(test) != 1)) {
	         /* handle error */
				 printf("Error reading stopping power data (format correct?)...");
				 return EXIT_FAILURE;
			 } else {
			 /* handle variables */
				 Z[i] = ZZ;
				 E[i] = EE;
				 S[i++] = SS;
			 }
	 	 }

	 printf("Read %d lines.\n", i--);


	 int* Z_array = (int*)malloc(i*sizeof(int));
	 float* E_array = (float*)malloc(i*sizeof(float));
	 float* S_array = (float*)malloc(i*sizeof(float));

	 memcpy(Z_array, Z, i*sizeof(int));
	 memcpy(E_array, E, i*sizeof(float));
	 memcpy(S_array, S, i*sizeof(float));

	 for(int j = 0; j < i; j++){
		 printf("line %d: %d %e %e\n", j, Z_array[j], E_array[j], S_array[j]);
	 }

	 return EXIT_SUCCESS;
};

