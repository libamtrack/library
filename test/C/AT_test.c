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

#include "AT_StoppingPower.h"

int main(){
	const double E_MeV_u   = 270.55;
	const long particle_no = 6012;
	const long material_no = 1;

	double test;

	test = AT_Stopping_Power_keV_um_single(FromFile,
			E_MeV_u,
			particle_no,
			material_no);

	printf("Ergebnis Bethe: %e\n", test);

	test = AT_Stopping_Power_keV_um_single(Bethe,
			E_MeV_u,
			particle_no,
			material_no);

	printf("Ergebnis Bethe: %e\n", test);

	test = AT_Stopping_Power_keV_um_single(PSTAR,
			E_MeV_u,
			particle_no,
			material_no);

	printf("Ergebnis PSTAR: %e\n", test);

	test = AT_Stopping_Power_keV_um_single(ICRU,
			E_MeV_u,
			particle_no,
			material_no);

	printf("Ergebnis ICRU: %e\n", test);

	return EXIT_SUCCESS;
};

