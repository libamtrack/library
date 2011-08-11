/**
 * @brief Dummy file to enable debugging...
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
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "AT_DataRange.h"

#include <assert.h>

int main(){

	
	const long n = 3;
	const long Z[3] = {1,6,8};
	const long A[3] = {1,12,16};
	const double weight_fraction[3] = {6, 3, 1};
	double average_A = 6666;
	double effective_Z = 6666;
    double exponent = 3.5;
	double el_dens = 0;
    double I_eV = 0;
    double density = 2.0;

	AT_average_A_from_composition( n,
	    A,
	    weight_fraction,
	    &average_A);

	AT_effective_Z_from_composition(n, Z, weight_fraction, exponent, &effective_Z);
	AT_I_eV_from_composition(n, A, Z, weight_fraction, &I_eV);
	AT_electron_density_m3_from_composition(n, density, Z, A, weight_fraction, &el_dens);

	return EXIT_SUCCESS;
};
