/**
 * @brief A simple file to check the library...
 */

/*
 *    AT_test.c
 *    ===================
 *
 *    Created on: 2011-08-12
 *    Creator: greilich
 *
 *    Copyright 2006, 2011 The libamtrack team
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
#include <assert.h>

#include "AT_Algorithms_CPP.h"

int main(){

	const double a = 2, b = 3;

	double result = AT_test_fun(a, b);

	if(result == a + b){
		printf("libamtrack is working.\n");
		return EXIT_SUCCESS;
	}else{
		printf("libamtrack is NOT working properly.\n");
		return EXIT_FAILURE;
	}
};

