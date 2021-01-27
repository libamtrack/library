/**
 * @brief Investigation of issue #120
 */

/*
 *    AT_test.c
 *    ===================
 *
 *    Created on: 2021-01-27
 *    Creator: Mateusz Baran
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
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
	#include <sys/malloc.h>
#else
	#include <malloc.h>
#endif

#include "AT_RDD.h"
#include "AT_StoppingPower.h"
#include "AT_PhysicsRoutines.h"
#include "AT_DataRange.h"

int main( int argc, char* argv[]){

	printf("Evaluated doses:\n");

	const int N = 2;
	double r_m_tab[] = {1e-13, 1e-8};
    double E_MeV_u = 150.0;
    long particle_no = AT_particle_no_from_particle_name_single("1H");
    long material_no = Water_Liquid;
    long rdd_model = 6;
    double KatzPoint_r_min_m = 1e-10;
    double a0_m = 1e-8;
    double d_min_Gy = 1e-80;
    double rdd_parameter[] = {KatzPoint_r_min_m, a0_m, d_min_Gy, 0.0};
    long er_model = ER_Waligorski;
    long stopping_power_source_no = 2;
    double D_RDD_Gy[2] = {0.0};

    long res = AT_D_RDD_Gy( N,
                     r_m_tab,
                     E_MeV_u,
                     particle_no,
                     material_no,
                     rdd_model,
                     rdd_parameter,
                     er_model,
                     stopping_power_source_no,
                     D_RDD_Gy);

    for(int i=0; i<N; ++i)
    {
        printf("%e : %e\n", r_m_tab[i], D_RDD_Gy[i]);
    }

	return EXIT_SUCCESS;
}
