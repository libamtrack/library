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

#include "AT_StoppingPower.h"
#include "AT_DataRange.h"

int main() {

    double Test = AT_CSDA_range_g_cm2_single(270.55,
                                             0.1,
                                             6012,
                                             1);

    printf("Range of 270.55 C-12 in water: %3.2f cm\n\n", Test);

    const double E_MeV_u[3] = {1, 10, 100};
    const long particle_no[3] = {6012, 6012, 6012};
    const long material_no = 2;

    double test[3];

    AT_Mass_Stopping_Power("Bethe",
                           3,
                           E_MeV_u,
                           particle_no,
                           material_no,
                           test);
    printf("Bethe: %e, %e, %e\n", test[0], test[1], test[2]);

    AT_Mass_Stopping_Power("PSTAR",
                           3,
                           E_MeV_u,
                           particle_no,
                           material_no,
                           test);
    printf("PSTAR: %e, %e, %e\n", test[0], test[1], test[2]);

    AT_Mass_Stopping_Power("ICRU",
                           3,
                           E_MeV_u,
                           particle_no,
                           material_no,
                           test);
    printf("ICRU: %e, %e, %e\n", test[0], test[1], test[2]);

    return EXIT_SUCCESS;
};

