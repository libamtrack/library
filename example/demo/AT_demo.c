/**
 * @brief Demo
 */

/*
 *    AT_demo.c
 *    ===================
 *
 *    Created on: 2010-10-06
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

#include "AT_PhysicsRoutines.h"

int main( int argc, char* argv[]){

	if( argc != 1){
		printf("Usage: %s\n", argv[0]);
		return EXIT_FAILURE;
	}

	double E_MeV_u = 60;
	double beta = AT_beta_from_E_single( E_MeV_u );
	printf("Relative speed of particle with energy %4.2f is equal %1.3f\n", E_MeV_u, beta);

	void AT_Vavilov_energy_loss_distribution( const long n,
											  const double energy_loss_keV[],
											  const double E_MeV_u,
											  const long particle_no,
											  const long material_no,
											  const double slab_thickness_um,
											  double fDdD[]);

	const long n = 100;
	double eloss_keV[100];
	double fDdD[100];

	for( int i = 0 ; i < 100 ; i++){
		eloss_keV[i] = (double)i*10;
		fDdD[i] = 0.0;
	}

	AT_Vavilov_energy_loss_distribution( n,
	eloss_keV,
	150.0,
	1001,
	1,
	1000.0,
	fDdD);

	for( int i = 0 ; i < n ; i++){
		printf("%g %g\n", eloss_keV[i], fDdD[i]);
	}


	return EXIT_SUCCESS;
}
