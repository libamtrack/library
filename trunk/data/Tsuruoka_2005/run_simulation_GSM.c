/**
* @file
* @brief ...
*/

/*
 *    run_simulation_GSM.c
 *    ==============
 *
 *    A simple example of grid summation algorithm
 *    Created on: 20.09.2009
 *    Author: grzanka
 *
 *    Copyright 2006, 2009 Steffen Greilich / the libamtrack team
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

#include "AmTrack.h"

int main(){

  long	n = 1;
  long	particle_no = 18;
  float	fluence_cm2 = -1.;
  long	material_no = 1;
  float	LET_MeV_cm2_g = 10*38;
  float	E_MeV_u;

  printf("LET = %g [MeV cm2 g]\n", LET_MeV_cm2_g);

  // convert LET to energy
  AT_E_MeV_from_LET(&n,&LET_MeV_cm2_g,&particle_no,&material_no,&E_MeV_u);
  printf("Energy = %g\n", E_MeV_u);

  long	RDD_model = 3;
  float	RDD_parameters[] = {1e-9};
  long	ER_model = 1;
  float	ER_parameters[] = {1.};
  long	gamma_model = 5;
  float	gamma_parameters[] = {0.536384,0.03793949,20.0};
  long	N_runs = 10;
  long  N2 = 10;
  float	fluence_factor = 1.0;
  bool	write_output = true;
  long	nX = 10;
  float	grid_size_m = 1e-6;
  bool	lethal_events_mode = true;
  float	result_GSM[10];

  float dose = 0.;
  for( dose = 1; dose < 6; dose += 1.){
    printf("\nDose: %g\n" , dose);
    fluence_cm2 = -dose;

    AT_GSM(	&n,
        &E_MeV_u,
        &particle_no,
        &fluence_cm2,
        &material_no,
        &RDD_model,
        RDD_parameters,
        &ER_model,
        ER_parameters,
        &gamma_model,
        gamma_parameters,
        &N_runs,
        &N2,
        &fluence_factor,
        &write_output,
        &nX,
        &grid_size_m,
        &lethal_events_mode,
        result_GSM);

    printf("Survival photon: %g\n" , result_GSM[3]);
    printf("Survival   ions: %g\n" , result_GSM[2]);

  }
  return EXIT_SUCCESS;
}
