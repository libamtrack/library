/*
 * gridSummation.c
 *
 *  Created on: 2009-09-14
 *      Author: grzanka
 */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "SGParticle.h"

int main(){

	long	n = 1;
	long	particle_no = 18;
	float	fluence_cm2 = -1.;
	long	material_no = 1;
	float	LET_MeV_cm2_g = 10*38;
	float	E_MeV_u = 60;

//	printf("Energy = %g\n", E_MeV_u);
//	SGP_LET_MeV_cm2_g(&n, &E_MeV_u, &particle_no,&material_no,&LET_MeV_cm2_g);
	printf("LET = %g [MeV cm2 g]\n", LET_MeV_cm2_g);

	SGP_E_MeV_from_LET(&n,&LET_MeV_cm2_g,&particle_no,&material_no,&E_MeV_u);
	printf("Energy = %g\n", E_MeV_u);

	//SGP_LET_keV_um(&n,&E_MeV_u)
	//E_MeV_u = -E_MeV_u;
	long	RDD_model = 3;
	float	RDD_parameters[] = {1e-9};
	long	ER_model = 1;
	float	ER_parameters[] = {1.};
	long	gamma_model = 5;
	float	gamma_parameters[] = {0.536384,0.03793949,20.0};
	long	N_runs = 10;
	float	fluence_factor = 1.0;
	bool	write_output = true;
	long	nX = 10;
	float	grid_size_m = 1e-6;
	bool	lethal_events_mode = false;
	float	results_f[10];
	float	results_t[10];

	float dose = 0.;
	for( dose = 1; dose < 6; dose += 1.){
		printf("\nDose: %g\n" , dose);
		fluence_cm2 = -dose;

		lethal_events_mode = false;
		SGP_efficiency_grid(	&n,
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
				&fluence_factor,
				&write_output,
				&nX,
				&grid_size_m,
				&lethal_events_mode,
				results_f);

		lethal_events_mode = true;
		SGP_efficiency_grid(	&n,
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
				&fluence_factor,
				&write_output,
				&nX,
				&grid_size_m,
				&lethal_events_mode,
				results_t);

		printf("Survival gamma: %g (f) %g (t)\n" , results_f[3], results_t[3]);
		printf("Survival HCP: %g (f) %g (t)\n" , results_f[2], results_t[2]);

	}
	return 1;
}
