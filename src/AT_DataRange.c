/**
 * @brief Range
 */

/*
 *    AT_DataRange.c
 *    ==============
 *
 *    Created on: 12.11.2010
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

#include "AT_DataRange.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>

double AT_Stopping_Power_Mass_MeV_cm2_g_int( double  E_MeV_u,
		void*   params){
	assert( params != NULL );
	AT_CSDA_range_parameters* int_params = (AT_CSDA_range_parameters*)params;
	double StPow = 0.00001;
        AT_Mass_Stopping_Power_with_no( PSTAR,
		1,
		&E_MeV_u,
		&(int_params->particle_no),
		int_params->material_no,
		&StPow);
 	return (1.0 / (StPow));
}

double AT_CSDA_range_g_cm2_single(	const double 	E_initial_MeV_u,
		const double 	E_final_MeV_u,
		const long 		particle_no,
		const long 		material_no){

	double range_cm2_g              = 0.0;
	AT_CSDA_range_parameters  params;
	params.material_no             = material_no;
	params.particle_no			   = 1001;          // Compute CSDA range for protons, then scale (see below)

	/* Initialize GSL integration workspace */
	gsl_set_error_handler_off();
	gsl_integration_workspace *w1   = gsl_integration_workspace_alloc (10000);
	gsl_function F;
	F.function                      = &AT_Stopping_Power_Mass_MeV_cm2_g_int;
	F.params                        = (void*)&params;

	/* Set integration limits */
	double   lower_lim_m            = E_final_MeV_u;
	double   upper_lim_m            = E_initial_MeV_u;
	double   error;

	/* Perform integration */
	int status      = gsl_integration_qags (        &F,
			lower_lim_m,
			upper_lim_m,
			1.0e-6,
			1.0e-3,
			10000,
			w1,
			&range_cm2_g,
			&error);
	if (status == GSL_EROUND){
#ifndef NDEBUG
		printf("Error in integration of CSDA range from Bethe formula: round-off error.\n");
#endif
		}
	if (status == GSL_ESING){
#ifndef NDEBUG
		printf("Error in integration of CSDA range from Bethe formula: singularity found!\n");
#endif
		}

	gsl_integration_workspace_free (w1);

	/* Scale range for ions Z,A >1 and return value */
	long     Z                      = AT_Z_from_particle_no_single(particle_no);
	long     A                      = AT_A_from_particle_no_single(particle_no);
	return(range_cm2_g * A / (Z * Z));
}

void AT_CSDA_range_g_cm2_multi(	const long    n,
		const double 	E_initial_MeV_u[],
		const double 	E_final_MeV_u[],
		const long 	particle_no[],
		const long 	material_no,
		double          CSDA_range_cm2_g[])
{
	long i;
	for (i = 0; i < n; i++){
		CSDA_range_cm2_g[i] = AT_CSDA_range_g_cm2_single(	E_initial_MeV_u[i],
				E_final_MeV_u[i],
				particle_no[i],
				material_no);
	}
}


double AT_CSDA_range_difference_solver( double  E_final_MeV_u,
		void*   params)
{
	assert( params != NULL );
	AT_CSDA_range_difference_parameters* solver_params = (AT_CSDA_range_difference_parameters*)params;

	double material_range_g_cm2 = AT_CSDA_range_g_cm2_single(	solver_params->E_initial_MeV_u,
			E_final_MeV_u,
			solver_params->particle_no,
			solver_params->material_no);

	return (material_range_g_cm2 - solver_params->range_g_cm2);
}


double AT_CSDA_energy_after_slab_E_MeV_u_single( const double E_initial_MeV_u,
		const long   particle_no,
		const long   material_no,
		const double slab_thickness_m)
{
	double slab_thickness_g_cm2 = (slab_thickness_m * m_to_cm) * AT_density_g_cm3_from_material_no(material_no);
	AT_CSDA_range_difference_parameters params;
	params.E_initial_MeV_u      = E_initial_MeV_u;
	params.particle_no          = particle_no;
	params.material_no          = material_no;
	params.range_g_cm2          = slab_thickness_g_cm2;

	const double  solver_accuracy  =  1e-6;
	// Set lower limit to 1 MeV (= 25 �m range error for protons in water, which should be tolerable) to
	// coincide with coded limits of Bethe formular routine. Otherwise, singularity error will appear
	// TODO: Generalize lower limit for all stopping power sources, e.g. lowest value for tabulated data
	const double  min_possible_value_E_final_MeV_u = 1.0;
	const double  max_possible_value_E_final_MeV_u = E_initial_MeV_u;

	double E_final_MeV_u =  zriddr(AT_CSDA_range_difference_solver,
			(void*)(&params),
			min_possible_value_E_final_MeV_u,
			max_possible_value_E_final_MeV_u,
			solver_accuracy);

	return (E_final_MeV_u);//E_final_MeV_u;
}

void AT_CSDA_energy_after_slab_E_MeV_u_multi( const long n,
		const double E_initial_MeV_u[],
		const long   particle_no[],
		const long   material_no,
		const double slab_thickness_m,
		double E_final_MeV_u[]){

	long i;
	for (i = 0; i < n; i++){
		E_final_MeV_u[i] = AT_CSDA_energy_after_slab_E_MeV_u_single(E_initial_MeV_u[i],
				particle_no[i],
				material_no,
				slab_thickness_m);
	}
}

double AT_WEPL_single(	const double 	E_MeV_u,
		const long 		particle_no,
		const long 		material_no,
		const double    slab_thickness_m){


	double E_final_MeV_u          = AT_CSDA_energy_after_slab_E_MeV_u_single(E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_m);

	double residual_range_g_cm2   = AT_CSDA_range_g_cm2_single( E_final_MeV_u,
			BETHE_LOWER_LIMIT_E_MEV_U,
			particle_no,
			Water_Liquid);

	double residual_range_m       = residual_range_g_cm2 / AT_density_g_cm3_from_material_no(Water_Liquid) / m_to_cm;

	double range_water_m          = AT_CSDA_range_g_cm2_single(E_MeV_u,
			BETHE_LOWER_LIMIT_E_MEV_U,
			particle_no,
			Water_Liquid) / m_to_cm;

	return (range_water_m - residual_range_m) / slab_thickness_m;
}

void AT_WEPL_multi(	const long    n,
		const double 	E_MeV_u[],
		const long 		particle_no[],
		const long 		material_no,
		const double    slab_thickness_m,
		double          WEPL[])
{
	long i;
	for (i = 0; i < n; i++){
		WEPL[i] = AT_WEPL_single(	E_MeV_u[i],
				particle_no[i],
				material_no,
				slab_thickness_m);
	}

}
