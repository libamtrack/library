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

double AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_int( double  E_MeV_u,
    void*   params){
  assert( params != NULL );
  AT_CSDA_range_Bethe_parameters* int_params = (AT_CSDA_range_Bethe_parameters*)params;

  const double E_restricted_keV = 0.0;
  double  StPow = AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_single(	E_MeV_u,
			int_params->particle_no,
			int_params->material_no,
			E_restricted_keV);
  return (1.0 / StPow);
}

double AT_CSDA_range_Bethe_cm2_g(	const double 	E_MeV_u,
		const long 		particle_no,
		const long 		material_no){

	double range_cm2_g              = 0.0;
	AT_CSDA_range_Bethe_parameters  params;
	params.material_no             = material_no;
	params.particle_no			   = particle_no;

	/* Initialize GSL integration workspace */
	gsl_set_error_handler_off();
	gsl_integration_workspace *w1   = gsl_integration_workspace_alloc (10000);
	gsl_function F;
	F.function                      = &AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_int;
	F.params                        = (void*)&params;

	/* Set integration limits */
	double   lower_lim_m            = 0.0;
	double   upper_lim_m            = E_MeV_u;
	double   error;

	/* Perform integration */
	int status      = gsl_integration_qags (        &F,
			lower_lim_m,
			upper_lim_m,
			1e-3,
			1e-3,
			10000,
			w1,
			&range_cm2_g,
			&error);
	if (status == GSL_EROUND || status == GSL_ESING){
		printf("Error in integration of CSDA range from Bethe formula!\n");
	}

	gsl_integration_workspace_free (w1);

	return(range_cm2_g);
}

double AT_WEPL_Bethe(	const double 	E_MeV_u,
		const long 		particle_no,
		const long 		material_no){
	const long material_no_water = Water_Liquid;
	return(AT_CSDA_range_Bethe_cm2_g(	E_MeV_u,
			particle_no,
			material_no) / AT_CSDA_range_Bethe_cm2_g(	E_MeV_u,
					particle_no,
					material_no_water));

};
