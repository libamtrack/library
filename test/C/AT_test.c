/**
 * @file
 * @brief Dummy file to enable debugging, to be changed by the user anyway they like.
 */

/*
 *    AT_test.c
 *    ===================
 *
 *    Created on: 2009-06-08
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
#include "AT_Wrapper_R.h"
#include "AT_DataMaterial.h"
#include "AT_DataLET.h"

int main(){

//	const int		n							=	1;
//	const float		E_MeV_u[]					=	{10};
//	const int		particle_no[]				=	{1001};
//	const float	fluence_cm2_or_dose_Gy[]	=	{-1.0};
//	const int		material_no					=	1;
//	const int		rdd_model					=	3;
//	const float	rdd_parameters[]			=	{50e-9};
//	const int		er_model					=	3;
//	const int		gamma_model					=	2;
//	const float	gamma_parameters[]			=	{1,10,1,1,0};
//	const float		saturation_cross_section_factor = 1.0;
//	const int		write_output				=	0;
//	float			results[10];
//
//	AT_run_IGK_method_R(  &n,
//	    E_MeV_u,
//	    particle_no,
//	    fluence_cm2_or_dose_Gy,
//	    &material_no,
//	    &rdd_model,
//	    rdd_parameters,
//	    &er_model,
//	    &gamma_model,
//	    gamma_parameters,
//	    &saturation_cross_section_factor,
//		&write_output,
//	    results);


//	const long		n				= 7;
//	const double 	E_MeV_u[]		= {1, 10, 100, 126, 10000, 100, 10000};
//	double			gamma[7], momentum_MeV_c_u[7], E_MeV_u_back[7];
//
//	AT_gamma_from_E(	n,
//						E_MeV_u,
//						gamma);
//	AT_momentum_from_E_MeV_c_u(	n,
//								E_MeV_u,
//								momentum_MeV_c_u);
//	AT_E_MeV_u_from_momentum(	n,
//								momentum_MeV_c_u,
//								E_MeV_u_back);

	char * test_name = (char*)calloc(MATERIAL_NAME_LENGTH, sizeof(char));
	const long test_no = 1;

	AT_material_name_from_number(test_no, test_name);

	free(test_name);
	return 0;
};
