/**
 * @brief Stopping Power
 */

/*
 *    AT_DataStoppingPower.c
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

#include "AT_StoppingPower.h"

/** ------ FUNCTIONS ------*/

/**
 * Main function to retrieve stopping powers
 *
 */
int AT_Mass_Stopping_Power( const char stopping_power_source[],
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_MeV_cm2_g[]){

	if (strcmp(stopping_power_source, "Bethe") == 0){
		AT_Mass_Stopping_Power_with_no( Bethe,
				n,
				E_MeV_u,
				particle_no,
				material_no,
				stopping_power_MeV_cm2_g);
		return AT_Success;
	}

	if (strcmp(stopping_power_source, "PSTAR") == 0){
		AT_Mass_Stopping_Power_with_no( PSTAR,
				n,
				E_MeV_u,
				particle_no,
				material_no,
				stopping_power_MeV_cm2_g);
		return AT_Success;
	}

	if (strcmp(stopping_power_source, "ICRU") == 0){
		AT_Mass_Stopping_Power_with_no( ICRU,
				n,
				E_MeV_u,
				particle_no,
				material_no,
				stopping_power_MeV_cm2_g);
		return AT_Success;
	}

	int result = AT_stopping_power_functions.function[FromFile](n,
			E_MeV_u,
			particle_no,
			material_no,
			stopping_power_source,
			stopping_power_MeV_cm2_g);

	return result;
}

int AT_Stopping_Power( const char stopping_power_source[],
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_keV_um[]){
	int result = AT_Mass_Stopping_Power( stopping_power_source,
			n,
			E_MeV_u,
			particle_no,
			material_no,
			stopping_power_keV_um);

	long i;
	double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);
	for(i = 0; i < n; i++){
		stopping_power_keV_um[i] *= material_density_g_cm3 / 10.0;
	}

	return (result);
}

int AT_Mass_Stopping_Power_with_no( const long stopping_power_source_no,
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_MeV_cm2_g[]){

	assert( stopping_power_source_no < STOPPING_POWER_SOURCE_N);
	assert( stopping_power_source_no >= 0 );

	int result = AT_stopping_power_functions.function[stopping_power_source_no](n,
			E_MeV_u,
			particle_no,
			material_no,
			NULL,
			stopping_power_MeV_cm2_g);

	return result;
}

int AT_Stopping_Power_with_no( const long stopping_power_source_no,
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_keV_um[]){

	int result = AT_Mass_Stopping_Power_with_no(stopping_power_source_no,
			n,
			E_MeV_u,
			particle_no,
			material_no,
			stopping_power_keV_um);

	long i;
	double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);
	for(i = 0; i < n; i++){
		stopping_power_keV_um[i] *= material_density_g_cm3 / 10.0;
	}

	return (result);
}



double AT_Energy_MeV_u_from_Stopping_Power_single( const long stopping_power_source_no,
		const double Stopping_Power_MeV_cm2_g,
		const long particle_no,
		const long material_no){

//	// should work for energies from 0.5 MeV till 1000 MeV
//
//	PSTAR_data_for_material_struct * source_for_given_material = NULL;
//
//	_AT_Stopping_Power_get_data( stopping_power_source_no, particle_no, material_no, &source_for_given_material);
//
//	assert( source_for_given_material != NULL );
//
//	const long n = source_for_given_material->number_of_data_points;
//
//	if( source_for_given_material->energy_and_stopping_power != NULL){
//
//		long lowest_index = locate_index_in_2d_table( source_for_given_material->energy_and_stopping_power, 0, n-1, 4.9, 0 );
//		long highest_index = locate_index_in_2d_table( source_for_given_material->energy_and_stopping_power, 0, n-1, 1000.0, 0 );
//
//		double tab[highest_index-lowest_index+1][2];
//		long i ;
//		for( i = lowest_index-1 ; i < highest_index ; i++ ){
//			double E_MeV_u                                =  source_for_given_material->energy_and_stopping_power[i][0];
//			double stopping_power_total_MeV_cm2_g_proton  =  source_for_given_material->energy_and_stopping_power[i][1];
//			tab[i - lowest_index + 1][0] = E_MeV_u;
//			if( particle_no != PARTICLE_PROTON_NUMBER ){
//				double Zeff_ion    =  AT_effective_charge_from_E_MeV_u_single( E_MeV_u, particle_no);
//				double Zeff_proton =  AT_effective_charge_from_E_MeV_u_single( E_MeV_u, PARTICLE_PROTON_NUMBER);
//				tab[i - lowest_index + 1][1] = stopping_power_total_MeV_cm2_g_proton * gsl_pow_2(Zeff_ion / Zeff_proton);
//			} else {
//				tab[i - lowest_index + 1][1] = stopping_power_total_MeV_cm2_g_proton;
//			}
//		}
//
//		if( (Stopping_Power_MeV_cm2_g < tab[highest_index - lowest_index][1] ) || (Stopping_Power_MeV_cm2_g > tab[0][1]) ){
//#ifndef NDEBUG
//			printf("Only energy region 5-1000 MeV supported\n");
//#endif
//			return -1;
//		}
//
//		double result = AT_get_interpolated_x_from_input_2d_table(
//				(const double (*)[2])tab,
//				0,
//				highest_index - lowest_index,
//				Stopping_Power_MeV_cm2_g);
//
//		return result;
//	} else {
//		char source_name[STOPPING_POWER_SOURCE_NAME_LENGTH];
//		AT_stopping_power_source_model_name_from_number(stopping_power_source_no,source_name);
//		char material_name[MATERIAL_NAME_LENGTH];
//		AT_material_name_from_number(material_no,material_name);
//#ifndef NDEBUG
//		printf("Missing data points for data source [%s] and material [%s]\n", source_name, material_name);
//#endif
//		}
	return -1;
}



