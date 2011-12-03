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

#include "AT_DataStoppingPower.h"


double AT_Bethe_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no){
	return AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_single(	E_MeV_u,
															particle_no,
															material_no,
															-1.0);
}

double AT_ShieldHit_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no){
	// TODO: Later - read-in data from table here
	if(material_no != 1){		/* Water data only */
		return 0.0;
	}

	long Z		= AT_Z_from_particle_no_single(particle_no);
	if( Z > 18){				/* Data for H ... Ar */
		return 0.0;
	}

	if (E_MeV_u < 0.025 || E_MeV_u > 1000){
		return 0.0;
	}

	double	StoppingPower_MeV_cm2_g = 0.0;

	StoppingPower_MeV_cm2_g = AT_get_interpolated_y_from_input_table(
			AT_stopping_power_ShieldHit_table.E_MeV_u_and_stopping_power_total_MeV_cm2_g[0],
			AT_stopping_power_ShieldHit_table.E_MeV_u_and_stopping_power_total_MeV_cm2_g[Z],
			AT_stopping_power_ShieldHit_table.number_of_data_points,
			E_MeV_u);
	return StoppingPower_MeV_cm2_g;
}

double AT_ICRU_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no){
	// TODO: Later - read-in data from table here
	if(material_no != 1){		/* Water data only */
		return 0.0;
	}

	long Z		= AT_Z_from_particle_no_single(particle_no);
	if( Z > 18){				/* Data for H ... Ar */
		return 0.0;
	}

	if (E_MeV_u < 0.025 || E_MeV_u > 1000){
		return 0.0;
	}

	if (Z == 2 && E_MeV_u > 250){   // Data for He only until 250 MeV/u = 1000 MeV !
		return 0.0;
	}

	double	StoppingPower_MeV_cm2_mg = 0.0;
	StoppingPower_MeV_cm2_mg = AT_get_interpolated_y_from_input_table(
			AT_stopping_power_ICRU_table.E_MeV_u_and_stopping_power_total_MeV_cm2_g[0],
			AT_stopping_power_ICRU_table.E_MeV_u_and_stopping_power_total_MeV_cm2_g[Z],
			AT_stopping_power_ICRU_table.number_of_data_points,
			E_MeV_u);
	return StoppingPower_MeV_cm2_mg * 1000;		// ICRU (73) gives stopping power in MeV*cm2/mg !
}

double _AT_Stopping_Power_get_data(const long stopping_power_source_no,
		const long 	    particle_no,
		const long 		material_no,
		AT_stopping_power_tabulated_source_for_given_material_struct ** source_for_given_material){

	assert( AT_stopping_power_tabulated_source.stopping_power_source_data_group != NULL );
	assert( AT_stopping_power_tabulated_source.number_of_sources == STOPPING_POWER_SOURCE_N);
	assert( stopping_power_source_no < AT_stopping_power_tabulated_source.number_of_sources);
	assert( stopping_power_source_no >= 0 );

	AT_stopping_power_tabulated_source_group_for_all_materials_struct * group_for_all_materials =
			(AT_stopping_power_tabulated_source_group_for_all_materials_struct *)AT_stopping_power_tabulated_source.stopping_power_source_data_group[stopping_power_source_no];

	if( group_for_all_materials != NULL ){

		assert( group_for_all_materials->material_no[material_no] == material_no);
		assert( group_for_all_materials->stopping_power_source_no == stopping_power_source_no);

		*source_for_given_material =	(AT_stopping_power_tabulated_source_for_given_material_struct *)group_for_all_materials->stopping_power_source_data[material_no];
		if( *source_for_given_material == NULL ){
			char source_name[STOPPING_POWER_SOURCE_NAME_LENGTH];
			AT_stopping_power_source_model_name_from_number(stopping_power_source_no,source_name);
			char material_name[MATERIAL_NAME_LENGTH];
			AT_material_name_from_number(material_no,material_name);
#ifndef NDEBUG
			printf("Missing data for data source [%s] and material [%s]\n", source_name, material_name);
#endif
			}
	} else {
		char source_name[STOPPING_POWER_SOURCE_NAME_LENGTH];
		AT_stopping_power_source_model_name_from_number(stopping_power_source_no,source_name);
#ifndef NDEBUG
		printf("Missing data for data source [%s]\n", source_name);
#endif
		return -1;
	}
	return -1;

}


double AT_Stopping_Power_data_interpolation(const long stopping_power_source_no,
		const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no){

	AT_stopping_power_tabulated_source_for_given_material_struct * source_for_given_material = NULL;

	_AT_Stopping_Power_get_data( stopping_power_source_no, particle_no, material_no, &source_for_given_material);

	assert( source_for_given_material != NULL );

	const long n = source_for_given_material->number_of_data_points;

	if( source_for_given_material->E_MeV_u_and_stopping_power_total_MeV_cm2_g != NULL){
		double result = AT_get_interpolated_y_from_input_2d_table(
				source_for_given_material->E_MeV_u_and_stopping_power_total_MeV_cm2_g,
				n,
				E_MeV_u);

		return result;
	} else {
		char source_name[STOPPING_POWER_SOURCE_NAME_LENGTH];
		AT_stopping_power_source_model_name_from_number(stopping_power_source_no,source_name);
		char material_name[MATERIAL_NAME_LENGTH];
		AT_material_name_from_number(material_no,material_name);
#ifndef NDEBUG
		printf("Missing data points for data source [%s] and material [%s]\n", source_name, material_name);
#endif
		}
	return -1;
}


double AT_Stopping_Power_MeV_cm2_g_single( const long stopping_power_source_no, const double E_MeV_u, const long particle_no, const long material_no){
	double Stopping_Power_MeV_cm2_g = 0.0;
	if( AT_stopping_power_analytical_source.access_function[stopping_power_source_no] != NULL){
		Stopping_Power_MeV_cm2_g = AT_stopping_power_analytical_source.access_function[stopping_power_source_no](E_MeV_u, particle_no, material_no);
	} else {
		Stopping_Power_MeV_cm2_g = AT_Stopping_Power_data_interpolation(stopping_power_source_no, E_MeV_u, particle_no, material_no);
	}

	// TODO: Move the following Z_eff scaling into separate PSTAR wrapper
	if( (particle_no != PARTICLE_PROTON_NUMBER) && (stopping_power_source_no == PSTAR)){
		double Zeff_ion    =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u, particle_no);
		double Zeff_proton =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u, PARTICLE_PROTON_NUMBER);
		Stopping_Power_MeV_cm2_g *= gsl_pow_2(Zeff_ion / Zeff_proton);
	}

	return Stopping_Power_MeV_cm2_g;
}


void AT_Stopping_Power_MeV_cm2_g_multi( const long stopping_power_source_no, const long number_of_particles, const double E_MeV_u[], const long particle_no[], const long material_no, double Stopping_Power_MeV_cm2_g[]){
	long i;
	for( i = 0 ; i < number_of_particles; i++){
		Stopping_Power_MeV_cm2_g[i] = AT_Stopping_Power_MeV_cm2_g_single(stopping_power_source_no, E_MeV_u[i], particle_no[i], material_no);
	}
}


double AT_Stopping_Power_keV_um_single( const long stopping_power_source_no, const double E_MeV_u, const long particle_no, const long material_no){

	double Stopping_Power_MeV_cm2_g = AT_Stopping_Power_MeV_cm2_g_single( stopping_power_source_no, E_MeV_u, particle_no, material_no);

    double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);

	return Stopping_Power_MeV_cm2_g * material_density_g_cm3 * 0.1;
}


void AT_Stopping_Power_keV_um_multi( const long stopping_power_source_no, const long number_of_particles, const double E_MeV_u[], const long particle_no[], const long material_no, double Stopping_Power_keV_um[]){
	long i;
	for( i = 0 ; i < number_of_particles; i++){
		Stopping_Power_keV_um[i] = AT_Stopping_Power_keV_um_single(stopping_power_source_no, E_MeV_u[i], particle_no[i], material_no);
	}
}


double AT_Energy_MeV_u_from_Stopping_Power_single( const long stopping_power_source_no,
		const double Stopping_Power_MeV_cm2_g,
		const long particle_no,
		const long material_no){

	// should work for energies from 0.5 MeV till 1000 MeV

	AT_stopping_power_tabulated_source_for_given_material_struct * source_for_given_material = NULL;

	_AT_Stopping_Power_get_data( stopping_power_source_no, particle_no, material_no, &source_for_given_material);

	assert( source_for_given_material != NULL );

	const long n = source_for_given_material->number_of_data_points;

	if( source_for_given_material->E_MeV_u_and_stopping_power_total_MeV_cm2_g != NULL){

		long lowest_index = locate_index_in_2d_table( source_for_given_material->E_MeV_u_and_stopping_power_total_MeV_cm2_g, 0, n-1, 4.9, 0 );
		long highest_index = locate_index_in_2d_table( source_for_given_material->E_MeV_u_and_stopping_power_total_MeV_cm2_g, 0, n-1, 1000.0, 0 );

		double tab[highest_index-lowest_index+1][2];
		long i ;
		for( i = lowest_index-1 ; i < highest_index ; i++ ){
			double E_MeV_u                                =  source_for_given_material->E_MeV_u_and_stopping_power_total_MeV_cm2_g[i][0];
			double stopping_power_total_MeV_cm2_g_proton  =  source_for_given_material->E_MeV_u_and_stopping_power_total_MeV_cm2_g[i][1];
			tab[i - lowest_index + 1][0] = E_MeV_u;
			if( particle_no != PARTICLE_PROTON_NUMBER ){
				double Zeff_ion    =  AT_effective_charge_from_E_MeV_u_single( E_MeV_u, particle_no);
				double Zeff_proton =  AT_effective_charge_from_E_MeV_u_single( E_MeV_u, PARTICLE_PROTON_NUMBER);
				tab[i - lowest_index + 1][1] = stopping_power_total_MeV_cm2_g_proton * gsl_pow_2(Zeff_ion / Zeff_proton);
			} else {
				tab[i - lowest_index + 1][1] = stopping_power_total_MeV_cm2_g_proton;
			}
		}

		if( (Stopping_Power_MeV_cm2_g < tab[highest_index - lowest_index][1] ) || (Stopping_Power_MeV_cm2_g > tab[0][1]) ){
#ifndef NDEBUG
			printf("Only energy region 5-1000 MeV supported\n");
#endif
			return -1;
		}

		double result = AT_get_interpolated_x_from_input_2d_table(
				(const double (*)[2])tab,
				0,
				highest_index - lowest_index,
				Stopping_Power_MeV_cm2_g);

		return result;
	} else {
		char source_name[STOPPING_POWER_SOURCE_NAME_LENGTH];
		AT_stopping_power_source_model_name_from_number(stopping_power_source_no,source_name);
		char material_name[MATERIAL_NAME_LENGTH];
		AT_material_name_from_number(material_no,material_name);
#ifndef NDEBUG
		printf("Missing data points for data source [%s] and material [%s]\n", source_name, material_name);
#endif
		}
	return -1;
}


int AT_stopping_power_source_model_name_from_number( const long source_no, char* source_name){

	assert( source_no >= 0);
	assert( source_no < STOPPING_POWER_SOURCE_N);
	assert(AT_stopping_power_sources.stopping_power_source_no[source_no] == source_no);
	if( source_no < 0)
		return -1;
	if( source_no >= STOPPING_POWER_SOURCE_N)
		return -1;	
    strcpy(source_name, AT_stopping_power_sources.stopping_power_source_name[source_no]);
    return AT_Success;
}


long AT_stopping_power_source_model_number_from_name( const char* source_name ){

	long  match;
	const long n_tmp = 1;

	assert( source_name != NULL);

	find_elements_char(  &source_name,
			n_tmp,
			AT_stopping_power_sources.stopping_power_source_name,
			AT_stopping_power_sources.n,
			&match);

	return match;
}


double AT_Stopping_Power_Bethe_Number(	const double 	E_MeV_u,
										const long 		particle_no,
										const long 		material_no,
										const double	E_restricted_keV)
{
	  const double beta2 	= gsl_pow_2(AT_beta_from_E_single(E_MeV_u));
	  const double I_eV		= AT_I_eV_from_material_no(material_no);
	  const double I_MeV	= I_eV * 1e-6;
	  double Wm_MeV			= AT_max_relativistic_E_transfer_MeV_single(E_MeV_u);

	  /* Restricted stopping number requested? */
	  bool restricted = false;
	  if(	(E_restricted_keV > 0.0) && (E_restricted_keV / 1000.0 < Wm_MeV))
		  restricted = true;

	  /* First part of stopping number */
	  double SN11			=  2.0 * electron_mass_MeV_c2 * beta2 / (1.0 - beta2);
	  assert( I_MeV > 0. );
	  SN11					/= I_MeV;

	  if(	restricted){
		  Wm_MeV				= E_restricted_keV * 1e-3;
	  }
	  double SN12			= Wm_MeV / I_MeV;

	  /* Second part of stopping number */
	  double SN2			= beta2;
	  if(	restricted){
		  SN2					/= 2;
		  SN2					+= (1.0 - beta2) * Wm_MeV / (4.0 * electron_mass_MeV_c2);
	  }

	  /* Third part of stopping number (density correction following Sternheimer, 1971) */
	  double delta				= 0.0;
	  long   phase              = AT_phase_from_material_no(material_no);
	  if( (phase =! phase_undefined) ){
		  double kinetic_variable	= AT_kinetic_variable_single(E_MeV_u);
		  double plasma_energy_J 	= AT_plasma_energy_J_from_material_no(material_no);
		  double I_J				= I_MeV * MeV_to_J;

		  double C					= 1.0 + 2.0 * log(I_J / plasma_energy_J);

		  // Find x_0 and x_1 dependent on phase, I-value and C
		  double x_0 = 0.0;
		  double x_1 = 0.0;
		  if( phase == phase_condensed){
			  if(I_eV < 100){
				  x_1	= 2.0;
				  if(C <= 3.681){
					  x_0 	= 0.2;
				  }else{
					  x_0	= 0.326 * C - 1.0;
				  }
			  }else{ // I_eV >= 100
				  x_1	= 3.0;
				  if(C <= 5.215){
					  x_0	= 0.2;
				  }else{
					  x_0	= 0.326 * C - 1.5;
				  }
			  }
		  }else{ // gaseous material
			  x_0	= 0.326 * C - 2.5;
			  x_1	= 5.0;
			  if(C < 10.0){
				  x_0	= 1.6;
				  x_1	= 4.0;
			  }
			  if(C >= 10.0 && C < 10.5){
				  x_0	= 1.7;
				  x_1	= 4.0;
			  }
			  if(C >= 10.5 && C < 11.0){
				  x_0	= 1.8;
				  x_1	= 4.0;
			  }
			  if(C >= 11.0 && C < 11.5){
				  x_0	= 1.9;
				  x_1	= 4.0;
			  }
			  if(C >= 11.5 && C < 12.25){
				  x_0	= 2.0;
				  x_1	= 4.0;
			  }
			  if(C >= 12.25 && C < 13.804){
				  x_0	= 2.0;
				  x_1	= 5.0;
			  }
		  }

		  double x_a				= C / 4.606;
		  double m					= 3.0;
		  double a					= 4.606 * (x_a - x_0) / pow(x_1 - x_0, m);

		  if( kinetic_variable >= x_0 && kinetic_variable <= x_1){
			  delta						= 4.606 * kinetic_variable - C + a * pow(x_1 - kinetic_variable, m);
		  }
		  if( kinetic_variable > x_1){
			  delta						= 4.606 * kinetic_variable - C;
		  }
	  }
	  double SN3		= delta;

	  /* Forth part of stopping number (shell correction) TODO: implement */


	  assert( SN11 > 0. );
	  assert( SN12 > 0. );

	  return (0.5 * log(SN11 * SN12) - SN2 - SN3);
}


double AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_single(	const double E_MeV_u,
														const long particle_no,
														const long material_no,
														const double E_restricted_keV)
{
	double result = 0;
	/* Compute only above 1.0 MeV, otherwise theory is too wrong
	 * below return zero */
	// TODO: Find smarter criterion because this may cause problems in the code (as it did
	// TODO: with the inappropriatly set lower limit for CSDA range integration (was 0, now 1.0 MeV)
	if(E_MeV_u >= BETHE_LOWER_LIMIT_E_MEV_U){
		const double beta2 		= gsl_pow_2(AT_beta_from_E_single(E_MeV_u));
		assert( beta2 > 0);
		const double Z			= AT_average_Z_from_material_no(material_no);
		const double A			= AT_average_A_from_material_no(material_no);
		assert( A > 0 );
		const double z_eff		= AT_effective_charge_from_E_MeV_u_single(E_MeV_u, particle_no);
		const double k_MeV_cm2_g	= 0.307075;												// ICRU49, p.6, after Cohen and Taylor (1986)
		const double SN			= AT_Stopping_Power_Bethe_Number(	E_MeV_u,
				particle_no,
				material_no,
				E_restricted_keV);
		result = k_MeV_cm2_g * (Z / A) * gsl_pow_2(z_eff) * SN / (beta2);
	}
    return result;

}


void AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_multi(	const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const double E_restricted_keV,
		double Mass_Stopping_Power_MeV_cm2_g[])
{
	long i;
	for (i = 0; i < n; i++){
		Mass_Stopping_Power_MeV_cm2_g[i]	=	AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_single(	E_MeV_u[i],
																								particle_no[i],
																								material_no,
																								E_restricted_keV);
	}
}

int AT_Rutherford_SDCS(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const long n,
		const double T_MeV[],
		double dsdT_m2_MeV[]){

	// Get particle data
	long I 			= AT_nuclear_spin_from_particle_no_single(particle_no);
	double Z_eff  	= AT_effective_charge_from_E_MeV_u_single(E_MeV_u, particle_no);
	long Z			= AT_Z_from_particle_no_single(particle_no);
	long A			= AT_A_from_particle_no_single(particle_no);

	// Get material Z
	double Z_material = AT_average_Z_from_material_no(material_no);

	// We follow the formulation of Sawakuchi, 2007
	double particle_mass_MeV_c2	= Z * proton_mass_MeV_c2 + (A-Z) * neutron_mass_MeV_c2;
	double total_E_MeV_u		= particle_mass_MeV_c2 / A + E_MeV_u;
	double Q_c_MeV_c2			= particle_mass_MeV_c2 * particle_mass_MeV_c2 / electron_mass_MeV_c2;
	double K_MeV_m2  			= 2 * M_PI * pow(classical_electron_radius_m, 2) * electron_mass_MeV_c2;

	double beta				= AT_beta_from_E_single(E_MeV_u);
	double beta2			= beta * beta;
	double T_max_MeV		= AT_max_E_transfer_MeV_single(E_MeV_u);

	double term_0, term_1, term_2;

	long i;
	for (i = 0; i < n; i++){
		if(T_MeV[i] > T_max_MeV){
			dsdT_m2_MeV[i]		= 0.0;
		}else{
			term_0				= K_MeV_m2 * Z_material * Z_eff * Z_eff / (beta2 * T_MeV[i] * T_MeV[i]);
			term_1				= 1.0 - beta2 * T_MeV[i] / T_max_MeV;
			if(I == 0.0){
				dsdT_m2_MeV[i]			= term_0 * term_1;
			}
			if(I == 0.5){
				term_2				= T_MeV[i] * T_MeV[i] / (2.0 * total_E_MeV_u * total_E_MeV_u);
				dsdT_m2_MeV[i]			= term_0 * (term_1 + term_2);
			}
			if(I == 1.0){
				term_2				=  (1.0 + T_MeV[i] / (2.0 * Q_c_MeV_c2));
				term_2				*= T_MeV[i] * T_MeV[i] / (3.0 * total_E_MeV_u * total_E_MeV_u);
				term_2				+= term_1 * (1.0 + T_MeV[i] / (3.0 * Q_c_MeV_c2));
				dsdT_m2_MeV[i]			= term_0 * term_2;
			}
		}
	}

	return EXIT_SUCCESS;
}
