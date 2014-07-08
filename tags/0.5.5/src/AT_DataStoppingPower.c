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
	return AT_Bethe_energy_loss_MeV_cm2_g_single(	E_MeV_u,
													particle_no,
													material_no,
													-1.0,
													true);
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

double AT_FLUKA_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no){
	if(material_no != 1){		/* Water data only */
		return 0.0;
	}

	long Z		= AT_Z_from_particle_no_single(particle_no);
	if( Z > 8){				/* Data for H ... O */
		return 0.0;
	}

	if (E_MeV_u < 0.0126 || E_MeV_u > 1000){
		return 0.0;
	}

	double	StoppingPower_MeV_cm2_g = 0.0;
	StoppingPower_MeV_cm2_g = AT_get_interpolated_y_from_input_table(
			AT_stopping_power_FLUKA_table.E_MeV_u_and_stopping_power_total_MeV_cm2_g[0],
			AT_stopping_power_FLUKA_table.E_MeV_u_and_stopping_power_total_MeV_cm2_g[Z],
			AT_stopping_power_FLUKA_table.number_of_data_points,
			E_MeV_u);

	return StoppingPower_MeV_cm2_g;	
}

double AT_ATIMA_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no){
	
	if(material_no != 6){		/* LiF data only */
		return 0.0;
	}

	long Z		= AT_Z_from_particle_no_single(particle_no);
	if( Z > 92){				/* Data for H ... U */
		return 0.0;
	}

	if (E_MeV_u < 0.01 || E_MeV_u > 10000){
		return 0.0;
	}

	double	StoppingPower_MeV_cm2_g = 0.0;
	StoppingPower_MeV_cm2_g = AT_get_interpolated_y_from_input_table(
			AT_stopping_power_ATIMA_table.E_MeV_u_and_stopping_power_total_MeV_cm2_g[0],
			AT_stopping_power_ATIMA_table.E_MeV_u_and_stopping_power_total_MeV_cm2_g[Z],
			AT_stopping_power_ATIMA_table.number_of_data_points,
			E_MeV_u);
	return StoppingPower_MeV_cm2_g;	
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





void get_table_value(
    const long    n,
    const double  x[],
    const long    subset_no,
    const double  x_table[],
    const double  y_table[],
    double        y[])
{
  // first: find those PSTAR entries that match the material number
  bool*    matches    =  (bool*)calloc(AT_PSTAR_Data.n, sizeof(bool));
  is_element_int(    subset_no,
      AT_PSTAR_Data.material_no,
      AT_PSTAR_Data.n,
      matches);

  long    n_matches  = 0;
  long    i;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      n_matches++;
    }
  }

  // allocate vectors for extracted LET entries
  double*  x_c  =  (double*)calloc(n_matches, sizeof(double));
  double*  y_c  =  (double*)calloc(n_matches, sizeof(double));

  // and get the values
  long     j  = 0;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      x_c[j]  = x_table[i];
      y_c[j]  = y_table[i];
      j++;
    }
  }
  for (i = 0; i < n; i++){
    // Get proton-LET for scaled energy from table E, L using linear interpolation (to be done on logscale TODO)
    y[i] = AT_get_interpolated_y_from_input_table(x_c, y_c, n_matches, x[i]);
  }

  free(x_c);
  free(y_c);
  free(matches);
}

double AT_CSDA_range_g_cm2_single(   const double  E_MeV_u,
    const long    particle_no,
    const long    material_no)
{
  double CSDA_range_g_cm2;
  const long number_of_particles  =  1;
  get_table_value(number_of_particles, &E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.range_cdsa_g_cm2, &CSDA_range_g_cm2);

  // Conversion CSDA_proton => CSDA_ion
  long Z = AT_Z_from_particle_no_single(particle_no);
  long A = AT_A_from_particle_no_single(particle_no);

  return CSDA_range_g_cm2  * (double)(A)/(double)(Z*Z);
}



void AT_CSDA_range_g_cm2(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    double        CSDA_range_g_cm2[])
{
  get_table_value(number_of_particles, E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.range_cdsa_g_cm2, CSDA_range_g_cm2);

  // Conversion CSDA_proton => CSDA_ion
  long*  Z  =  (long*)calloc(number_of_particles, sizeof(long));
  long*  A  =  (long*)calloc(number_of_particles, sizeof(long));
  AT_Z_from_particle_no(        number_of_particles,
                                particle_no,
                                Z);
  AT_A_from_particle_no(        number_of_particles,
                                particle_no,
                                A);
  long i = 0;
  for (i = 0; i < number_of_particles; i++){
    if (particle_no[i] != 1){
      CSDA_range_g_cm2[i]  *=   (double)(A[i])/(double)((Z[i])*(Z[i]));
    }
  }

  free(Z);
  free(A);
}


double AT_CSDA_range_m_single(  const double  E_MeV_u,
		const long    particle_no,
		const long    material_no){

	double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);
	double CSDA_range_g_cm2 = AT_CSDA_range_g_cm2_single(E_MeV_u, particle_no, material_no);

	return CSDA_range_g_cm2 / (material_density_g_cm3 * 100.0);
}


void AT_CSDA_range_m(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    double        CSDA_range_m[])
{
  // Get material density
  double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);

  // Get mass-norm. CSDA range
  AT_CSDA_range_g_cm2(  number_of_particles,
      E_MeV_u,
      particle_no,
      material_no,
      CSDA_range_m);

  long  i;
  for (i = 0; i < number_of_particles; i++){
    CSDA_range_m[i]  /=  material_density_g_cm3 * 100.0;
  }

}


