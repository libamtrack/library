/**
 * @brief Material data
 */

/*
 *    AT_DataMaterial.c
 *    ==============
 *
 *    Created on: 09.01.2010
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

#include "AT_DataMaterial.h"

extern int AT_establish_LET_data( AT_single_material_data_struct* material);

///////////////////////////////////////////// SINGLE MATERIAL IMPLEMENTATION ///////////////////////////////////////////////////////


long AT_material_index_from_material_number( const long material_number ){
  long  index                =  -1;
  long  number_of_materials  =  1;
  find_elements_int(  &material_number,
      number_of_materials,
      AT_Material_Data.material_no,
      AT_Material_Data.n,
      &index);   // TODO replace call to find_elements_int by call to simpler function which will find the index
  return index;
}


void AT_material_name_from_number(
    const long material_no,
    char* material_name){

  long  index = AT_material_index_from_material_number( material_no );

  if( index != -1){
    strcpy(material_name, AT_Material_Data.material_name[index]);
  } else {
    strcpy(material_name,"*** invalid choice ***");
  }
}


long AT_material_number_from_name(
    const char* material_name){

  // find look-up index for material name in material data table
  long  match;
  const long n_tmp = 1;

  assert( material_name != NULL);

  find_elements_char(  &material_name,
      n_tmp,
      AT_Material_Data.material_name,
      AT_Material_Data.n,
      &match);

  if( match != -1){
    return AT_Material_Data.material_no[match];
  } else {
    return -1;
  }
}



double AT_density_g_cm3_from_material_no( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  return AT_Material_Data.density_g_cm3[index];
}



double AT_I_eV_from_material_no( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  return AT_Material_Data.I_eV[index];
}


double AT_alpha_g_cm2_MeV_from_material_no( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  return AT_Material_Data.alpha_g_cm2_MeV[index];
}

double AT_p_MeV_from_material_no( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  return AT_Material_Data.p_MeV[index];
}


double AT_m_g_cm2_from_material_no( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  return AT_Material_Data.m_g_cm2[index];
}

double AT_average_A_from_material_no( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  return AT_Material_Data.average_A[index];
}


double AT_average_Z_from_material_no( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  return AT_Material_Data.average_Z[index];
}

long AT_phase_from_material_no( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  return AT_Material_Data.phase[index];
}
void AT_get_material_data(     const long  material_no,
    double*  density_g_cm3,
    double*  I_eV,
    double*  alpha_g_cm2_MeV,
    double*  p_MeV,
    double*  m_g_cm2,
    double*  average_A,
    double*  average_Z){

  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    *density_g_cm3        =  0.0;
    *I_eV                 =  0.0;
    *alpha_g_cm2_MeV      =  0.0;
    *p_MeV                =  0.0;
    *m_g_cm2              =  0.0;
    *average_A            =  0.0;
    *average_Z            =  0.0;
    return ;
  }
  if( density_g_cm3 != NULL ){
      *density_g_cm3            =  AT_Material_Data.density_g_cm3[index];
  }
  if( I_eV != NULL ){
      *I_eV                     =  AT_Material_Data.I_eV[index];
  }
  if( alpha_g_cm2_MeV != NULL ){
      *alpha_g_cm2_MeV          =  AT_Material_Data.alpha_g_cm2_MeV[index];
  }
  if( p_MeV != NULL ){
      *p_MeV                    =  AT_Material_Data.p_MeV[index];
  }
  if( m_g_cm2 != NULL ){
      *m_g_cm2                  =  AT_Material_Data.m_g_cm2[index];
  }
  if( average_A != NULL ){
      *average_A                =  AT_Material_Data.average_A[index];
  }
  if( average_Z != NULL ){
      *average_Z                =  AT_Material_Data.average_Z[index];
  }
}


///////////////////////////////////////////// MULTIPLE MATERIAL IMPLEMENTATION ///////////////////////////////////////////////////////


void AT_get_materials_data( const long  number_of_materials,
    const long  material_no[],
    double  density_g_cm3[],
    double  I_eV[],
    double  alpha_g_cm2_MeV[],
    double  p_MeV[],
    double  m_g_cm2[],
    double  average_A[],
    double  average_Z[])
{
  long*  match  =  (long*)calloc(number_of_materials, sizeof(long));
  find_elements_int(  material_no,
      number_of_materials,
      AT_Material_Data.material_no,
      AT_Material_Data.n,
      match);

  long i;
  if( density_g_cm3 != NULL ){
    for(i = 0; i < number_of_materials; i++){
      density_g_cm3[i]      = AT_Material_Data.density_g_cm3[match[i]];
    }
  }
  if( I_eV != NULL ){
    for(i = 0; i < number_of_materials; i++){
      I_eV[i]      = AT_Material_Data.I_eV[match[i]];
    }
  }
  if( alpha_g_cm2_MeV != NULL ){
    for(i = 0; i < number_of_materials; i++){
      alpha_g_cm2_MeV[i]      = AT_Material_Data.alpha_g_cm2_MeV[match[i]];
    }
  }
  if( p_MeV != NULL ){
    for(i = 0; i < number_of_materials; i++){
      p_MeV[i]      = AT_Material_Data.p_MeV[match[i]];
    }
  }
  if( m_g_cm2 != NULL ){
    for(i = 0; i < number_of_materials; i++){
      m_g_cm2[i]      = AT_Material_Data.m_g_cm2[match[i]];
    }
  }
  if( average_A != NULL ){
    for(i = 0; i < number_of_materials; i++){
      average_A[i]      = AT_Material_Data.average_A[match[i]];
    }
  }
  if( average_Z != NULL ){
    for(i = 0; i < number_of_materials; i++){
      average_Z[i]      = AT_Material_Data.average_Z[match[i]];
    }
  }
  free(match);
}

double AT_plasma_energy_J_from_material_no(const long material_no){
	  long  index = AT_material_index_from_material_number( material_no );
	  if( index == -1){
#ifndef NDEBUG
	    printf("Material no %ld not found\n", material_no);
#endif
	    return 0.0;
	  }
	  double electron_density_m3 = AT_electron_density_m3_from_material_no_single(material_no);
	  return(AT_plasma_energy_J_single( electron_density_m3));
}

double AT_plasma_energy_J_single(const double electron_density_m3){
	  return(sqrt(4 * M_PI * electron_density_m3 * classical_electron_radius_m) * Dirac_constant_J_s * c_m_s);
}

double AT_electron_density_m3_from_material_no_single( const long   material_no )
{
  long  index = AT_material_index_from_material_number( material_no );
  if( index == -1){
#ifndef NDEBUG
    printf("Material no %ld not found\n", material_no);
#endif
    return 0.0;
  }
  double density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);
  double average_Z     = AT_average_Z_from_material_no(material_no);
  double average_A     = AT_average_A_from_material_no(material_no);
  return(AT_electron_density_m3_single( density_g_cm3,
    average_Z,
    average_A));
}

void AT_electron_density_m3_from_material_no_multi( const long n,
		const long   material_no[],
		double electron_density_m3[])
{
  long i;
  for (i = 0; i < n; i++){
	electron_density_m3[i]  = 	AT_electron_density_m3_from_material_no_single(material_no[i]);
  }
}


double AT_electron_density_m3_single( const double density_g_cm3,
    const double average_Z,
    const double average_A){

	// TODO: Again, number of nucleon rather than atomic mass is used
	double molar_mass_g_mol          = average_A;
	double number_of_atoms_per_g     = Avogadro_constant_1_mol / molar_mass_g_mol;
	double number_of_electrons_per_g = number_of_atoms_per_g * average_Z;
	double electron_density_per_cm3  = number_of_electrons_per_g * density_g_cm3;
	return(electron_density_per_cm3 * m_to_cm * m_to_cm * m_to_cm);
}

void AT_electron_density_m3_multi( const long n,
    const double density_g_cm3[],
    const double average_Z[],
    const double average_A[],
    double electron_density_m3[]){

	long i;
	for (i = 0; i < n; i++){
		electron_density_m3[i]   = AT_electron_density_m3_single(density_g_cm3[i],
				average_Z[i],
				average_A[i]);
	}
}

void AT_electron_density_m3_from_composition( const long n,
    const double density_g_cm3,
    const long Z[],
    const long A[],
    const double weight_fraction[],
    double* electron_density_m3)
{
	double* normalized_weight_fraction = (double*)calloc(n, sizeof(double));
	AT_normalize(n, weight_fraction, normalized_weight_fraction);

	*electron_density_m3 = 0.0;
	long i;
	for (i = 0; i < n; i++){
		*electron_density_m3 += normalized_weight_fraction[i] * Avogadro_constant_1_mol * Z[i] / A[i];
	}

	*electron_density_m3 *= density_g_cm3* m_to_cm * m_to_cm * m_to_cm;

	free(normalized_weight_fraction);
}


void AT_average_A_from_composition( const long n,
    const long A[],
    const double weight_fraction[],
    double* average_A)
{
	double* normalized_weight_fraction = (double*)calloc(n, sizeof(double));
	AT_normalize(n, weight_fraction, normalized_weight_fraction);

	*average_A = 0.0;
	long i;
	for (i = 0; i < n; i++){
		*average_A += normalized_weight_fraction[i] * A[i];
	}

	free(normalized_weight_fraction);
}

void AT_average_Z_from_composition( const long n,
    const long Z[],
    const double weight_fraction[],
    double* average_Z)
{
	double* normalized_weight_fraction = (double*)calloc(n, sizeof(double));
	AT_normalize(n, weight_fraction, normalized_weight_fraction);

	*average_Z = 0.0;
	long i;
	for (i = 0; i < n; i++){
		*average_Z += normalized_weight_fraction[i] * Z[i];
	}

	free(normalized_weight_fraction);
}

void AT_effective_Z_from_composition( const long n,
    const long Z[],
    const double weight_fraction[],
    const double electron_densities_cm3[],
    const double exponent,
    double* effective_Z)
{
	double* normalized_weight_fraction = (double*)calloc(n, sizeof(double));
	AT_normalize(n, weight_fraction, normalized_weight_fraction);

	double* normalized_electron_densities = (double*)calloc(n, sizeof(double));
    bool use_electron_densities = false;
    if(AT_sum(n, electron_densities_cm3) > 0){
    	AT_normalize(n, electron_densities_cm3, normalized_electron_densities);
    	use_electron_densities = true;
    }

	double effective_Z_sum = 0.0;
	long i;
	for (i = 0; i < n; i++){
		double tmp = normalized_weight_fraction[i] * pow(Z[i], exponent);
		if(use_electron_densities){
			tmp *= normalized_electron_densities[i];
		}
		effective_Z_sum += tmp;
	}
	*effective_Z = pow(effective_Z_sum, 1/exponent);

	free(normalized_electron_densities);
	free(normalized_weight_fraction);
}

void AT_I_eV_from_composition( const long n,
	const long A[],
    const long Z[],
    const double weight_fraction[],
    double* I_eV)
{
	double* normalized_weight_fraction = (double*)calloc(n, sizeof(double));
	AT_normalize(n, weight_fraction, normalized_weight_fraction);

	double* I_eV_elements   = (double*)calloc(n, sizeof(double));
	long*   particle_nos    = (long*)calloc(n, sizeof(long));
	AT_particle_no_from_Z_and_A(n, Z, A, particle_nos);
	AT_I_eV_from_particle_no(n, particle_nos, I_eV_elements);
	long i;
	double numerator = 0.0;
	double denominator = 0.0;

	for (i = 0; i < n; i++){
			numerator	   += (normalized_weight_fraction[i] * (Z[i])/A[i]) * log(I_eV_elements[i]);
	        denominator	   += normalized_weight_fraction[i] * (Z[i]/A[i]);
	        *I_eV 	= exp(numerator / denominator);
	}

	free(I_eV_elements);
	free(normalized_weight_fraction);
}

void AT_set_user_material( const double density_g_cm3,
		const double I_eV,
		const double average_A,
		const double average_Z,
		long* status)
{
	long  index = User_Defined_Material;
	AT_Material_Data.density_g_cm3[index]     = density_g_cm3;
	AT_Material_Data.I_eV[index]              = I_eV;
	AT_Material_Data.average_A[index]         = average_A;
	AT_Material_Data.average_Z[index]         = average_Z;
	AT_Material_Data.ready[index]             = true;
	*status = EXIT_SUCCESS;
}

void AT_set_user_material_from_composition( const long n,
		const double density_g_cm3,
		const long A[],
	    const long Z[],
	    const double weight_fraction[],
		long* status)
{

	double* normalized_weight_fraction = (double*)calloc(n, sizeof(double));
	AT_normalize(n, weight_fraction, normalized_weight_fraction);

	double I_eV			= 0.0;
	double average_Z	= 0.0;
	double average_A	= 0.0;
	AT_set_user_material(density_g_cm3,
			I_eV,
			average_A,
			average_Z,
			status);

	free(normalized_weight_fraction);
}
