/**
* @file
* @brief Wrapper functions
*
* C functions which are called from R cannot have input
* integer parameters of type "long". Only "int" type is
* accepted. Floating point parameters should have only
* "float" arguments.
* This file contains set of wrapper functions,
* which are casting int arguments to long if necessary
* and also from double to float.
*/

/*
 *    AT_Wrapper_R.c
 *    ==============
 *
 *    R Wrapper
 *    Created on: 20.02.2010
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

#include "AT_Wrapper_R.h"

void AT_D_Gy_R( const int* n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    float*        D_Gy)
{
  long i;
  const long n_long             = (long)(*n);
  const long material_no_long   = (long)(*material_no);
  long * particle_no_long       = (long*)calloc(n_long,sizeof(long));
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]                 = (long)particle_no[i];
  }

  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double   = (double*)calloc(n_long,sizeof(double));
  double * D_Gy_double   		= (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]                   = (double)E_MeV_u[i];
    fluence_cm2_double[i]               = (double)fluence_cm2[i];
  }

  AT_D_Gy( 	n_long,
			E_MeV_u_double,
			particle_no_long,
			fluence_cm2_double,
			material_no_long,
			D_Gy_double);

  for(i = 0 ; i < n_long ; i++){
    D_Gy[i]                   = (float)D_Gy_double[i];
  }		
			
			
  free(fluence_cm2_double);
  free(E_MeV_u_double);
  free(particle_no_long);
  free(D_Gy_double);
}


void AT_D_RDD_Gy_R( const int*  n,
    const float*  r_m,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    float*        D_RDD_Gy){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long er_model_long = (long)(*er_model);
  const long rdd_model_long = (long)(*rdd_model);
  const long particle_no_long = (long)(*particle_no);

  /* float -> double conversion */
  double * r_m_double = (double*)calloc(n_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_long ; i++){
    r_m_double[i] = (double)r_m[i];
  }
  double E_MeV_u_double = (double)(*E_MeV_u);
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameter[i];
  }

  /* place for results */
  double * D_RDD_Gy_double = (double*)calloc(n_long,sizeof(double));

  AT_D_RDD_Gy( n_long,
      r_m_double,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      D_RDD_Gy_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    D_RDD_Gy[i] = (float)D_RDD_Gy_double[i];
  }

  free(r_m_double);
  free(D_RDD_Gy_double);
}


void AT_RDD_f1_parameters_mixed_field_R( const int* n,
    const float* E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    const int*    rdd_model,
    const float*    rdd_parameters,
    const int*    er_model,
    float*    f1_parameters){

  long i;

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long er_model_long = (long)(*er_model);
  const long rdd_model_long = (long)(*rdd_model);
  long * particle_no_long       = (long*)calloc(n_long,sizeof(long));
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]                 = (long)particle_no[i];
  }

  /* float -> double conversion */
  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]                   = (double)E_MeV_u[i];
  }

  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
	  rdd_parameter_double[i] = (double)rdd_parameters[i];
  }

  /* place for results */
  double * f1_parameters_double = (double*)calloc(n_long*AT_SC_F1_PARAMETERS_SINGLE_LENGTH,sizeof(double));

  AT_RDD_f1_parameters_mixed_field( n_long,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      f1_parameters_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long*AT_SC_F1_PARAMETERS_SINGLE_LENGTH; i++){
	  f1_parameters[i] = (float)f1_parameters_double[i];
  }

  free(E_MeV_u_double);
  free(particle_no_long);
  free(f1_parameters_double);
}


void AT_fluence_cm2_R( const int* n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*        D_Gy,
    const int*    material_no,
    float*  fluence_cm2)
{
  long i;
  const long n_long             = (long)(*n);
  const long material_no_long   = (long)(*material_no);
  long * particle_no_long       = (long*)calloc(n_long,sizeof(long));
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]                 = (long)particle_no[i];
  }

  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double   = (double*)calloc(n_long,sizeof(double));
  double * D_Gy_double   		= (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]               = (double)E_MeV_u[i];
    D_Gy_double[i]               	= (double)D_Gy[i];
  }

  AT_fluence_cm2( 	n_long,
					E_MeV_u_double,
					particle_no_long,
					D_Gy_double,
					material_no_long,
					fluence_cm2_double);

  for(i = 0 ; i < n_long ; i++){
	  fluence_cm2[i]                   = (float)fluence_cm2_double[i];
  }


  free(fluence_cm2_double);
  free(E_MeV_u_double);
  free(particle_no_long);
  free(D_Gy_double);
}


void AT_gamma_response_R( const int*  n,
    const float*  d_Gy,
    const int*    gamma_model,
    const float*  gamma_parameter,
    const int*    lethal_events_mode,
	float*        S){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long gamma_model_long = (long)(*gamma_model);

  /* int -> bool conversion */
  const bool lethal_events_mode_bool = (bool)(*lethal_events_mode);

  /* float -> double conversion */
  double * d_Gy_double = (double*)calloc(n_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_long ; i++){
    d_Gy_double[i] = (double)d_Gy[i];
  }
  long number_of_gamma_parameters = 0;
  if ( gamma_model_long == GR_GeneralTarget ){
    while  (gamma_parameter[number_of_gamma_parameters] != 0){
      number_of_gamma_parameters  += 4;
    }
    number_of_gamma_parameters++; /* to include also trailing zero */
  } else {
    number_of_gamma_parameters = AT_Gamma_number_of_parameters( gamma_model_long );
  }
  double * gamma_parameter_double  = (double*)calloc(number_of_gamma_parameters,sizeof(double));
  for(i = 0 ; i < number_of_gamma_parameters ; i++){
    gamma_parameter_double[i] = (double)gamma_parameter[i];
  }

  /* place for results */
  double * S_double = (double*)calloc(n_long,sizeof(double));

  AT_gamma_response(n_long,
      d_Gy_double,
      gamma_model_long,
      gamma_parameter_double,
      lethal_events_mode_bool,
      S_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    S[i] = (float)S_double[i];
  }

  free(d_Gy_double);
  free(S_double);
}


void AT_LET_MeV_cm2_g_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    float*        LET_MeV_cm2_g){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  long * particle_no_long = (long*)calloc(n_long,sizeof(long));
  long i;
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
  }

  /* place for results */
  double * LET_MeV_cm2_g_double = (double*)calloc(n_long,sizeof(double));

  AT_LET_MeV_cm2_g(  n_long,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      LET_MeV_cm2_g_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    LET_MeV_cm2_g[i] = (float)LET_MeV_cm2_g_double[i];
  }

  free(E_MeV_u_double);
  free(LET_MeV_cm2_g_double);
}

void AT_LET_keV_um_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    float*        LET_keV_um){

	  /* int -> long conversion */
	  const long n_long = (long)(*n);
	  const long material_no_long = (long)(*material_no);
	  long * particle_no_long = (long*)calloc(n_long,sizeof(long));
	  long i;
	  for(i = 0 ; i < n_long ; i++){
	    particle_no_long[i] = (long)particle_no[i];
	  }

	  /* float -> double conversion */
	  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
	  for(i = 0 ; i < n_long ; i++){
	    E_MeV_u_double[i] = (double)E_MeV_u[i];
	  }

	  /* place for results */
	  double * LET_keV_um_double = (double*)calloc(n_long,sizeof(double));

	  AT_LET_keV_um(  n_long,
	      E_MeV_u_double,
	      particle_no_long,
	      material_no_long,
	      LET_keV_um_double);

	  /* double -> float conversion (results) */
	  for(i = 0 ; i < n_long ; i++){
	    LET_keV_um[i] = (float)LET_keV_um_double[i];
	  }

	  free(E_MeV_u_double);
	  free(LET_keV_um_double);
}



void AT_max_E_transfer_MeV_R(  const int*  n,
    const float*  E_MeV_u,
    float*        max_E_transfer_MeV){

  /* int -> long conversion */
  long  n_long           = (long)*n;

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
  }

  /* place for results */
  double * max_E_transfer_MeV_double = (double*)calloc(n_long,sizeof(double));

  AT_max_E_transfer_MeV(  n_long,
      E_MeV_u_double,
      max_E_transfer_MeV_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    max_E_transfer_MeV[i] = (float)max_E_transfer_MeV_double[i];
  }

  free( max_E_transfer_MeV_double );
  free( E_MeV_u_double );
}


void AT_max_electron_ranges_m_R(  const int*  number_of_particles,
    const float*  E_MeV_u,
    const int*    material_no,
    const int*    er_model,
    float*        max_electron_range_m)
{
  /* int -> long conversion */
  const long number_of_particles_long = (long)(*number_of_particles);

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(*number_of_particles,sizeof(double));
  long i;
  for(i = 0 ; i < number_of_particles_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
  }

  /* place for results */
  double * max_electron_range_m_double = (double*)calloc(*number_of_particles,sizeof(double));

  AT_max_electron_ranges_m( number_of_particles_long,
      E_MeV_u_double,
      *material_no,
      *er_model,
      max_electron_range_m_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < number_of_particles_long ; i++){
    max_electron_range_m[i] = (float)max_electron_range_m_double[i];
  }

  free( max_electron_range_m_double );
  free( E_MeV_u_double );
}

void AT_particle_no_from_Z_and_A_R  ( const int*  n,
    const int*    Z,
    const int*    A,
    int*          particle_no)
{
  long  n_long = (long)*n;
  long* Z_long = (long*)calloc(n_long,sizeof(long));
  long* A_long = (long*)calloc(n_long,sizeof(long));
  long* particle_no_long = (long*)calloc(n_long,sizeof(long));

  long i;
  for(i = 0 ; i < n_long ; i++){
    Z_long[i] = (long)Z[i];
    A_long[i] = (long)A[i];
  }

  AT_particle_no_from_Z_and_A( n_long,
      Z_long,
      A_long,
      particle_no_long);

  for(i = 0 ; i < n_long ; i++){
    particle_no[i] = (int)particle_no_long[i];
  }

  free(Z_long);
  free(A_long);
  free(particle_no_long);
}

void AT_r_RDD_m_R  ( const int*  n,
    const float*  D_RDD_Gy,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    float*        r_RDD_m){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long er_model_long = (long)(*er_model);
  const long rdd_model_long = (long)(*rdd_model);
  const long particle_no_long = (long)(*particle_no);

  /* float -> double conversion */
  double * D_RDD_Gy_double = (double*)calloc(n_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_long ; i++){
    D_RDD_Gy_double[i] = (double)D_RDD_Gy[i];
  }
  double E_MeV_u_double = (double)(*E_MeV_u);
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameter[i];
  }

  /* place for results */
  double * r_RDD_m_double = (double*)calloc(n_long,sizeof(double));

  AT_r_RDD_m( n_long,
      D_RDD_Gy_double,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      r_RDD_m_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    r_RDD_m[i] = (float)r_RDD_m_double[i];
  }

  free( r_RDD_m_double );
  free( D_RDD_Gy_double );
}


void AT_KatzModel_inactivation_probability_R( const int* n,
    const float*   r_m,
    const float*   E_MeV_u,
    const int*     particle_no,
    const int*     material_no,
    const int*     rdd_model,
    const float*   rdd_parameters,
    const int*     er_model,
    const float*   gamma_parameters,
    float*         inactivation_probability){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long material_no_long = (long)(*material_no);
  const long particle_no_long = (long)(*particle_no);

  /* float -> double conversion */
  double * r_m_double = (double*)calloc(n_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_long ; i++){
    r_m_double[i] = (double)r_m[i];
  }
  double E_MeV_u_double = (double)(*E_MeV_u);
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameters[i];
  }
  double gamma_parameter_double[5];
  for(i = 0 ; i < 5 ; i++){
    gamma_parameter_double[i] = (double)gamma_parameters[i];
  }

  /* place for results */
  double * inactivation_probability_double = (double*)calloc(n_long,sizeof(double));

  AT_KatzModel_inactivation_probability( n_long,
      r_m_double,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      gamma_parameter_double,
      inactivation_probability_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    inactivation_probability[i] = (float)inactivation_probability_double[i];
  }

  free( r_m_double );
}


void AT_KatzModel_inactivation_cross_section_m2_R( const int* n,
    const float*   E_MeV_u,
    const int*     particle_no,
    const int*     material_no,
    const int*     rdd_model,
    const float*   rdd_parameters,
    const int*     er_model,
    const float*   gamma_parameters,
    float*         inactivation_cross_section_m2){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long material_no_long = (long)(*material_no);
  const long particle_no_long = (long)(*particle_no);

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
  }
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameters[i];
  }
  double gamma_parameter_double[5];
  for(i = 0 ; i < 5 ; i++){
    gamma_parameter_double[i] = (double)gamma_parameters[i];
  }

  /* place for results */
  double * inactivation_cross_section_m2_double = (double*)calloc(n_long,sizeof(double));

  AT_KatzModel_inactivation_cross_section_m2( n_long,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      gamma_parameter_double,
      inactivation_cross_section_m2_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    inactivation_cross_section_m2[i] = (float)inactivation_cross_section_m2_double[i];
  }

  free( E_MeV_u_double );
}

void AT_particle_name_from_particle_no_R(const int* particle_no,
    char** particle_name){

  const long n = 1;
  const long particle_no_long = (long)(*particle_no);
  char particle_name_str[1][6];

  strcpy(particle_name_str[0], *particle_name);

  AT_particle_name_from_particle_no( n,
      &particle_no_long,
      particle_name_str);

  strcpy(*particle_name, particle_name_str[0]);
}


void AT_particle_no_from_particle_name_R(const char** particle_name,
    int* particle_no){
  const long n = 1;
  char * particle_name_str = (char*)calloc(PARTICLE_NAME_NCHAR, sizeof(char));
  strcpy(particle_name_str, *particle_name);

  long particle_no_long = 0;

  AT_particle_no_from_particle_name( n,
      &particle_name_str,
      &particle_no_long);

  free( particle_name_str);

  *particle_no = (int)particle_no_long;
}

void AT_material_name_from_number_R( const int* material_no,
    char** material_name){

	const long material_no_long = (long)(*material_no);
	char * material_name_str = (char*)calloc(MATERIAL_NAME_LENGTH, sizeof(char));

	AT_material_name_from_number( material_no_long,
	      material_name_str);

	strcpy(*material_name, material_name_str);

	free(material_name_str);
}

void AT_material_number_from_name_R( const char** material_name,
	int* material_no){

	  char * material_name_str = (char*)calloc(MATERIAL_NAME_LENGTH, sizeof(char));
	  strcpy(material_name_str, *material_name);

	  long material_no_long = 0;

	  material_no_long	= AT_material_number_from_name( material_name_str);

	  free( material_name_str);

	  *material_no = (int)material_no_long;
}

void AT_run_GSM_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2_or_dose_Gy,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameters,
    const int*    er_model,
    const int*    gamma_model,
    const float*  gamma_parameters,
    const int*    N_runs,
    const int*    write_output,
    const int*    nX,
    const float*  voxel_size_m,
    const int*    lethal_events_mode,
    float*        results){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long gamma_model_long = (long)(*gamma_model);
  const long material_no_long = (long)(*material_no);
  long * particle_no_long = (long*)calloc(n_long,sizeof(long));
  long i;
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }
  const long N_runs_long = (long)(*N_runs);
  const long nX_long = (long)(*nX);

  /* int -> bool conversion */
  const bool write_output_bool = (bool)(*write_output);
  const bool lethal_events_mode_bool = (bool)(*lethal_events_mode);

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_or_dose_Gy_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
    fluence_cm2_or_dose_Gy_double[i] = (double)fluence_cm2_or_dose_Gy[i];
  }
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameters[i];
  }
  long number_of_gamma_parameters = 0;
  if ( gamma_model_long == GR_GeneralTarget ){
    while  (gamma_parameters[number_of_gamma_parameters] != 0){
      number_of_gamma_parameters  += 4;
    }
    number_of_gamma_parameters++; /* to include also trailing zero */
  } else {
    number_of_gamma_parameters = AT_Gamma_number_of_parameters( gamma_model_long );
  }
  double * gamma_parameter_double  = (double*)calloc(number_of_gamma_parameters,sizeof(double));
  for(i = 0 ; i < number_of_gamma_parameters ; i++){
    gamma_parameter_double[i] = (double)gamma_parameters[i];
  }
  double voxel_size_m_double = (double)(*voxel_size_m);

  /* place for results */
  double results_double[10];

  AT_run_GSM_method(n_long,
      E_MeV_u_double,
      particle_no_long,
      fluence_cm2_or_dose_Gy_double,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      gamma_model_long,
      gamma_parameter_double,
      N_runs_long,
      write_output_bool,
      nX_long,
      voxel_size_m_double,
      lethal_events_mode_bool,
      results_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < 10 ; i++){
    results[i] = (float)results_double[i];
  }
  free(particle_no_long);
  free(E_MeV_u_double);
  free(fluence_cm2_or_dose_Gy_double);
}


void AT_run_SPIFF_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2_or_dose_Gy,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameters,
    const int*    er_model,
    const int*    gamma_model,
    const float*  gamma_parameters,
    int*          N2,
    const float*  fluence_factor,
    const int*    write_output,
    const int*    shrink_tails,
    const float*  shrink_tails_under,
    const int*    adjust_N2,
    const int*    lethal_events_mode,
    float*       results){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long gamma_model_long = (long)(*gamma_model);
  const long material_no_long = (long)(*material_no);
  long * particle_no_long = (long*)calloc(n_long,sizeof(long));
  long i;
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }
  long N2_long = (long)(*N2);

  /* int -> bool conversion */
  const bool write_output_bool = (bool)(*write_output);
  const bool shrink_tails_bool = (bool)(*shrink_tails);
  const bool adjust_N2_bool = (bool)(*adjust_N2);
  const bool lethal_events_mode_bool = (bool)(*lethal_events_mode);

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_or_dose_Gy_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
    fluence_cm2_or_dose_Gy_double[i] = (double)fluence_cm2_or_dose_Gy[i];
  }
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameters[i];
  }
  long number_of_gamma_parameters = 0;
  if ( gamma_model_long == GR_GeneralTarget ){
    while  (gamma_parameters[number_of_gamma_parameters] != 0){
      number_of_gamma_parameters  += 4;
    }
    number_of_gamma_parameters++; /* to include also trailing zero */
  } else {
    number_of_gamma_parameters = AT_Gamma_number_of_parameters( gamma_model_long );
  }
  double * gamma_parameter_double  = (double*)calloc(number_of_gamma_parameters,sizeof(double));
  for(i = 0 ; i < number_of_gamma_parameters ; i++){
    gamma_parameter_double[i] = (double)gamma_parameters[i];
  }
  double fluence_factor_double = (double)(*fluence_factor);
  double shrink_tails_under_double = (double)(*shrink_tails_under);

  /* place for results */
  double results_double[10];

  AT_run_SPIFF_method(  n_long,
      E_MeV_u_double,
      particle_no_long,
      fluence_cm2_or_dose_Gy_double,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      gamma_model_long,
      gamma_parameter_double,
      N2_long,
      fluence_factor_double,
      write_output_bool,
      shrink_tails_bool,
      shrink_tails_under_double,
      adjust_N2_bool,
      lethal_events_mode_bool,
      results_double);

  /* long -> int conversion (results) */
   *N2 = (int)N2_long;

   /* double -> float conversion (results) */
   for(i = 0 ; i < 10 ; i++){
     results[i] = (float)results_double[i];
   }
   free(particle_no_long);
   free(E_MeV_u_double);
   free(fluence_cm2_or_dose_Gy_double);
}


void AT_run_IGK_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2_or_dose_Gy,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameters,
    const int*    er_model,
    const int*    gamma_model,
    const float*  gamma_parameters,
    const float*  saturation_cross_section_factor,
	const int*    write_output,
    float*       results){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long gamma_model_long = (long)(*gamma_model);
  const long material_no_long = (long)(*material_no);
  long * particle_no_long = (long*)calloc(n_long,sizeof(long));
  long i;
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  /* int -> bool conversion */
  const bool write_output_bool = (bool)(*write_output);

  /* float -> double conversion */
  const double saturation_cross_section_factor_double = (const double)(*saturation_cross_section_factor);
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_or_dose_Gy_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
    fluence_cm2_or_dose_Gy_double[i] = (double)fluence_cm2_or_dose_Gy[i];
  }
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameters[i];
  }
  long number_of_gamma_parameters = 0;
  if ( gamma_model_long == GR_GeneralTarget ){
    while  (gamma_parameters[number_of_gamma_parameters] != 0){
      number_of_gamma_parameters  += 4;
    }
    number_of_gamma_parameters++; /* to include also trailing zero */
  } else {
    number_of_gamma_parameters = AT_Gamma_number_of_parameters( gamma_model_long );
  }
  double * gamma_parameter_double  = (double*)calloc(number_of_gamma_parameters,sizeof(double));
  for(i = 0 ; i < number_of_gamma_parameters ; i++){
    gamma_parameter_double[i] = (double)gamma_parameters[i];
  }

  /* place for results */
  double results_double[10];

  AT_run_IGK_method(  n_long,
      E_MeV_u_double,
      particle_no_long,
      fluence_cm2_or_dose_Gy_double,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      gamma_model_long,
      gamma_parameter_double,
      saturation_cross_section_factor_double,
      write_output_bool,
      results_double);

    /* double -> float conversion (results) */
   for(i = 0 ; i < 10 ; i++){
     results[i] = (float)results_double[i];
   }
   free(particle_no_long);
   free(E_MeV_u_double);
   free(fluence_cm2_or_dose_Gy_double);
}


void  AT_SC_get_f1_array_size_R(
    const int*    n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    const int*    N2,
    int*          n_bins_f1,
    float*        f1_parameters){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long material_no_long = (long)(*material_no);
  const long N2_long = (long)(*N2);
  long * particle_no_long = (long*)calloc(n_long,sizeof(long));
  long i;
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }
  long n_bins_f1_long;

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
  }
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameter[i];
  }

  /* place for results */
  double * f1_parameters_double = (double*)calloc( n_long * AT_SC_F1_PARAMETERS_SINGLE_LENGTH, sizeof(double));

  n_bins_f1_long = AT_SC_get_f1_array_size( n_long,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      N2_long);

  AT_RDD_f1_parameters_mixed_field( n_long,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      f1_parameters_double);

  /* long -> int conversion (results) */
  *n_bins_f1 = (int)n_bins_f1_long;

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long * AT_SC_F1_PARAMETERS_SINGLE_LENGTH ; i++){
    f1_parameters[i] = (float)f1_parameters_double[i];
  }
  free(f1_parameters_double);
  free(particle_no_long);
  free(E_MeV_u_double);
}


void  AT_SC_get_f1_R(
    const int*    n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2_or_dose_Gy,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    const int*    N2,
    const int*    n_bins_f1,
    const float*  f1_parameters,
    float*        f1_d_Gy,
    float*        f1_dd_Gy,
    float*        f1){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long material_no_long = (long)(*material_no);
  const long N2_long = (long)(*N2);
  const long n_bins_f1_long = (long)(*n_bins_f1);
  long * particle_no_long = (long*)calloc(n_long,sizeof(long));
  long i;
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_or_dose_Gy_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
    fluence_cm2_or_dose_Gy_double[i] = (double)fluence_cm2_or_dose_Gy[i];
  }
  double rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
    rdd_parameter_double[i] = (double)rdd_parameter[i];
  }
  double * f1_parameters_double = (double*)calloc(n_long * AT_SC_F1_PARAMETERS_SINGLE_LENGTH,sizeof(double));
  for(i = 0 ; i < n_long * AT_SC_F1_PARAMETERS_SINGLE_LENGTH ; i++){
    f1_parameters_double[i] = (double)f1_parameters[i];
  }

  /* place for results */
  double * f1_d_Gy_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  double * f1_dd_Gy_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  double * f1_double = (double*)calloc(n_bins_f1_long,sizeof(double));

  AT_SC_get_f1( n_long,
      E_MeV_u_double,
      particle_no_long,
      fluence_cm2_or_dose_Gy_double,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      N2_long,
      n_bins_f1_long,
      f1_parameters_double,
      f1_d_Gy_double,
      f1_dd_Gy_double,
      f1_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_bins_f1_long ; i++){
    f1_d_Gy[i] = (float)f1_d_Gy_double[i];
    f1_dd_Gy[i] = (float)f1_dd_Gy_double[i];
    f1[i] = (float)f1_double[i];
  }
  free(particle_no_long);
  free(E_MeV_u_double);
  free(fluence_cm2_or_dose_Gy_double);
  free(f1_parameters_double);
  free(f1_d_Gy_double);
  free(f1_dd_Gy_double);
  free(f1_double);
}


void  AT_SC_get_f_array_size_R(
    const float*  u,
    const float*  fluence_factor,
    const int*    N2,
    const int*    n_bins_f1,
    const float*  f1_d_Gy,
    const float*  f1_dd_Gy,
    const float*  f1,
    int*          n_bins_f,
    float*        u_start,
    int*          n_convolutions){

  /* int -> long conversion */
  const long N2_long = (long)(*N2);
  const long n_bins_f1_long = (long)(*n_bins_f1);

  /* float -> double conversion */
  double u_double = (double)(*u);
  double fluence_factor_double = (double)(*fluence_factor);
  double * f1_d_Gy_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  double * f1_dd_Gy_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  double * f1_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_bins_f1_long ; i++){
    f1_d_Gy_double[i]   =  (double)f1_d_Gy[i];
    f1_dd_Gy_double[i]  =  (double)f1_dd_Gy[i];
    f1_double[i]        =  (double)f1[i];
  }

  /* place for results */
  long n_bins_f_long;
  long n_convolutions_long;
  double u_start_double;

  AT_SC_get_f_array_size(  u_double,
      fluence_factor_double,
      N2_long,
      n_bins_f1_long,
      f1_d_Gy_double,
      f1_dd_Gy_double,
      f1_double,
      &n_bins_f_long,
      &u_start_double,
      &n_convolutions_long);

  /* long -> int conversion (results) */
  *n_bins_f = (int)n_bins_f_long;
  *n_convolutions = (int)n_convolutions_long;

  /* double -> float conversion (results) */
  *u_start = (float)u_start_double;

  free(f1_d_Gy_double);
  free(f1_dd_Gy_double);
  free(f1_double);
}


void  AT_SC_get_f_start_R( const int*    n_bins_f1,
    const int*    N2,
    const float*  f1_d_Gy,
    const float*  f1_dd_Gy,
    const float*  f1,
    const int*    n_bins_f,
    float*        f_d_Gy,
    float*        f_dd_Gy,
    float*        f_start){

  /* int -> long conversion */
  const long n_bins_f1_long = (long)(*n_bins_f1);
  const long N2_long = (long)(*N2);
  const long n_bins_f_long = (long)(*n_bins_f);

  /* float -> double conversion */
  double * f1_d_Gy_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  double * f1_dd_Gy_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  double * f1_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_bins_f1_long ; i++){
    f1_d_Gy_double[i]   =  (double)f1_d_Gy[i];
    f1_dd_Gy_double[i]  =  (double)f1_dd_Gy[i];
    f1_double[i]        =  (double)f1[i];
  }

  /* place for results */
  double * f_d_Gy_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double * f_dd_Gy_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double * f_start_double = (double*)calloc(n_bins_f_long,sizeof(double));

  AT_SC_get_f_start( n_bins_f1_long,
      N2_long,
      f1_d_Gy_double,
      f1_dd_Gy_double,
      f1_double,
      n_bins_f_long,
      f_d_Gy_double,
      f_dd_Gy_double,
      f_start_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_bins_f_long ; i++){
    f_d_Gy[i] = (float)f_d_Gy_double[i];
    f_dd_Gy[i] = (float)f_dd_Gy_double[i];
    f_start[i] = (float)f_start_double[i];
  }
  free(f1_d_Gy_double);
  free(f1_dd_Gy_double);
  free(f1_double);
  free(f_d_Gy_double);
  free(f_dd_Gy_double);
  free(f_start_double);
}


void AT_SC_get_gamma_response_R(  const int* number_of_bins,
    const float*   d_Gy,
    const float*   dd_Gy,
    const float*   f,
    const float*   f0,
    const int*     gamma_model,
    const float*   gamma_parameter,
    const int*     lethal_events_mode,
    float*         S,
    float*         S_HCP,
    float*         S_gamma,
    float*         efficiency){

  const long number_of_bins_long        = (long)(*number_of_bins);
  const double f0_double                = (const double)*f0;
  const long gamma_model_long           = (long)*gamma_model;
  const bool lethal_events_mode_bool    = (const bool)*lethal_events_mode;

  double * d_Gy_double                  = (double*)calloc(number_of_bins_long,sizeof(double));
  double * dd_Gy_double                 = (double*)calloc(number_of_bins_long,sizeof(double));
  double * f_double                     = (double*)calloc(number_of_bins_long,sizeof(double));
  long i;
  for(i = 0 ; i < number_of_bins_long ; i++){
    d_Gy_double[i]   =  (double)d_Gy[i];
    dd_Gy_double[i]  =  (double)dd_Gy[i];
    f_double[i]      =  (double)f[i];
  }

  const long n_gamma_parameters         = 5; //AT_Gamma_number_of_parameters( gamma_model_long);
  double * gamma_parameters_double      = (double*)calloc(n_gamma_parameters,sizeof(double));
  for(i = 0 ; i < n_gamma_parameters ; i++){
    gamma_parameters_double[i]   =  (double)gamma_parameter[i];
  }

  /* place for results */
  double * S_double                     = (double*)calloc(number_of_bins_long,sizeof(double));
  double S_HCP_double, S_gamma_double, efficiency_double;

  AT_get_gamma_response(  number_of_bins_long,
    d_Gy_double,
    dd_Gy_double,
    f_double,
    f0_double,
    gamma_model_long,
    gamma_parameters_double,
    lethal_events_mode_bool,
    // return
    S_double,
    &S_HCP_double,
    &S_gamma_double,
    &efficiency_double);

  for(i = 0 ; i < number_of_bins_long ; i++){
    S[i]   =  (float)S_double[i];
  }

  *S_HCP      = (float)S_HCP_double;
  *S_gamma    = (float)S_gamma_double;
  *efficiency = (float)efficiency_double;

  free(S_double);
  free(gamma_parameters_double);
  free(f_double);
  free(dd_Gy_double);
  free(d_Gy_double);
}


void AT_SuccessiveConvolutions_R( const float*  u,
    const int*    n_bins_f,
    int*          N2,
    int*          n_bins_f_used,
    float*        f_d_Gy,
    float*        f_dd_Gy,
    float*        f,
    float*        f0,
    float*        fdd,
    float*        dfdd,
    float*        d,
    const int*    write_output,
    const int*    shrink_tails,
    const float*  shrink_tails_under,
    const int*    adjust_N2){

  /* int -> long conversion */
  const long n_bins_f_long = (long)(*n_bins_f);
  long N2_long = (long)(*N2);
  long n_bins_f_used_long = (long)(*n_bins_f_used);

  /* int -> bool conversion */
  const bool write_output_bool = (bool)(*write_output);
  const bool shrink_tails_bool = (bool)(*shrink_tails);
  const bool adjust_N2_bool = (bool)(*adjust_N2);

  /* float -> double conversion */
  double u_double = (double)(*u);
  double shrink_tails_under_double = (double)(*shrink_tails_under);

  /* place for results */
  double * f_d_Gy_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double * f_dd_Gy_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double * f_double = (double*)calloc(n_bins_f_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_bins_f_long ; i++){
    f_d_Gy_double[i]   =  (double)f_d_Gy[i];
    f_dd_Gy_double[i]  =  (double)f_dd_Gy[i];
    f_double[i]        =  (double)f[i];
  }
  double f0_double;
  double * fdd_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double * dfdd_double = (double*)calloc(n_bins_f_long,sizeof(double));
  double d_double;

  AT_SuccessiveConvolutions( u_double,
      n_bins_f_long,
      &N2_long,
      &n_bins_f_used_long,
      f_d_Gy_double,
      f_dd_Gy_double,
      f_double,
      &f0_double,
      fdd_double,
      dfdd_double,
      &d_double,
      write_output_bool,
      shrink_tails_bool,
      shrink_tails_under_double,
      adjust_N2_bool);

  /* long -> int conversion (results) */
  *N2   = (int)N2_long;
  *n_bins_f_used = (int)n_bins_f_used_long;

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_bins_f_long ; i++){
    f_d_Gy[i] = (float)f_d_Gy_double[i];
    f_dd_Gy[i] = (float)f_dd_Gy_double[i];
    f[i] = (float)f_double[i];
    fdd[i] = (float)fdd_double[i];
    dfdd[i] = (float)dfdd_double[i];
  }
  *f0 = (float)f0_double;
  *d = (float)d_double;
}


void AT_total_D_Gy_R( const int* n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    float*        total_D_Gy)
{
  long i;
  const long n_long             = (long)(*n);
  const long material_no_long   = (long)(*material_no);
  long * particle_no_long       = (long*)calloc(n_long,sizeof(long));
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]                 = (long)particle_no[i];
  }

  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double   = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]                   = (double)E_MeV_u[i];
    fluence_cm2_double[i]               = (double)fluence_cm2[i];
  }

  *total_D_Gy   = (float)AT_total_D_Gy( n_long,
                                        E_MeV_u_double,
                                        particle_no_long,
                                        fluence_cm2_double,
                                        material_no_long);

  free(fluence_cm2_double);
  free(E_MeV_u_double);
  free(particle_no_long);
}


void AT_total_fluence_cm2_R( const int* n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    float*        total_fluence_cm2)
{
  long i;
  const long n_long             = (long)(*n);
  const long material_no_long   = (long)(*material_no);
  long * particle_no_long       = (long*)calloc(n_long,sizeof(long));
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]                 = (long)particle_no[i];
  }

  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double   = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]                   = (double)E_MeV_u[i];
    fluence_cm2_double[i]               = (double)fluence_cm2[i];
  }

  *total_fluence_cm2 = (float)AT_total_fluence_cm2( n_long,
                                        E_MeV_u_double,
                                        particle_no_long,
                                        fluence_cm2_double,
                                        material_no_long);

  free(fluence_cm2_double);
  free(E_MeV_u_double);
  free(particle_no_long);
}


void AT_fluence_weighted_E_MeV_u_R( const int* n,
    const float*  E_MeV_u,
    const float*  fluence_cm2,
    float* fluence_weighted_E_MeV_u)
{
  long i;
  const long n_long             = (long)(*n);
  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double   = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]                   = (double)E_MeV_u[i];
    fluence_cm2_double[i]               = (double)fluence_cm2[i];
  }

  *fluence_weighted_E_MeV_u = (float)AT_fluence_weighted_E_MeV_u( n_long,
                                        E_MeV_u_double,
                                        fluence_cm2_double);

  free(fluence_cm2_double);
  free(E_MeV_u_double);
}


void AT_dose_weighted_E_MeV_u_R( const int* n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    float*        dose_weighted_E_MeV_u)
{
  long i;
  const long n_long             = (long)(*n);
  const long material_no_long   = (long)(*material_no);
  long * particle_no_long       = (long*)calloc(n_long,sizeof(long));
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]                 = (long)particle_no[i];
  }

  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double   = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]                   = (double)E_MeV_u[i];
    fluence_cm2_double[i]               = (double)fluence_cm2[i];
  }

  *dose_weighted_E_MeV_u = (float)AT_dose_weighted_E_MeV_u( n_long,
                                        E_MeV_u_double,
                                        particle_no_long,
                                        fluence_cm2_double,
                                        material_no_long);

  free(fluence_cm2_double);
  free(E_MeV_u_double);
  free(particle_no_long);
}


void AT_fluence_weighted_LET_MeV_cm2_g_R( const int* n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    float*        fluence_weighted_LET_MeV_cm2_g)
{
  long i;
  const long n_long             = (long)(*n);
  const long material_no_long   = (long)(*material_no);
  long * particle_no_long       = (long*)calloc(n_long,sizeof(long));
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]                 = (long)particle_no[i];
  }

  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double   = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]                   = (double)E_MeV_u[i];
    fluence_cm2_double[i]               = (double)fluence_cm2[i];
  }

  *fluence_weighted_LET_MeV_cm2_g = (float)AT_fluence_weighted_LET_MeV_cm2_g( n_long,
                                        E_MeV_u_double,
                                        particle_no_long,
                                        fluence_cm2_double,
                                        material_no_long);

  free(fluence_cm2_double);
  free(E_MeV_u_double);
  free(particle_no_long);
}


void AT_dose_weighted_LET_MeV_cm2_g_R( const int* n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    float*        dose_weighted_LET_MeV_cm2_g)
{
  long i;
  const long n_long             = (long)(*n);
  const long material_no_long   = (long)(*material_no);
  long * particle_no_long       = (long*)calloc(n_long,sizeof(long));
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]                 = (long)particle_no[i];
  }

  double * E_MeV_u_double       = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double   = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]                   = (double)E_MeV_u[i];
    fluence_cm2_double[i]               = (double)fluence_cm2[i];
  }

  *dose_weighted_LET_MeV_cm2_g = (float)AT_dose_weighted_LET_MeV_cm2_g( n_long,
                                        E_MeV_u_double,
                                        particle_no_long,
                                        fluence_cm2_double,
                                        material_no_long);

  free(fluence_cm2_double);
  free(E_MeV_u_double);
  free(particle_no_long);
}


void AT_get_materials_data_R( const int*  number_of_materials,
    const int*  material_no,
    float*  density_g_cm3,
    float*  electron_density_m3,
    float*  I_eV,
    float*  alpha_g_cm2_MeV,
    float*  p_MeV,
    float*  m_g_cm2,
    float*  average_A,
    float*  average_Z)
{
  const long number_of_materials_long           = (const long)*number_of_materials;
  long* material_no_long                        = (long*)calloc(number_of_materials_long,sizeof(long));
  double* density_g_cm3_double                  = (double*)calloc(number_of_materials_long,sizeof(double));
  double* electron_density_m3_double            = (double*)calloc(number_of_materials_long,sizeof(double));
  double* I_eV_double                           = (double*)calloc(number_of_materials_long,sizeof(double));
  double* alpha_g_cm2_MeV_double                = (double*)calloc(number_of_materials_long,sizeof(double));
  double* p_MeV_double                          = (double*)calloc(number_of_materials_long,sizeof(double));
  double* m_g_cm2_double                        = (double*)calloc(number_of_materials_long,sizeof(double));
  double* average_A_double                      = (double*)calloc(number_of_materials_long,sizeof(double));
  double* average_Z_double                      = (double*)calloc(number_of_materials_long,sizeof(double));

  long i;
  for(i = 0 ; i < number_of_materials_long ; i++){
    material_no_long[i]            = (long)material_no[i];
   }

  AT_get_materials_data( number_of_materials_long,
      material_no_long,
      density_g_cm3_double,
      electron_density_m3_double,
      I_eV_double,
      alpha_g_cm2_MeV_double,
      p_MeV_double,
      m_g_cm2_double,
      average_A_double,
      average_Z_double);

  for(i = 0 ; i < number_of_materials_long ; i++){
    density_g_cm3[i]        = (float)density_g_cm3_double[i];
    electron_density_m3[i]  = (float)electron_density_m3_double[i];
    I_eV[i]                 = (float)I_eV_double[i];
    alpha_g_cm2_MeV[i]      = (float)alpha_g_cm2_MeV_double[i];
    p_MeV[i]                = (float)p_MeV_double[i];
    m_g_cm2[i]              = (float)m_g_cm2_double[i];
    average_A[i]            = (float)average_A_double[i];
    average_Z[i]            = (float)average_Z_double[i];
  }

  free(material_no_long);
  free(density_g_cm3_double);
  free(electron_density_m3_double);
  free(I_eV_double);
  free(alpha_g_cm2_MeV_double);
  free(p_MeV_double);
  free(m_g_cm2_double);
  free(average_A_double);
  free(average_Z_double);
}


void AT_CSDA_range_g_cm2_R(  const int* number_of_particles,
    const float*   E_MeV_u,
    const int*     particle_no,
    const int*     material_no,
    float*         CSDA_range_g_cm2)
{
  const long number_of_particles_long          = (const long)*number_of_particles;
  const long material_no_long                  = (const long)*material_no;

  long* particle_no_long                  = (long*)calloc(number_of_particles_long,sizeof(long));
  double* E_MeV_u_double                  = (double*)calloc(number_of_particles_long,sizeof(double));
  double* CSDA_range_g_cm2_double         = (double*)calloc(number_of_particles_long,sizeof(double));

  long i;
  for(i = 0 ; i < number_of_particles_long ; i++){
    particle_no_long[i]            = (long)particle_no[i];
    E_MeV_u_double[i]              = (double)E_MeV_u[i];
   }

  AT_CSDA_range_g_cm2(  number_of_particles_long,
    E_MeV_u_double,
    particle_no_long,
    material_no_long,
    CSDA_range_g_cm2_double);

  for(i = 0 ; i < number_of_particles_long ; i++){
    CSDA_range_g_cm2[i]        = (float)CSDA_range_g_cm2_double[i];
   }

  free(particle_no_long);
  free(E_MeV_u_double);
  free(CSDA_range_g_cm2_double);
}


void AT_A_from_particle_no_R( const int*  n,
    const int* particle_no,
    int*  A)
{
  long i;
  for (i = 0; i < (const long)*n; i++){
    A[i]  =  (int)AT_A_from_particle_no_single((const long)particle_no[i]);
  }
}


void AT_Z_from_particle_no_R( const int*  n,
    const int* particle_no,
    int*  Z)
{
  long i;
  for (i = 0; i < (const long)*n; i++){
    Z[i]  =  (int)AT_Z_from_particle_no_single((const long)particle_no[i]);
  }
}


void AT_total_u_R(    const int * n,
                const float * E_MeV_u,
                const int   * particle_no,
                const float * fluence_cm2,
                const int   * material_no,
                const int   * er_model,
                float *       u)
{

  //TODO solve negative fluence problem !

  /* int -> long conversion */
  const long n_long = (long)(*n);
  long * particle_no_long  = (long*)calloc(n_long,sizeof(long));
  long i;
  for(i = 0 ; i < n_long ; i++){
    particle_no_long[i]     =  (long)particle_no[i];
  }
  long material_no_long = (long)(*material_no);
  long er_model_long = (long)(*er_model);

  double * E_MeV_u_double  = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i]     =  (double)E_MeV_u[i];
  }

  double*  fluence_cm2_local    =  (double*)calloc(n_long, sizeof(double));
  if(fluence_cm2[0] < 0){
    double*  dose_Gy_local      =  (double*)calloc(n_long, sizeof(double));
    for (i = 0; i < n_long; i++){
      dose_Gy_local[i] = -1.0 * fluence_cm2[i];
    }
    // convert dose to fluence
    AT_fluence_cm2(  n_long,
        E_MeV_u_double,
        particle_no_long,
        dose_Gy_local,
        material_no_long,
        fluence_cm2_local);
    free( dose_Gy_local );
  }else{
    for (i = 0; i < n_long; i++){
      fluence_cm2_local[i] = fluence_cm2[i];
    }
  }

  /* float -> double conversion */
  double * fluence_cm2_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    fluence_cm2_double[i] =  (double)fluence_cm2_local[i];
  }

  free( fluence_cm2_local );

  *u = (float)AT_total_u(   n_long,
                  E_MeV_u_double,
                  particle_no_long,
                  fluence_cm2_double,
                  material_no_long,
                  er_model_long);

  free( particle_no_long );
  free( E_MeV_u_double );
  free( fluence_cm2_double );
}

void AT_convert_beam_parameters_R(  const int*  n,
    float* fluence_cm2,
    float* sigma_cm,
    float* N,
    float* FWHM_mm)
{
	  long i;

	  const long n_long = (long)(*n);

	  double * fluence_cm2_double  	= (double*)calloc(n_long,sizeof(double));
	  double * sigma_cm_double  	= (double*)calloc(n_long,sizeof(double));
	  double * N_double  			= (double*)calloc(n_long,sizeof(double));
	  double * FWHM_mm_double  		= (double*)calloc(n_long,sizeof(double));
	  for(i = 0 ; i < n_long ; i++){
		  fluence_cm2_double[i]     	=  (double)fluence_cm2[i];
		  sigma_cm_double[i]     		=  (double)sigma_cm[i];
		  N_double[i]     				=  (double)N[i];
		  FWHM_mm_double[i]     		=  (double)FWHM_mm[i];
	  }

	AT_convert_beam_parameters(  n_long,
			fluence_cm2_double,
			sigma_cm_double,
			N_double,
			FWHM_mm_double);

	  for(i = 0 ; i < n_long ; i++){
		  fluence_cm2[i]     	=  (float)fluence_cm2_double[i];
		  sigma_cm[i]    		=  (float)sigma_cm_double[i];
		  N[i]     				=  (float)N_double[i];
		  FWHM_mm[i]    		=  (float)FWHM_mm_double[i];
	  }

	free(fluence_cm2_double);
	free(sigma_cm_double);
	free(N_double);
	free(FWHM_mm_double);
}



void AT_GSM_calculate_dose_histogram_R( const int*  number_of_field_components,
    const float*   E_MeV_u,
    const float*   fluence_cm2,
    const int*     particle_no,
    const int*     material_no,
    const int*     rdd_model,
    const float*   rdd_parameter,
    const int*     er_model,
    const int*     nX,
    const float*   pixel_size_m,
    const int*		N_runs,
    const int*     number_of_bins,
    const float*   dose_bin_centers_Gy,
    float *       zero_dose_fraction,
    float *       dose_frequency_Gy){

	  long i;

	  const long n_long = (long)(*number_of_field_components);

	  double * 	E_MeV_u_double  			= (double*)calloc(n_long,sizeof(double));
	  double * 	fluence_cm2_double  		= (double*)calloc(n_long,sizeof(double));
	  long*		particle_no_long			= (long*)calloc(n_long,sizeof(long));

	  for(i = 0 ; i < n_long ; i++){
		  E_MeV_u_double[i]     		=  (double)E_MeV_u[i];
		  fluence_cm2_double[i]     	=  (double)fluence_cm2[i];
		  particle_no_long[i]			=  (long)particle_no[i];
	  }

	  const long		material_no_long			= (const long)(*material_no);
	  const long		rdd_model_long				= (const long)(*rdd_model);
	  const long		er_model_long				= (const long)(*er_model);
	  const long		nX_long						= (const long)(*nX);

	  double 		rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
	  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
	    rdd_parameter_double[i] = (double)rdd_parameter[i];
	  }

	  double	pixel_size_m_double			= (double)(*pixel_size_m);
	  long		number_of_bins_long			= (long)(*number_of_bins);
	  double * 	dose_bin_centers_Gy_double 	= (double*)calloc(number_of_bins_long,sizeof(double));
	  double * 	dose_frequency_Gy_double	= (double*)calloc(number_of_bins_long,sizeof(double));
	  double * 	dose_frequency_Gy_tmp		= (double*)calloc(number_of_bins_long,sizeof(double));

	  for(i = 0 ; i < number_of_bins_long ; i++){
		  dose_bin_centers_Gy_double[i]     		=  (double)dose_bin_centers_Gy[i];
	  }

	  double       zero_dose_fraction_double		= 0.0;
	  double       zero_dose_fraction_tmp			= 0.0;

	  long			N_runs_long						= (long)*N_runs;
	  long			j;

	  /* Create and initialize random number generator */
	  gsl_rng * rng  								= gsl_rng_alloc(gsl_rng_taus);
	  gsl_rng_set(rng, 137);
	  unsigned long  random_number_generator_seed	= gsl_rng_get(rng);

	  for (i = 0; i < N_runs_long; i++){
		  AT_GSM_calculate_dose_histogram( n_long,
			E_MeV_u_double,
			fluence_cm2_double,
			particle_no_long,
			material_no_long,
			rdd_model_long,
			rdd_parameter_double,
			er_model_long,
			nX_long,
			pixel_size_m_double,
			number_of_bins_long,
			dose_bin_centers_Gy_double,
			&random_number_generator_seed,
			&zero_dose_fraction_tmp,
			dose_frequency_Gy_tmp);

		  zero_dose_fraction_double		+= zero_dose_fraction_tmp;
		  zero_dose_fraction_tmp		 = 0.0;

		  for (j = 0; j < number_of_bins_long; j++){
			  dose_frequency_Gy_double[j] += dose_frequency_Gy_tmp[j];
			  dose_frequency_Gy_tmp[j]	   = 0;
		  }
	  }

	  *zero_dose_fraction			= (float)(zero_dose_fraction_double / (double)N_runs_long);

	  for(i = 0 ; i < number_of_bins_long ; i++){
		  dose_frequency_Gy[i]     		=  (float)(dose_frequency_Gy_double[i] / (double)N_runs_long);
	  }

	  free(E_MeV_u_double);
	  free(fluence_cm2_double);
	  free(particle_no_long);

	  free(dose_bin_centers_Gy_double);
	  free(dose_frequency_Gy_double);
	  free(dose_frequency_Gy_tmp);

}

void AT_GSM_calculate_multiple_dose_histograms_R( const int*  number_of_field_components,
    const float*   	E_MeV_u,
    const float*   	fluence_cm2,
    const int*     	particle_no,
    const int*     	material_no,
    const int*     	rdd_model,
    const float*   	rdd_parameter,
    const int*     	er_model,
    const int*     	nX,
    const float*   	pixel_size_m,
    const int*		N_runs,
    const int*		N_repetitions,
    const int*     	number_of_bins,
    const float*   	dose_bin_centers_Gy,
    float *   		dose_bin_width_Gy,
    float *       	mean_d_check_Gy,
    float *       	sd_d_check_Gy,
    float *       	mean_zero_dose_fraction,
    float *       	sd_zero_dose_fraction,
    float *       	mean_dose_frequency_Gy,
    float *       	sd_dose_frequency_Gy){

	  long i;

	  const long n_long = (long)(*number_of_field_components);

	  double * 	E_MeV_u_double  			= (double*)calloc(n_long,sizeof(double));
	  double * 	fluence_cm2_double  		= (double*)calloc(n_long,sizeof(double));
	  long*		particle_no_long			= (long*)calloc(n_long,sizeof(long));

	  for(i = 0 ; i < n_long ; i++){
		  E_MeV_u_double[i]     		=  (double)E_MeV_u[i];
		  fluence_cm2_double[i]     	=  (double)fluence_cm2[i];
		  particle_no_long[i]			=  (long)particle_no[i];
	  }

	  const long		material_no_long			= (const long)(*material_no);
	  const long		rdd_model_long				= (const long)(*rdd_model);
	  const long		er_model_long				= (const long)(*er_model);
	  const long		nX_long						= (const long)(*nX);

	  double 		rdd_parameter_double[RDD_MAX_NUMBER_OF_PARAMETERS];
	  for(i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++){
	    rdd_parameter_double[i] = (double)rdd_parameter[i];
	  }

	  double	pixel_size_m_double				= (double)(*pixel_size_m);
	  long		number_of_bins_long				= (long)(*number_of_bins);
	  double * 	dose_bin_centers_Gy_double 		= (double*)calloc(number_of_bins_long,sizeof(double));
	  double * 	dose_bin_widths_Gy_double 		= (double*)calloc(number_of_bins_long,sizeof(double));

	  for(i = 0 ; i < number_of_bins_long ; i++){
		  dose_bin_centers_Gy_double[i]     		=  (double)dose_bin_centers_Gy[i];
	  }


	  long			N_runs_long						= (long)*N_runs;
	  long			N_repetitions_long				= (long)*N_repetitions;

	  double		mean_d_check_Gy_double			= 0;
	  double		sd_d_check_Gy_double			= 0;
	  double		mean_zero_dose_fraction_double	= 0;
	  double		sd_zero_dose_fraction_double	= 0;

	  double * 		mean_dose_frequency_Gy_double 	= (double*)calloc(number_of_bins_long,sizeof(double));
	  double * 		sd_dose_frequency_Gy_double		= (double*)calloc(number_of_bins_long,sizeof(double));

	AT_GSM_calculate_multiple_dose_histograms( n_long,
	    E_MeV_u_double,
	    fluence_cm2_double,
	    particle_no_long,
	    material_no_long,
	    rdd_model_long,
	    rdd_parameter_double,
	    er_model_long,
	    nX_long,
	    pixel_size_m_double,
	    N_runs_long,
	    N_repetitions_long,
	    number_of_bins_long,
	    dose_bin_centers_Gy_double,
	    dose_bin_widths_Gy_double,
	    &mean_d_check_Gy_double,
	    &sd_d_check_Gy_double,
	    &mean_zero_dose_fraction_double,
	    &sd_zero_dose_fraction_double,
	    mean_dose_frequency_Gy_double,
	    sd_dose_frequency_Gy_double);

	  *mean_d_check_Gy					= (float)mean_d_check_Gy_double;
	  *sd_d_check_Gy					= (float)sd_d_check_Gy_double;
	  *mean_zero_dose_fraction			= (float)mean_zero_dose_fraction_double;
	  *sd_zero_dose_fraction			= (float)sd_zero_dose_fraction_double;

	  for(i = 0 ; i < number_of_bins_long ; i++){
		  dose_bin_width_Gy[i]     		=  	(float)dose_bin_widths_Gy_double[i];
		  mean_dose_frequency_Gy[i]		=	(float)mean_dose_frequency_Gy_double[i];
		  sd_dose_frequency_Gy[i]		=	(float)sd_dose_frequency_Gy_double[i];
	  }

	  free(E_MeV_u_double);
	  free(fluence_cm2_double);
	  free(particle_no_long);

	  free(dose_bin_centers_Gy_double);
	  free(dose_bin_widths_Gy_double);
	  free(mean_dose_frequency_Gy_double);
	  free(sd_dose_frequency_Gy_double);

}

void AT_fluence_weighted_stopping_power_ratio_R( const int*     n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    const int*	  reference_material_no,
    float*		  fluence_weighted_stopping_power_ratio){

	long i;

	const long n_long = (long)(*n);

	double * 	E_MeV_u_double  			= (double*)calloc(n_long,sizeof(double));
	double * 	fluence_cm2_double  		= (double*)calloc(n_long,sizeof(double));
	long*		particle_no_long			= (long*)calloc(n_long,sizeof(long));

	for(i = 0 ; i < n_long ; i++){
		E_MeV_u_double[i]     		=  (double)E_MeV_u[i];
		fluence_cm2_double[i]     	=  (double)fluence_cm2[i];
		particle_no_long[i]			=  (long)particle_no[i];
	}

	const long		material_no_long			= (const long)(*material_no);
	const long		reference_material_no_long	= (const long)(*reference_material_no);

	double fluence_weighted_stopping_power_ratio_double = AT_fluence_weighted_stopping_power_ratio( n_long,
	    E_MeV_u_double,
	    particle_no_long,
	    fluence_cm2_double,
	    material_no_long,
	    reference_material_no_long);

	*fluence_weighted_stopping_power_ratio = (float)fluence_weighted_stopping_power_ratio_double;

	free(E_MeV_u_double);
	free(fluence_cm2_double);
	free(particle_no_long);
}

void AT_Bethe_Mass_Stopping_Power_MeV_cm2_g_R(	const int* n,
		const float* E_MeV_u,
		const int* particle_no,
		const int* material_no,
		const float* E_restricted_keV,
		float* Mass_Stopping_Power_MeV_cm2_g){

	long i;

	const long n_long = (long)(*n);

	double * 	E_MeV_u_double  			= (double*)calloc(n_long,sizeof(double));
	long*		particle_no_long			= (long*)calloc(n_long,sizeof(long));

	for(i = 0 ; i < n_long ; i++){
		E_MeV_u_double[i]     		=  (double)E_MeV_u[i];
		particle_no_long[i]			=  (long)particle_no[i];
	}

	const long		material_no_long			= (const long)(*material_no);
	const double	E_restricted_keV_double		= (const double)(*E_restricted_keV);

	double * 		Mass_Stopping_Power_MeV_cm2_g_double  	= (double*)calloc(n_long,sizeof(double));

	AT_Bethe_Mass_Stopping_Power_MeV_cm2_g( n_long,
	    E_MeV_u_double,
	    particle_no_long,
	    material_no_long,
	    E_restricted_keV_double,
	    Mass_Stopping_Power_MeV_cm2_g_double);

	for(i = 0 ; i < n_long ; i++){
		Mass_Stopping_Power_MeV_cm2_g[i]     		=  (float)Mass_Stopping_Power_MeV_cm2_g_double[i];
	}

	free(E_MeV_u_double);
	free(particle_no_long);
	free(Mass_Stopping_Power_MeV_cm2_g_double);
}
