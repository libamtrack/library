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
 *    Author: kongruencja
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

#include "AT_Wrapper_R.h"


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


void AT_gamma_response_R( const int*  n,
    const float*  d_Gy,
    const int*    gamma_model,
    const float*  gamma_parameter,
    float*        S){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long gamma_model_long = (long)(*gamma_model);

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

void AT_run_GSM_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameters,
    const int*    er_model,
    const int*    gamma_model,
    const float*  gamma_parameters,
    const int*    N_runs,
    const int*    N2,
    const float*  fluence_factor,
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
  const long N2_long = (long)(*N2);
  const long nX_long = (long)(*nX);

  /* int -> bool conversion */
  const bool write_output_bool = (bool)(*write_output);
  const bool lethal_events_mode_bool = (bool)(*lethal_events_mode);

  /* float -> double conversion */
  double * E_MeV_u_double = (double*)calloc(n_long,sizeof(double));
  double * fluence_cm2_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
    fluence_cm2_double[i] = (double)fluence_cm2[i];
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
  double voxel_size_m_double = (double)(*voxel_size_m);

  /* place for results */
  double results_double[10];

  AT_run_GSM_method(n_long,
      E_MeV_u_double,
      particle_no_long,
      fluence_cm2_double,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      gamma_model_long,
      gamma_parameter_double,
      N_runs_long,
      N2_long,
      fluence_factor_double,
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
  free(fluence_cm2_double);
}


void AT_run_SPIFF_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
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
  double * fluence_cm2_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
    fluence_cm2_double[i] = (double)fluence_cm2[i];
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
      fluence_cm2_double,
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
   free(fluence_cm2_double);
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

  AT_SC_get_f1_array_size( n_long,
      E_MeV_u_double,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      N2_long,
      &n_bins_f1_long,
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
    const float*  fluence_cm2,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    const int*    N2,
    const int*    n_bins_f1,
    const float*  f1_parameters,
    float*        norm_fluence,
    float*        dose_contribution_Gy,
    float*        f_parameters,
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
  double * fluence_cm2_double = (double*)calloc(n_long,sizeof(double));
  for(i = 0 ; i < n_long ; i++){
    E_MeV_u_double[i] = (double)E_MeV_u[i];
    fluence_cm2_double[i] = (double)fluence_cm2[i];
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
  double * norm_fluence_double = (double*)calloc(n_long,sizeof(double));
  double * dose_contribution_Gy_double = (double*)calloc(n_long,sizeof(double));
  double f_parameters_double[AT_SC_F_PARAMETERS_LENGTH];
  double * f1_d_Gy_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  double * f1_dd_Gy_double = (double*)calloc(n_bins_f1_long,sizeof(double));
  double * f1_double = (double*)calloc(n_bins_f1_long,sizeof(double));

  AT_SC_get_f1( n_long,
      E_MeV_u_double,
      particle_no_long,
      fluence_cm2_double,
      material_no_long,
      rdd_model_long,
      rdd_parameter_double,
      er_model_long,
      N2_long,
      n_bins_f1_long,
      f1_parameters_double,
      norm_fluence_double,
      dose_contribution_Gy_double,
      f_parameters_double,
      f1_d_Gy_double,
      f1_dd_Gy_double,
      f1_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    norm_fluence[i] = (float)norm_fluence_double[i];
    dose_contribution_Gy[i] = (float)dose_contribution_Gy_double[i];
  }
  for(i = 0 ; i < AT_SC_F_PARAMETERS_LENGTH ; i++){
    f_parameters[i] = (float)f_parameters_double[i];
  }
  for(i = 0 ; i < n_bins_f1_long ; i++){
    f1_d_Gy[i] = (float)f1_d_Gy_double[i];
    f1_dd_Gy[i] = (float)f1_dd_Gy_double[i];
    f1[i] = (float)f1_double[i];
  }
  free(particle_no_long);
  free(E_MeV_u_double);
  free(fluence_cm2_double);
  free(norm_fluence_double);
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

