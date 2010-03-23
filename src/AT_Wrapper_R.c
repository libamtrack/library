/**
 * @file
 * @brief Wrapper functions
 *
 * C functions which are called from R cannot have input
 * integer parameters of type "long". Only "int" type is
 * accepted. This file contains set of wrapper functions,
 * which are casting int arguments to long if necessary.
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
    /* radiation field parameters */
    const float*  E_MeV_u,
    const int*  particle_no,
    /* detector parameters */
    const int*  material_no,
    /* radial dose distribution model */
    const int*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const int*  er_model,
    const float*  er_parameter,
    float*  D_RDD_Gy){

  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long er_model_long = (long)(*er_model);
  const long rdd_model_long = (long)(*rdd_model);
  const long particle_no_long = (long)(*particle_no);

  AT_D_RDD_Gy( n_long,
      r_m,
      *E_MeV_u,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter,
      er_model_long,
      er_parameter,
      D_RDD_Gy);

}

void AT_D_RDD_ExtendedTarget_Gy_R( const int*  n,
    const float* r_m,
    const float*  a0_m,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const int*    particle_no,
    /* detector parameters */
    const int*    material_no,
    /* radial dose distribution model */
    const int*    rdd_model,
    const float* rdd_parameter,
    /* electron range model */
    const int*    er_model,
    const float* er_parameter,
    float*       D_RDD_Gy){

  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long er_model_long = (long)(*er_model);
  const long rdd_model_long = (long)(*rdd_model);
  const long particle_no_long = (long)(*particle_no);

  AT_RDD_ExtendedTarget_Gy(n_long,r_m,*a0_m,*E_MeV_u,particle_no_long,material_no_long,rdd_model_long,rdd_parameter,er_model_long,er_parameter,D_RDD_Gy);

}

void AT_LET_MeV_cm2_g_R(  const int*  n,
    const float*  E_MeV_u,
    const int*   particle_no,
    const int*   material_no,
    float*       LET_MeV_cm2_g){

  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  long i;
  long * particle_no_long = (long*)calloc(*n,sizeof(long));
  for(i = 0 ; i < *n ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  AT_LET_MeV_cm2_g(  n_long,
      E_MeV_u,
      particle_no_long,
      material_no_long,
      LET_MeV_cm2_g);
}

void AT_gamma_response_R( const int*  n,
    const float*  d_Gy,
    const int*  gamma_model,
    const float*  gamma_parameter,
    float*  S){

  const long n_long = (long)(*n);
  const long gamma_model_long = (long)(*gamma_model);

  AT_gamma_response(n_long,d_Gy,gamma_model_long,gamma_parameter,S);

}

void AT_max_E_transfer_MeV_R(  const int*  n,
    const float*  E_MeV_u,
    // results
    float*  max_E_transfer_MeV){

  long  n_R           = (long)*n;
  AT_max_E_transfer_MeV(  n_R, E_MeV_u, max_E_transfer_MeV);

}

void AT_max_electron_ranges_m_R(  const int*  number_of_particles,
    const float*  E_MeV_u,
    const int*    material_no,
    const int*    er_model,
    float*        max_electron_range_m)
{
  const long number_of_particles_long = (long)(*number_of_particles);

  AT_max_electron_ranges_m( number_of_particles_long,
      E_MeV_u,
      *material_no,
      *er_model,
      max_electron_range_m);

}

void AT_r_RDD_m_R  ( const int*  n,
    const float*  D_RDD_Gy,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const int*  particle_no,
    /* detector parameters */
    const int*  material_no,
    /* radial dose distribution model */
    const int*  rdd_model,
    const float*  rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
    /* electron range model */
    const int*  er_model,
    const float*  er_parameter,
    float*  r_RDD_m){

  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long er_model_long = (long)(*er_model);
  const long rdd_model_long = (long)(*rdd_model);
  const long particle_no_long = (long)(*particle_no);

  AT_r_RDD_m( n_long,
      D_RDD_Gy,
      *E_MeV_u,
      particle_no_long,
      material_no_long,
      rdd_model_long,
      rdd_parameter,
      er_model_long,
      er_parameter,
      r_RDD_m);

}

void  AT_SC_get_f1_array_size_R(  /* radiation field parameters */
    const int*  n,
    const float*  E_MeV_u,
    const int*  particle_no,
    /* detector parameters */
    const int*  material_no,
    const int*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const int*  er_model,
    const float*  er_parameter,
    /* algorithm parameters*/
    const int*  N2,
    // from here: return values
    int*  n_bins_f1,
    float*  f1_parameters){

  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long material_no_long = (long)(*material_no);
  const long N2_long = (long)(*N2);

  long i;
  long * particle_no_long = (long*)calloc(*n,sizeof(long));
  for(i = 0 ; i < *n ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  long n_bins_f1_long;

  AT_SC_get_f1_array_size( &n_long,
      E_MeV_u,
      particle_no_long,
      &material_no_long,
      &rdd_model_long,
      rdd_parameter,
      &er_model_long,
      er_parameter,
      &N2_long,
      &n_bins_f1_long,
      f1_parameters);

  *n_bins_f1 = (int)n_bins_f1_long;

  free(particle_no_long);
}

void  AT_SC_get_f1_R(  /* radiation field parameters */
    const int*  n,
    const float*  E_MeV_u,
    const int*  particle_no,
    const float*  fluence_cm2,
    /* detector parameters */
    const int*  material_no,
    const int*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const int*  er_model,
    const float*  er_parameter,
    /* algorithm parameters*/
    const int*  N2,
    const int*  n_bins_f1,
    /* f1 parameters*/
    const float*  f1_parameters,
    // from here: return values
    float*  norm_fluence,
    float*  dose_contribution_Gy,
    float*  f_parameters,
    float*  f1_d_Gy,
    float*  f1_dd_Gy,
    float*  f1){

  const long n_long = (long)(*n);
  const long rdd_model_long = (long)(*rdd_model);
  const long er_model_long = (long)(*er_model);
  const long material_no_long = (long)(*material_no);
  const long N2_long = (long)(*N2);
  const long n_bins_f1_long = (long)(*n_bins_f1);

  long i;
  long * particle_no_long = (long*)calloc(*n,sizeof(long));
  for(i = 0 ; i < *n ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  AT_SC_get_f1( &n_long,
      E_MeV_u,
      particle_no_long,
      fluence_cm2,
      &material_no_long,
      &rdd_model_long,
      rdd_parameter,
      &er_model_long,
      er_parameter,
      &N2_long,
      &n_bins_f1_long,
      f1_parameters,
      norm_fluence,
      dose_contribution_Gy,
      f_parameters,
      f1_d_Gy,
      f1_dd_Gy,
      f1);

  free(particle_no_long);
}
