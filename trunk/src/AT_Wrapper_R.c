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

void AT_max_E_transfer_MeV_R(  const int*  n,
    const float*  E_MeV_u,
    // results
    float*  max_E_transfer_MeV){

  long  n_R           = (long)*n;

  AT_max_E_transfer_MeV(  &n_R, E_MeV_u, max_E_transfer_MeV);

}

void AT_max_electron_range_m_R(  const int*  n,
    const float*  E_MeV_u,
    const int*  material_no,
    const int*   er_model,
    float*  max_electron_range_m)
{
  const long n_long = (long)(*n);

  AT_max_electron_range_m( n_long,
      E_MeV_u,
      *material_no,
      *er_model,
      max_electron_range_m);

}


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
  long i;
  long * particle_no_long = (long*)calloc(*n,sizeof(long));
  for(i = 0 ; i < *n ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  AT_D_RDD_Gy( &n_long,
      r_m,
      E_MeV_u,
      particle_no_long,
      &material_no_long,
      &rdd_model_long,
      rdd_parameter,
      &er_model_long,
      er_parameter,
      D_RDD_Gy);

  free(particle_no_long);
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
  long i;
  long * particle_no_long = (long*)calloc(*n,sizeof(long));
  for(i = 0 ; i < *n ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  AT_r_RDD_m( &n_long,
      D_RDD_Gy,
      E_MeV_u,
      particle_no_long,
      &material_no_long,
      &rdd_model_long,
      rdd_parameter,
      &er_model_long,
      er_parameter,
      r_RDD_m);

  free(particle_no_long);
}

void AT_RDD_ExtendedTarget_Gy_R( const int*  n,
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
