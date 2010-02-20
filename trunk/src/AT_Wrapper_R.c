/**
 * @file
 * @brief Wrapper functions
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
    const int*  particle_no,
    // results
    float*  max_E_transfer_MeV){

  long  n_R           = (long)*n;
  long* particle_no_R = (long*)calloc(n_R, sizeof(long));

  int i;
  for( i = 0 ; i < *n ; i++ ){
    particle_no_R[i] = particle_no[i];
  }

  AT_max_E_transfer_MeV(  &n_R, E_MeV_u, particle_no_R , max_E_transfer_MeV);

  free(particle_no_R);
}

void AT_max_electron_range_m_R(  const int*  n,
    const float*  E_MeV_u,
    const int*  particle_no,
    const int*  material_no,
    const int*   er_model,
    float*  max_electron_range_m)
{
  const long n_long = (long)(*n);
  const long material_no_long = (long)(*material_no);
  const long er_model_long = (long)(*er_model);

  long i;
  long * particle_no_long = (long*)calloc(*n,sizeof(long));
  for(i = 0 ; i < *n ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  AT_max_electron_range_m( &n_long,
      E_MeV_u,
      particle_no_long,
      &material_no_long,
      &er_model_long,
      max_electron_range_m);

  free(particle_no_long);

}
