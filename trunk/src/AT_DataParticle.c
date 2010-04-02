/**
 * @file
 * @brief ...
 */


/*
*    AT_DataParticle.c
*    ==============
*
*    Created on: 09.01.2010
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

#include "AT_DataParticle.h"


inline long AT_A_from_particle_no_single(  const long  particle_no ){
  // TODO maybe we could use modulo division here ? A = particle_no % 1000 ?
  long A = particle_no / 1000;
  return particle_no - A * 1000;
}


int AT_A_from_particle_no( const long  n,
    const long  particle_no[],
    long  A[])
{
  long i;
  for (i = 0; i < n; i++){
    A[i]  =  AT_A_from_particle_no_single(particle_no[i]);
  }
  return 0;
}


inline long AT_Z_from_particle_no_single(  const long  particle_no ){
  return particle_no / 1000;
}


int AT_Z_from_particle_no( const long  n,
    const long  particle_no[],
    long  Z[])
{
  long i;
  for (i = 0; i < n; i++){
    Z[i]  =  AT_Z_from_particle_no_single(particle_no[i]);
  }
  return 0;
}


int AT_atomic_weight_from_particle_no( const long  n,
    const long  particle_no[],
    float  atomic_weight[])
{
  long i;
  long*  matches  =  (long*)calloc(n, sizeof(long));
  long*  Z        =  (long*)calloc(n, sizeof(long));

  AT_Z_from_particle_no( n,
      particle_no,
      Z);

  long   n_particle_data = PARTICLE_DATA_N;

  find_elements_int(  Z,
      n,
      AT_Particle_Data.Z,
      n_particle_data,
      matches);

  for (i = 0; i < n; i++){
    atomic_weight[i]    = AT_Particle_Data.atomic_weight[matches[i]];
  }

  free(Z);
  free(matches);
  return 0;
}

/*
void AT_Particle_Properties(  const long*  n,
    const long*  particle_no,
    long*   Z,
    long*   A,
    char**  element_name,
    char**  element_acronym,
    float** density_g_cm3,
    float** I_eV)
{

  // find look-up indices for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));

  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  // loop over n to find data for all given particles
  long  i;
  if( particle_name != NULL ){
    for(i = 0; i < *n; i++){
      strcpy(particle_name[i], AT_Particle_Data.particle_name[matches[i]]);
    }
  }
  if( USRTRACK_name != NULL ){
    for(i = 0; i < *n; i++){
      strcpy(USRTRACK_name[i], AT_Particle_Data.USRTRACK_name[matches[i]]);
    }
  }
  if( element_name != NULL ){
    for(i = 0; i < *n; i++){
      strcpy(element_name[i], AT_Particle_Data.element_name[matches[i]]);
    }
  }
  if( mass != NULL ){
    for(i = 0; i < *n; i++){
      mass[i] = AT_Particle_Data.mass[matches[i]];
    }
  }
  if( A != NULL ){
    for(i = 0; i < *n; i++){
      A[i] = AT_Particle_Data.A[matches[i]];
    }
  }
  if( Z != NULL ){
    for(i = 0; i < *n; i++){
      Z[i] = AT_Particle_Data.Z[matches[i]];
    }
  }
  free(matches);
}
*/

