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

void AT_mass_from_particle_no( const long*  n,
    const long*  particle_no,
    float*  mass)
{
  AT_Particle_Properties(n,particle_no,NULL,NULL,NULL,NULL,NULL,mass);
}

void AT_A_from_particle_no( const long*  n,
    const long*  particle_no,
    long*  A)
{
  AT_Particle_Properties(n,particle_no,NULL,NULL,NULL,NULL,A,NULL);
}

void AT_Z_from_particle_no( const long*  n,
    const long*  particle_no,
    long*  Z)
{
  AT_Particle_Properties(n,particle_no,NULL,NULL,NULL,Z,NULL,NULL);
}

void AT_Particle_Properties(  const long*  n,
    const long*  particle_no,
    /* return values*/
    char**  particle_name,
    char**  USRTRACK_name,
    char**  element_name,
    long*  Z,
    long*  A,
    float*  mass)
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
      strcpy(particle_name[i], AT_Particle_Data.particle_name[i]);
    }
  }
  if( USRTRACK_name != NULL ){
    for(i = 0; i < *n; i++){
      strcpy(USRTRACK_name[i], AT_Particle_Data.USRTRACK_name[i]);
    }
  }
  if( element_name != NULL ){
    for(i = 0; i < *n; i++){
      strcpy(element_name[i], AT_Particle_Data.element_name[i]);
    }
  }
  if( mass != NULL ){
    for(i = 0; i < *n; i++){
      mass[i] = AT_Particle_Data.mass[i];
    }
  }
  if( A != NULL ){
    for(i = 0; i < *n; i++){
      A[i] = AT_Particle_Data.A[i];
    }
  }
  if( Z != NULL ){
    for(i = 0; i < *n; i++){
      Z[i] = AT_Particle_Data.Z[i];
    }
  }
  free(matches);
}


