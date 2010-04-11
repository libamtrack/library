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
  long A = particle_no % 1000;
  if( (1 <= A) && (A <= 300)){
    return A;
  } else {
    printf( "Wrong particle number %ld, please provide it in correct format (XXXYYY, where XXX is Z and YYY is A)\n", particle_no);
    return 1;
  }
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
  long Z = particle_no / 1000;
  if( (1 <= Z) && (Z <= 118) ){
    return Z;
  } else {
    printf( "Wrong particle number %ld, please provide it in correct format (XXXYYY, where XXX is Z and YYY is A)\n", particle_no);
    return 1;
  }
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
    double  atomic_weight[])
{
  long i;
  long*  matches  =  (long*)calloc(n, sizeof(long));
  long*  Z        =  (long*)calloc(n, sizeof(long));

  AT_Z_from_particle_no( n,
      particle_no,
      Z);

  find_elements_int(  Z,
      n,
      AT_Particle_Data.Z,
      AT_Particle_Data.n,
      matches);

  for (i = 0; i < n; i++){
    atomic_weight[i]    = AT_Particle_Data.atomic_weight[matches[i]];
  }

  free(Z);
  free(matches);
  return 0;
}



int AT_particle_name_from_particle_no( const long  n,
    const long  particle_no[],
    char  particle_name[][PARTICLE_NAME_NCHAR])
{
  long i;
  long*  matches  =  (long*)calloc(n, sizeof(long));
  long*  Z        =  (long*)calloc(n, sizeof(long));
  long*  A        =  (long*)calloc(n, sizeof(long));

  AT_Z_from_particle_no( n,
      particle_no,
      Z);

  AT_A_from_particle_no( n,
      particle_no,
      A);

  find_elements_int(  Z,
      n,
      AT_Particle_Data.Z,
      AT_Particle_Data.n,
      matches);

  for (i = 0; i < n; i++){
    ltoa(A[i], particle_name[i], 10);
    strcat(particle_name[i], AT_Particle_Data.element_acronym[matches[i]]);
  }



  free(A);
  free(Z);
  free(matches);
  return 0;
}
