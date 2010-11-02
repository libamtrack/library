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

inline long AT_particle_no_from_Z_and_A_single(  const long  Z,
    const long  A){
  return 1000 * Z + A;
}

int AT_particle_no_from_Z_and_A( const long  n,
    const long  Z[],
    const long  A[],
    long  particle_no[])
{
  long i;
  for (i = 0; i < n; i++){
    particle_no[i]  =  AT_particle_no_from_Z_and_A_single(Z[i], A[i]);
  }
  return 0;
}

inline long AT_A_from_particle_no_single(  const long  particle_no ){
  long A = particle_no % 1000;
  if( (1 <= A) && (A <= 300)){
    return A;
  } else {
    printf( "Wrong particle number %ld, please provide it in correct format \n(XXXYYY, where XXX is Z (from 1 to 118) and YYY is A (from 1 to 300)\n", particle_no);
    return -1;
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
    printf( "Wrong particle number %ld, please provide it in correct format \n(XXXYYY, where XXX is Z (from 1 to 118) and YYY is A (from 1 to 300)\n", particle_no);    return -1;
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
    char        particle_name[][PARTICLE_NAME_NCHAR])
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
    sprintf(particle_name[i], "%ld", A[i]);
    if( matches[i] >= 0 ){
      strcat(particle_name[i], AT_Particle_Data.element_acronym[matches[i]]);
    } else {
      const char * unknown_acronym = "??";
      strcat(particle_name[i], unknown_acronym);
    }
  }

  free(A);
  free(Z);
  free(matches);
  return 0;
}


// TODO single particle method needed
// TODO some comments needed
int AT_particle_no_from_particle_name( const long  n,
    char * particle_name[],
    long particle_no[]){

  long i;
  char   any_character[30] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  for (i = 0; i < n; i++){
    long split_pos = strcspn(particle_name[i], any_character);

    char A_str[4] = "\0\0\0\0";
    strncpy(A_str, particle_name[i], split_pos);
    long A        = atol(A_str);

    char element_acronym_str[1][PARTICLE_NAME_NCHAR];
    strncpy(element_acronym_str[0], &particle_name[i][split_pos], PARTICLE_NAME_NCHAR - split_pos);
    long  match = 0;
    char** test = (char**)malloc(sizeof(char*));
    *test = (char*)element_acronym_str[0];

    find_elements_char( (const char**)test,
        1,
        AT_Particle_Data.element_acronym,
        PARTICLE_DATA_N,
        &match);

    // TODO check if match != -1 !

    long Z = AT_Particle_Data.Z[match];

    particle_no[i] = AT_particle_no_from_Z_and_A_single(Z, A);

    free(test);
  }
  return 0;
}
