/**
 * @brief Particle
 */


/*
 *    AT_DataParticle.c
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

#include "AT_DataParticle.h"

 long AT_particle_no_from_Z_and_A_single(  const long  Z,
		const long  A){
	assert((1 <= A) && (A <= 300));
	assert((1 <= Z) && (Z <= 118));
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
  return AT_Success;
}

 long AT_A_from_particle_no_single(  const long  particle_no ){
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
    if( A[i] < 0 ){
      return AT_Particle_Not_Defined;
    }
  }
  return AT_Success;
}


 long AT_Z_from_particle_no_single(  const long  particle_no ){
  long Z = particle_no / 1000;
  if( (1 <= Z) && (Z <= 118) ){
    return Z;
  } else {
    printf( "Wrong particle number %ld, please provide it in correct format \n(XXXYYY, where XXX is Z (from 1 to 118) and YYY is A (from 1 to 300)\n", particle_no);
    return -1;
  }
}


int AT_Z_from_particle_no( const long  n,
    const long  particle_no[],
    long  Z[])
{
  long i;
  for (i = 0; i < n; i++){
    Z[i]  =  AT_Z_from_particle_no_single(particle_no[i]);
    if( Z[i] < 0){
      return AT_Particle_Not_Defined;
    }
  }
  return AT_Success;
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
  return AT_Success;
}

int AT_I_eV_from_particle_no( const long  n,
    const long  particle_no[],
    double  I_eV[])
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
    I_eV[i]    = AT_Particle_Data.I_eV_per_Z[matches[i]] * Z[i];
  }

  free(Z);
  free(matches);
  return AT_Success;
}

int AT_nuclear_spin_from_particle_no_multi( const long  n,
    const long  particle_no[],
    double  I[]){

	long i;
	for( i = 0; i < n; i++){
		I[i]	= AT_nuclear_spin_from_particle_no_single(particle_no[i]);
	}

	return EXIT_SUCCESS;
}

double AT_nuclear_spin_from_particle_no_single( const long  particle_no){

	return AT_nuclear_spin_from_Z_and_A( AT_Z_from_particle_no_single(particle_no),
			AT_A_from_particle_no_single(particle_no));
}

double AT_nuclear_spin_from_Z_and_A( const long  Z,
    const long  A){

	if(A % 2 == 0){
		if(Z % 2 == 0){
			return 0.0;		// If even number of neutrons AND protons -> I = 0
		}else{
			return 1.0;     // If odd number of n AND p but even number of nucleons -> I = 1
		}
	}else{
		return 0.5;			// If odd number of nucleons -> I = 1/2
	}
}

int AT_particle_name_from_particle_no_single( const long  particle_no,
    char * particle_name){

	  long  Z = AT_Z_from_particle_no_single(  particle_no );
	  long  A = AT_A_from_particle_no_single(  particle_no );

	  long  match;

	  find_elements_int(  &Z,
			  1,
			  AT_Particle_Data.Z,
			  AT_Particle_Data.n,
			  &match);

	  sprintf(particle_name, "%ld", A);
	  if( match >= 0 ){
		  strcat(particle_name, AT_Particle_Data.element_acronym[match]);
	  } else {
		  const char * unknown_acronym = "??";
		  strcat(particle_name, unknown_acronym);
	  }

	  return AT_Success;
}


long AT_particle_no_from_particle_name_single( const char particle_name[PARTICLE_NAME_NCHAR]){
	assert( particle_name != NULL);

	char * literal_part = 0;
	long A = strtol(particle_name,&literal_part,10);
	if( (A == 0) && (*literal_part != 0)){
		return -1;
	}
	if( A == 0){
		return -1;
	}

	long match;
	find_elements_char( (const char**)(&literal_part),
			1,
			AT_Particle_Data.element_acronym,
			PARTICLE_DATA_N,
			&match);

	if( match == -1){
		return -1;
	}

	long Z = AT_Particle_Data.Z[match];

	return AT_particle_no_from_Z_and_A_single(Z, A);
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
  return AT_Success;
}


// TODO some comments needed
int AT_particle_no_from_particle_name( const long  n,
    char * particle_name[],
    long particle_no[]){

  assert( particle_name != NULL);
  long i;
  for (i = 0; i < n; i++){
	particle_no[i] = AT_particle_no_from_particle_name_single( particle_name[i] );
	if( particle_no[i] < 0)
		return AT_Particle_Not_Defined;
  }
  return AT_Success;
}
