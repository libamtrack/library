/**
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


void AT_Z_from_particle_no(  long*  n,
    long*  particle_no,
    long*  Z)
{
  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));

  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  // loop over n to find Z for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    Z[i]  = AT_Particle_Data.Z[matches[i]];}

  free(matches);
}
