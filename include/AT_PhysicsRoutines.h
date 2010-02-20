#ifndef AT_PHYSICSROUTINES_H_
#define AT_PHYSICSROUTINES_H_

/**
 * @file
 * @brief Header file for Physics related routines
 */

/*
*    AT_PhysicsRoutines.h
*    ==================
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

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "AT_Constants.h"
#include "AT_DataMaterial.h"
#include "AT_DataParticle.h"
#include "AT_DataLET.h"
#include "AT_ElectronRange.h"
#include "AT_NumericalRoutines.h"


void AT_beta_from_mass( const long*  n,
    const float*  E_MeV_u,
    const float*  mass,
    float*  beta);

void AT_E_from_beta_and_mass(  const long*  n,
    const float*  beta,
    const float*  mass,
    float*  E_MeV_u);

void AT_effective_charge_from_beta(  const long*  n,
    const float*  beta,
    const long*  Z,
    float*  effective_charge);

void AT_beta_from_particle_no(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    // results
    float*  beta);

/**
 * Get Bohr's energy spread (Wilson, 1947, Phys Rev 71, 385)
 */
void AT_Bohr_Energy_Straggling_g_cm2(  const long*  n,
    const long*  material_no,
    float*  dsE2dz);

void AT_effective_charge_from_particle_no(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    float*  effective_charge);

void AT_scaled_energy( const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    float*  scaled_energy);

void AT_E_MeV_u_from_scaled_energy(  const long*  n,
    const float*  scaled_energy,
    const long*  particle_no,
    float*  E_MeV_u);

void AT_max_E_transfer_MeV(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    // results
    float*  max_E_transfer_MeV);

/**
* Returns dose in Gy for each given particle
* @param  n      number of particle types in the mixed particle field (pointer to single variable)
* @param  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
* @param  fluence_cm2     fluence for each particle type (pointer to array of size n)
* @param  particle_no    type of the particles in the mixed particle field (pointer to array of size n)
* @see          AT_DataParticle.h for definition
* @param  material_no  material index
* @see          AT_DataMaterial.h for definition
* @param  D_Gy  pointer to vector of size n to be allocated by the user which will be used to return the results
* @return  none
*/
void AT_D_Gy(  const long*  n,
    const float*  E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    // results
     float* D_Gy);


void AT_interparticleDistance_m(       const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  fluence_cm2,
    float*  results_m
);

void AT_inv_interparticleDistance_Gy(  const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  distance_m,
    float*  results_Gy
);

void AT_inv_interparticleDistance_cm2( const long*   n,
    const float*  distance_m,
    float*  results_cm2
);


#endif /* AT_PHYSICSROUTINES_H_ */
