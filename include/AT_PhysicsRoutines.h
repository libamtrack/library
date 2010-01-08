#ifndef AT_PHYSICSROUTINES_H_
#define AT_PHYSICSROUTINES_H_

/**
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
#include "AT_Data.h"
#include "AT_Utils.h"
#include "AT_ElectronRange.h"

void AT_beta_from_mass(  long*  n,
              float*  E_MeV_u,
              float*  mass,
              float*  beta);

void AT_E_from_beta_and_mass(  long*  n,
    float*  beta,
    float*  mass,
    float*  E_MeV_u);

void AT_effective_charge_from_beta(  long*  n,
                    float*  beta,
                    long*  Z,
                    float*  effective_charge);

void AT_beta_from_particle_no(  long*  n,
                  float*  E_MeV_u,
                  long*  particle_no,
                  float*  beta);


void AT_Bohr_Energy_Straggling_g_cm2(  long*  n,
                    char**  material_name,
                    float*  dsE2dz);

void AT_effective_charge_from_particle_no(  long*  n,
                        float*  E_MeV_u,
                        long*  particle_no,
                        float*  effective_charge);

void AT_scaled_energy(  long*  n,
            float*  E_MeV_u,
            long*  particle_no,
            float*  effective_charge);

void AT_E_MeV_u_from_scaled_energy(  long*  n,
            float*  E_MeV_u,
            long*  particle_no,
            float*  scaled_energy);

void AT_max_E_transfer_MeV(  long*  n,
                float*  E_MeV_u,
                long*  particle_no,
                float*  max_E_transfer_MeV);


void AT_Z_from_particle_no(  long*  n,                // should rather be AT_ParticleProperties(particle.index)
                long*  particle_no,
                long*  Z);


#endif /* AT_PHYSICSROUTINES_H_ */
