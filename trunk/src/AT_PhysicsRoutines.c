/**
*    AT_PhysicsRoutines.c
*    ==============
*
*    Created on: 8.01.2010
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

#include "AT_PhysicsRoutines.h"


void AT_beta_from_mass( const long*  n,
    const float*  E_MeV_u,
    const float*  mass,
    float*  beta)
{
  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    float  E_MeV    =  E_MeV_u[i] * mass[i];

    // Get rest energy
    float  E0_MeV    =  (float)proton_mass_MeV_c2 * mass[i];

    // Return relativistic speed
    beta[i]        =  (float)sqrt(1 - 1/((1 + E_MeV/E0_MeV)*(1 + E_MeV/E0_MeV)));
  }
}

void AT_E_from_beta_and_mass(  const long*  n,
    const float*  beta,
    const float*  mass,
    float*  E_MeV_u)
{
  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    // Get rest energy
    float  E0_MeV    =  (float)proton_mass_MeV_c2 * mass[i];

    E_MeV_u[i]      =  E0_MeV * (1.0f / (1 - beta[i]*beta[i]) - 1);

    E_MeV_u[i]      /=  mass[i];
  }
}


void AT_effective_charge_from_beta( const long*  n,
    const float*  beta,
    const long*  Z,
    float*  effective_charge)
{
  // loop over n
  long  i;
  for (i = 0; i < *n; i++){
    // Return effective charge according to Barkas-Bethe-approximation (but not for protons!)
    //    if (Z[i]!=1){
    effective_charge[i]  = (float)(Z[i]) * (1 - (float)exp(-125.0f * beta[i] / (pow(Z[i], 2.0f/3.0f))));//}
    //    else{
    //      effective_charge[i]  = (float)(Z[i]);}
  }
}

void AT_beta_from_particle_no(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    // results
    float*  beta)
{

  //TODO change to another routine (avoid pmatchi), like get_mass_of_particle
  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  float*  mass  =  (float*)calloc(*n, sizeof(float));

  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    mass[i]  = AT_Particle_Data.mass[matches[i]];
  }

  AT_beta_from_mass(  n,
      E_MeV_u,
      mass,
      beta);
  free(mass);
  free(matches);
}


void AT_E_from_beta_and_particle_no( const  long*  n,
    const float*  beta,
    const long*  particle_no,
    float*  E_MeV_u)
{
  //TODO change to another routine (avoid pmatchi), like get_mass_of_particle
  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  float*  mass  =  (float*)calloc(*n, sizeof(float));

  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    mass[i]  = AT_Particle_Data.mass[matches[i]];}

  AT_E_from_beta_and_mass(  n,
      beta,
      mass,
      E_MeV_u);
  free(mass);
  free(matches);
}


void AT_effective_charge_from_particle_no( const  long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    float*  effective_charge)
{
  // get relativistic speeds for all given particles and energies
  float*  beta  =  (float*)calloc(*n, sizeof(float));
  long*  Z    =  (long*)calloc(*n, sizeof(long));

  AT_beta_from_particle_no(  n,
      E_MeV_u,
      particle_no,
      beta);

  //TODO change to another routine (avoid pmatchi), like get_Z_of_particle
  // find look-up indices for Z's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  long i;
  for (i = 0; i < *n; i++){
    Z[i]  =  AT_Particle_Data.Z[matches[i]];
  }

  AT_effective_charge_from_beta(  n,
      beta,
      Z,
      effective_charge);

  free(beta);
  free(Z);
  free(matches);
}



void AT_scaled_energy(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    float*  scaled_energy)
{
  //TODO change to another routine (avoid pmatchi), like get_mass_of_particle and get_A_of_particle
  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    // total kinetic energy
    float  E_MeV    =  E_MeV_u[i] * AT_Particle_Data.A[matches[i]];

    float mass_fraction = AT_Particle_Data.mass[matches[i]] / AT_Particle_Data.mass[0];

    // Return mass-scaled energy
    scaled_energy[i]  =  E_MeV / mass_fraction ;
  }

  free(matches);
}

void AT_E_MeV_u_from_scaled_energy(  const long*  n,
    const float*  scaled_energy,
    const long*  particle_no,
    float*  E_MeV_u)
{
  //TODO change to another routine (avoid pmatchi), like get_mass_of_particle and get_A_of_particle

  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    float mass_fraction = AT_Particle_Data.mass[matches[i]] / AT_Particle_Data.mass[0];

    float  E_MeV  = scaled_energy[i] * mass_fraction;

    // Return energy per nucleon
    E_MeV_u[i]  =  E_MeV / AT_Particle_Data.A[matches[i]];
  }

  free(matches);
}


void AT_max_E_transfer_MeV(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    // results
    float*  max_E_transfer_MeV)
{
  /* if E_MeV_u < 0:    use non-relativistic formula
   * if E_MeV_u > 0:    use relativistic formula
   */

  float* E_MeV_u_copy   =  (float*)calloc(*n, sizeof(float));
  memccpy(E_MeV_u_copy,E_MeV_u,*n,sizeof(float));

  int*  relativistic    =  (int*)calloc(*n, sizeof(int));
  long  i;
  for (i = 0; i < *n; i++){
    if(E_MeV_u[i] >= 0){
      relativistic[i]      = 1;
    }else{
      relativistic[i]      = 0;
      E_MeV_u_copy[i]     *= -1.0f;
    }
  }

  // get relativistic speeds for all given particles and energies
  float*  beta  =  (float*)calloc(*n, sizeof(float));
  AT_beta_from_particle_no(  n,
      E_MeV_u_copy,
      particle_no,
      beta);

  for (i = 0; i < *n; i++){
    if(relativistic[i] == 0){
      max_E_transfer_MeV[i]  =  4.0f * electron_mass_MeV_c2 / proton_mass_MeV_c2* E_MeV_u[i];
    }else{
      max_E_transfer_MeV[i]  =  2.0f * electron_mass_MeV_c2 * beta[i] * beta[i] / (1.0f - beta[i] * beta[i]);
    }
  }

  free(E_MeV_u_copy);
  free(relativistic);
  free(beta);
}


// Get Bohr's energy spread (Wilson, 1947, Phys Rev 71, 385)
void AT_Bohr_Energy_Straggling_g_cm2(  const long*  n,
    const char**  material_name,
    float*  dsE2dz)
{
  long  i;
  double tmp;
  float  electron_density_m3;
  long  n_dummy = 1;
  for (i = 0; i < *n; i++){
    AT_electron_density_m3(  &n_dummy,
        material_name[i],
        &electron_density_m3);
    tmp                  =  e_C * e_C * e_C * e_C * electron_density_m3;
    tmp                 /=  4.0 * pi * e0_F_m * e0_F_m;
    tmp                 /=  MeV_to_J * MeV_to_J * m_to_cm;
    dsE2dz[i]            =  (float)tmp;
  }
}

