/**
 * @file
 * @brief Source file for Physics related routines
 */

/*
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

void AT_beta_from_E( const long*  n,
    const float*  E_MeV_u,
    float*  beta)
{
  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){

    // Return relativistic speed
    // E_MeV/E0_MeV = E_MeV_u[i] / proton_mass_MeV_c2
    // beta = sqrt(1. - 1/((1 + E_MeV/E0_MeV)*(1 + E_MeV/E0_MeV)))
    beta[i]        =  (float)sqrt(1.0f - 1.0f/gsl_pow_2(1 + E_MeV_u[i]/proton_mass_MeV_c2));
  }
}

void AT_E_from_beta(  const long*  n,
    const float*  beta,
    float*  E_MeV_u)
{
  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    // Get rest energy
    //float  E0_MeV    =  (float)proton_mass_MeV_c2 * mass[i];
    //E_MeV_u[i]      =  E0_MeV * (1.0f / (1 - beta[i]*beta[i]) - 1);
    //E_MeV_u[i]      /=  mass[i];

    E_MeV_u[i]      =  proton_mass_MeV_c2 * (1.0f / (1.0f - gsl_pow_2(beta[i])) - 1.0f);
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
    if (Z[i]!=1){
      effective_charge[i]  = (float)(Z[i]) * (1 - (float)exp(-125.0f * beta[i] / (pow(Z[i], 2.0f/3.0f))));//}
    }else{
      effective_charge[i]  = (float)(Z[i]);
    }
  }
}

void AT_effective_charge_from_particle_no( const  long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    float*  effective_charge)
{
  // get relativistic speeds for all given particles and energies
  float*  beta  =  (float*)calloc(*n, sizeof(float));
  long*  Z    =  (long*)calloc(*n, sizeof(long));

  AT_beta_from_E(  n,
      E_MeV_u,
      beta);

  AT_Z_from_particle_no(n,particle_no,Z);

  AT_effective_charge_from_beta(  n,
      beta,
      Z,
      effective_charge);

  free(beta);
  free(Z);
}



void AT_scaled_energy(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    float*  scaled_energy)
{
  long*  A    =  (long*)calloc(*n, sizeof(long));
  float*  mass  =  (float*)calloc(*n, sizeof(float));

  AT_Particle_Properties(n,particle_no,NULL,NULL,NULL,NULL,A,mass);

  long dummy_n = 1;
  long proton_particle_no = 1;
  float proton_mass;

  AT_mass_from_particle_no(&dummy_n,&proton_particle_no,&proton_mass);

  // loop over n to find beta for all given particles and energies
  long  i;
  float  E_MeV;
  float mass_fraction;
  for(i = 0; i < *n; i++){
    // total kinetic energy
    E_MeV    =  E_MeV_u[i] * A[i];

    mass_fraction = mass[i] / proton_mass;

    // Return mass-scaled energy
    scaled_energy[i]  =  E_MeV / mass_fraction ;
  }

  free(A);
  free(mass);
}

void AT_E_MeV_u_from_scaled_energy(  const long*  n,
    const float*  scaled_energy,
    const long*  particle_no,
    float*  E_MeV_u)
{
  long*  A    =  (long*)calloc(*n, sizeof(long));
  float*  mass  =  (float*)calloc(*n, sizeof(float));

  AT_Particle_Properties(n,particle_no,NULL,NULL,NULL,NULL,A,mass);

  long dummy_n = 1;
  long proton_particle_no = 1;
  float proton_mass;

  AT_mass_from_particle_no(&dummy_n,&proton_particle_no,&proton_mass);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    float mass_fraction = mass[i] / proton_mass;

    float  E_MeV  = scaled_energy[i] * mass_fraction;

    // Return energy per nucleon
    E_MeV_u[i]  =  E_MeV / A[i];
  }

  free(A);
  free(mass);
}


void AT_max_E_transfer_MeV(  const long*  n,
    const float*  E_MeV_u,
    // results
    float*  max_E_transfer_MeV)
{
  /**
   *  if E_MeV_u < 0:    use non-relativistic formula
   * if E_MeV_u > 0:    use relativistic formula
   */

  float* E_MeV_u_copy   =  (float*)calloc(*n, sizeof(float));
  memcpy(E_MeV_u_copy,E_MeV_u,(*n)*sizeof(float));

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
  AT_beta_from_E(  n,
      E_MeV_u_copy,
      beta);

  for (i = 0; i < *n; i++){
    if(relativistic[i] == 0){
      max_E_transfer_MeV[i]  =  4.0f * electron_mass_MeV_c2 / proton_mass_MeV_c2 * E_MeV_u[i];
    }else{
      max_E_transfer_MeV[i]  =  2.0f * electron_mass_MeV_c2 * gsl_pow_2(beta[i]) / (1.0f - gsl_pow_2(beta[i]));
    }
  }

  free(E_MeV_u_copy);
  free(relativistic);
  free(beta);
}


void AT_Bohr_Energy_Straggling_g_cm2(  const long*  n,
    const long*  material_no,
    float*  dsE2dz)
{
  long  i;
  double tmp;
  float  electron_density_m3;
  long  n_dummy = 1;
  for (i = 0; i < *n; i++){

    AT_getMaterialData(  &n_dummy,
        material_no,
        NULL,
        &electron_density_m3,
        NULL, NULL, NULL, NULL, NULL, NULL);

    tmp                  =  e_C * e_C * e_C * e_C * electron_density_m3;
    tmp                 /=  4.0 * pi * e0_F_m * e0_F_m;
    tmp                 /=  MeV_to_J * MeV_to_J * m_to_cm;
    dsE2dz[i]            =  (float)tmp;
  }
}

// TODO: Use this routine in AT_SC_get_f1... rather than local particle field to dose conversions
void AT_D_Gy(  const long*  n,
    const float*        E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* D_Gy)
{
  // Get LET (write already into D_Gy)
  AT_LET_MeV_cm2_g(n,
      E_MeV_u,
      particle_no,
      material_no,
      D_Gy);
  // Multiply by fluence, convert from MeV/g to Gy
  long  i;
  for (i = 0; i < *n; i++){
    D_Gy[i]     =       D_Gy[i] * fluence_cm2[i] * MeV_g_to_J_kg;
  }
}

void AT_interparticleDistance_m( const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  fluence_cm2,
    float*  results_m
){
  long i;
  float fluence;
  for( i = 0 ; i < *n ; i++ ){
    if( fluence_cm2[i] > 0 ){
      results_m[i] = 2.0f / sqrt(M_PI*1e4*fluence_cm2[i]);
    } else {
      fluence = (-fluence_cm2[i]) / (LET_MeV_cm2_g[i] * MeV_g_to_J_kg);
      results_m[i] = 2.0f / sqrt(M_PI*1e4*fluence);
    }
  }
}

void AT_convert_beam_parameters(  const long*  n,
    float* fluence_cm2,
    float* sigma_cm,
    float* N,
    float* FWHM_mm)
{
  long  i;
  if(sigma_cm[0] == 0.0f){
    for (i = 0; i < *n; i++){
      sigma_cm[i]    = FWHM_mm[i] / (2.354820046f * cm_to_mm);                                // 2 * sqrt(2*ln(2))
    }
  }else{
    for (i = 0; i < *n; i++){
      FWHM_mm[i]     = sigma_cm[i] * (2.354820046f * cm_to_mm);
    }
  }

  if(fluence_cm2[0] == 0.0f){
    for (i = 0; i < *n; i++){
      if(sigma_cm[i] != 0.0f){
        fluence_cm2[i] = N[i] / (sigma_cm[i] * sigma_cm[i] * 2.0f * pi);
      }
    }
  }else{
    for (i = 0; i < *n; i++){
      N[i]           = fluence_cm2[i] * sigma_cm[i] * sigma_cm[i] * 2.0f * pi;
    }
  }
}

void AT_inv_interparticleDistance_Gy( const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  distance_m,
    float*  results_Gy
){
  long i;
  float fluence;
  for( i = 0 ; i < *n ; i++ ){
    fluence = (2.0f/distance_m[i])*(2.0f/distance_m[i])*M_1_PI*1e-4;
    results_Gy[i] = fluence * (LET_MeV_cm2_g[i] * MeV_g_to_J_kg);
  }
}

void AT_inv_interparticleDistance_cm2( const long*   n,
    const float*  distance_m,
    float*  results_cm2
){
  long i;
  for( i = 0 ; i < *n ; i++ ){
    results_cm2[i] = (2.0f/distance_m[i])*(2.0f/distance_m[i])*M_1_PI*1e-4;
  }
}
