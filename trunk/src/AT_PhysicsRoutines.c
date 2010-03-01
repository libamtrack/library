/**
 * @file
 * @brief Physics related routines
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

  const long dummy_n = 1;
  const long proton_particle_no = 1;
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

  const long dummy_n = 1;
  const long proton_particle_no = 1;
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
   * if E_MeV_u < 0:    use non-relativistic formula
   * if E_MeV_u > 0:    use relativistic formula
   * // TODO instead of using negative values of the energy switch parameter "relativistic" should be added to argument list
   */

  float* E_MeV_u_copy   =  (float*)calloc(*n, sizeof(float));
  memcpy(E_MeV_u_copy,E_MeV_u,(*n)*sizeof(float));

  int*  relativistic    =  (int*)calloc(*n, sizeof(int));
  long  i;
  for (i = 0; i < *n; i++){
    if(E_MeV_u_copy[i] >= 0){
      relativistic[i]      = 1;
    }else{
      relativistic[i]      = 0;
      E_MeV_u_copy[i]     *= -1.0f; // E_MeV_u_copy is negative so we set it back to positive value
    }
  }

  // get relativistic speeds for all given particles and energies
  float*  beta  =  (float*)calloc(*n, sizeof(float));
  AT_beta_from_E(  n,
      E_MeV_u_copy,
      beta);

  for (i = 0; i < *n; i++){
    if(relativistic[i] == 0){
      max_E_transfer_MeV[i]  =  4.0f * electron_mass_MeV_c2 / proton_mass_MeV_c2 * E_MeV_u_copy[i];
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

    tmp                  =  gsl_pow_4(e_C) * electron_density_m3;
    tmp                 /=  4.0 * M_PI * gsl_pow_2(e0_F_m);
    tmp                 /=  gsl_pow_2(MeV_to_J) * m_to_cm;
    dsE2dz[i]            =  (float)tmp;
  }
}

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

void AT_fluence_cm2(  const long*  n,
    const float*  E_MeV_u,
    const long* particle_no,
    const float* D_Gy,
    const long* material_no,
    float* fluence_cm2)
{
  // Get LET (write already into fluence_cm2)
  AT_LET_MeV_cm2_g(n,
      E_MeV_u,
      particle_no,
      material_no,
      fluence_cm2);
  // Divide by dose, convert from Gy to MeV/g
  long  i;
  for (i = 0; i < *n; i++){
    fluence_cm2[i]     =       fluence_cm2[i] / (D_Gy[i] * MeV_g_to_J_kg);
  }
}

void AT_interparticleDistance_m( const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  fluence_cm2,
    float*  results_m)
{
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

// TODO shall it be split into separate functions, like FWHM_to_sigma_cm and so on... ?
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
        fluence_cm2[i] = N[i] / (gsl_pow_2(sigma_cm[i]) * 2.0f * M_PI);
      }
    }
  }else{
    for (i = 0; i < *n; i++){
      N[i]           = fluence_cm2[i] * gsl_pow_2(sigma_cm[i]) * 2.0f * M_PI;
    }
  }
}

void AT_inv_interparticleDistance_Gy( const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  distance_m,
    float*  results_Gy)
{
  long i;
  float fluence;
  for( i = 0 ; i < *n ; i++ ){
    fluence = gsl_pow_2(2.0f/distance_m[i]) * M_1_PI * 1e-4;
    results_Gy[i] = fluence * (LET_MeV_cm2_g[i] * MeV_g_to_J_kg);
  }
}

void AT_inv_interparticleDistance_cm2( const long*   n,
    const float*  distance_m,
    float*  results_cm2)
{
  long i;
  for( i = 0 ; i < *n ; i++ ){
    results_cm2[i] = gsl_pow_2(2.0f/distance_m[i]) * M_1_PI * 1e-4;
  }
}

// TODO: Allow for passing negative values for fluence_cm2 meaning dose then?
void AT_single_impact_fluence_cm2( const long* n,
    const float* E_MeV_u,
    const long* material_no,
    const long* er_model,
    float* single_impact_fluence_cm2)
{
  float* max_electron_range_m    = (float*)calloc(*n, sizeof(float));

  // get max. electron ranges
  void AT_max_electron_range_m( n,
      E_MeV_u,
      material_no,
      er_model,
      max_electron_range_m);

  long i;
  for( i = 0 ; i < *n ; i++ ){
    single_impact_fluence_cm2[i] = M_1_PI / gsl_pow_2( max_electron_range_m[i] * m_to_cm ); // pi * r_max_m^2 = Track area -> single_impact_fluence [1/cm2]
   }

  free(max_electron_range_m);
}

void AT_total_D_Gy( const long* n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* total_dose_Gy)
{

  float*  single_doses_Gy        =  (float*)calloc(*n, sizeof(float));

  AT_D_Gy(      n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      single_doses_Gy);

  *total_dose_Gy = 0.0f;
  long i;
  for (i = 0; i < *n; i++){
    *total_dose_Gy       += single_doses_Gy[i];
  }
  free(single_doses_Gy);
}

void AT_total_fluence_cm2( const long* n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* D_Gy,
    const long* material_no,
    float* total_fluence_cm2)
{

  float*  single_fluences_cm2        =  (float*)calloc(*n, sizeof(float));

  AT_fluence_cm2(      n,
      E_MeV_u,
      particle_no,
      D_Gy,
      material_no,
      single_fluences_cm2);

  *total_fluence_cm2 = 0.0f;
  long i;
  for (i = 0; i < *n; i++){
    *total_fluence_cm2       += single_fluences_cm2[i];
  }
  free(single_fluences_cm2);
}

void AT_fluenceweighted_E_MeV_u( const long*     n,
    const float* E_MeV_u,
    const float* fluence_cm2,
    float* average_E_MeV_u)
 {
  long i;

  float total_fluence_cm2 = 0.0f;
  for (i = 0; i < *n; i++){
    total_fluence_cm2 += fluence_cm2[i];
  }

  *average_E_MeV_u      = 0.0f;
  for (i = 0; i < *n; i++){
     *average_E_MeV_u += fluence_cm2[i] * E_MeV_u[i];
   }

   *average_E_MeV_u /= total_fluence_cm2;
 }


void AT_doseweighted_E_MeV_u( const long*     n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
   float* doseweighted_E_MeV_u)
 {
  long i;

  float*  single_doses_Gy        =  (float*)calloc(*n, sizeof(float));

  AT_D_Gy(      n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      single_doses_Gy);

  float total_dose_Gy = 0.0f;

  for (i = 0; i < *n; i++){
    total_dose_Gy       += single_doses_Gy[i];
  }

  *doseweighted_E_MeV_u      = 0.0f;
  for (i = 0; i < *n; i++){
     *doseweighted_E_MeV_u += single_doses_Gy[i] * E_MeV_u[i];
   }

   *doseweighted_E_MeV_u /= total_dose_Gy;

   free(single_doses_Gy);
}

void AT_fluenceweighted_LET_MeV_cm2_g( const long*     n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* average_LET_MeV_cm2_g)
 {
  long i;

  float*  single_LETs_MeV_cm2_g        =  (float*)calloc(*n, sizeof(float));

  float total_fluence_cm2 = 0.0f;
  for (i = 0; i < *n; i++){
    total_fluence_cm2 += fluence_cm2[i];
  }

  AT_LET_MeV_cm2_g(  n,
      E_MeV_u,
      particle_no,
      material_no,
      single_LETs_MeV_cm2_g);

  *average_LET_MeV_cm2_g      = 0.0f;
  for (i = 0; i < *n; i++){
     *average_LET_MeV_cm2_g += fluence_cm2[i] * single_LETs_MeV_cm2_g[i];
   }

   *average_LET_MeV_cm2_g /= total_fluence_cm2;

   free(single_LETs_MeV_cm2_g);
}


void AT_doseweighted_LET_MeV_cm2_g( const long*     n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* doseweighted_LET_MeV_cm2_g)
 {
  long i;

  float*  single_LETs_MeV_cm2_g  =  (float*)calloc(*n, sizeof(float));
  float*  single_doses_Gy        =  (float*)calloc(*n, sizeof(float));

  AT_LET_MeV_cm2_g(  n,
      E_MeV_u,
      particle_no,
      material_no,
      single_LETs_MeV_cm2_g);

  AT_D_Gy(      n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      single_doses_Gy);

  float total_dose_Gy = 0.0f;

  for (i = 0; i < *n; i++){
    total_dose_Gy       += single_doses_Gy[i];
  }

  *doseweighted_LET_MeV_cm2_g      = 0.0f;
  for (i = 0; i < *n; i++){
     *doseweighted_LET_MeV_cm2_g += single_doses_Gy[i] * single_LETs_MeV_cm2_g[i];
   }

   *doseweighted_LET_MeV_cm2_g /= total_dose_Gy;

   free(single_LETs_MeV_cm2_g);
   free(single_doses_Gy);
 }




