#ifndef AT_PHYSICSROUTINES_H_
#define AT_PHYSICSROUTINES_H_

/**
 * @file
 * @brief Physics related routines
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

/**
 *  Return relativistic speed
 *
 * @param[in] n                        number of particles
 * @param[in] E_MeV_u                  vector of energies of particle per nucleon [MeV]
 * @param[out] beta                    vector of relative particle speed beta = v/c
 */
void AT_beta_from_E( const long*  n,
    const float*  E_MeV_u,
    float*  beta);

/**
 *  Return energy per nucleon of particle with relative speed beta
 *
 * @param[in] n                        number of particles
 * @param[in] beta                     vector of relative particle speed beta = v/c
 * @param[out] E_MeV_u                 vector of energies of particle per nucleon [MeV]
 */
void AT_E_from_beta(  const long*  n,
    const float*  beta,
    float*  E_MeV_u);

/**
 * Effective charge
 * @param[in] n                        number of particles
 * @param[in] beta                     vector of relative particle speed beta = v/c
 * @param[in] Z                        TODO
 * @param[out] effective_charge
 */
void AT_effective_charge_from_beta(  const long*  n,
    const float*  beta,
    const long*  Z,
    float*  effective_charge);

/**
 * Get Bohr's energy spread (Wilson, 1947, Phys Rev 71, 385)
 * @param[in] n                        number of particles
 * @param[in] material_no              TODO
 * @param[out] dsE2dz                  TODO
 */
void AT_Bohr_Energy_Straggling_g_cm2(  const long*  n,
    const long*  material_no,
    float*  dsE2dz);

/**
 * Effective charge
 * @param[in] n                        number of particles
 * @param[in] E_MeV_u                  TODO
 * @param[in] particle_no              TODO
 * @param[out] effective_charge
 */
void AT_effective_charge_from_particle_no(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    float*  effective_charge);

/**
 * Scaled energy TODO
 * @param[in] n                        number of particles
 * @param[in] E_MeV_u                  TODO
 * @param[in] particle_no              TODO
 * @param[out] scaled_energy
 */
void AT_scaled_energy( const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    float*  scaled_energy);

/**
 * Scaled energy TODO
 * @param[in] n                        number of particles
 * @param[in] scaled_energy            TODO
 * @param[in] particle_no              TODO
 * @param[out] E_MeV_u
 */
void AT_E_MeV_u_from_scaled_energy(  const long*  n,
    const float*  scaled_energy,
    const long*  particle_no,
    float*  E_MeV_u);

/**
 * Max energy transfer TODO
 * @param[in] n                        number of particles
 * @param[in] E_MeV_u                  TODO
 * @param[out] max_E_transfer_MeV
 */
void AT_max_E_transfer_MeV(  const long*  n,
    const float*  E_MeV_u,
    float*  max_E_transfer_MeV);

/**
 * Returns dose in Gy for each given particle
 * @param[in]  n            number of particle types in the mixed particle field (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  fluence_cm2  fluence for each particle type (pointer to array of size n)
 * @param[in]  particle_no  type of the particles in the mixed particle field (pointer to array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out] D_Gy         pointer to vector of size n to be allocated by the user which will be used to return the results
 */
void AT_D_Gy(  const long*  n,
    const float*  E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* D_Gy);

/**
 * Returns fluence in 1/cm2 for each given particle
 * @param[in]  n            number of particle types in the mixed particle field (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  D_Gy         dose / Gy for each particle type (pointer to array of size n)
 * @param[in]  particle_no  type of the particles in the mixed particle field (pointer to array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out] fluence_cm2         pointer to vector of size n to be allocated by the user which will be used to return the results
 */
void AT_fluence_cm2(  const long*  n,
    const float*  E_MeV_u,
    const long* particle_no,
    const float* D_Gy,
    const long* material_no,
    float* fluence_cm2);

/**
 * Converts pair-wise physical beam parameters of a symmetric, double Gaussian shape beam (lateral), i.e.
 * central (peak) fluence / width (std.dev.)
 * and accelerator parameters, i.e.
 * total number of particle / FWHM
 *
 * The routine completes the missing data, e.g. FWHM if sigma given, fluence_cm2 (if set 0) if N given etc.
 * If both sigma_cm and FWHM or fluence_cm2 and N are given the physical parameters are taken and the acc. reevaluated
 *
 * @param[in]      n             length of vectors for parameters (pointer to single variable)
 * @param[in,out]  fluence_cm2   fluence in beam center (pointer to array of size n)
 * @param[in,out]  sigma_cm      beam width stdev (pointer to array of size n)
 * @param[in,out]  N             pointer to vector of size n to be allocated by the user which will be used to return the results (absolute particle number)
 * @param[in,out]  FWHM_mm       pointer to vector of size n to be allocated by the user which will be used to return the results (in mm)
 */
void AT_convert_beam_parameters(  const long*  n,
    float* fluence_cm2,
    float* sigma_cm,
    float* N,
    float* FWHM_mm);

/**
 * Interparticle distance TODO
 * @param[in]      n               length of vectors for parameters (pointer to single variable)
 * @param[in]      LET_MeV_cm2_g   TODO
 * @param[in]      fluence_cm2     TODO
 * @param[out]     results_m       TODo
 */
void AT_interparticleDistance_m(       const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  fluence_cm2,
    float*  results_m
);

/**
 * Interparticle distance TODO
 * @param[in]      n               length of vectors for parameters (pointer to single variable)
 * @param[in]      LET_MeV_cm2_g   TODO
 * @param[in]      distance_m      TODO
 * @param[out]     results_Gy      TODO
 */
void AT_inv_interparticleDistance_Gy(  const long*   n,
    const float*  LET_MeV_cm2_g,
    const float*  distance_m,
    float*  results_Gy
);

/**
 * Interparticle distance TODO
 * @param[in]      n               length of vectors for parameters (pointer to single variable)
 * @param[in]      distance_m      TODO
 * @param[out]     results_cm2     TODO
 */
void AT_inv_interparticleDistance_cm2( const long*   n,
    const float*  distance_m,
    float*  results_cm2
);

/**
 * Computes the fluences at which (for a given material and electron-range model) every
 * point of the detector lies within the area ONE track only
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  er_model     index of electron-range model
 * @see          AT_ElectronRange.h for definition
 * @param[out] single_impact_fluence_cm2  results (one for each entry in the parameter vectors)
  */
void AT_single_impact_fluence_cm2( const long* n,
    const float* E_MeV_u,
    const long* material_no,
    const long* er_model,
    float* single_impact_fluence_cm2);

/**
 * Computes the total dose of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  particle_no  particle index (pointer to array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out] total_dose_Gy  result (pointer to float)
  */
void AT_total_D_Gy( const long* n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* total_dose_Gy);

/**
 * Computes the total fluence of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  particle_no  particle index (pointer to array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  D_Gy  doses of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out] total_fluence_cm  result (pointer to float)
  */
void AT_total_fluence_cm2( const long* n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* D_Gy,
    const long* material_no,
    float* total_fluence_cm2);

/**
 * Computes the fluence-weighted average energy of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (pointer to array of size n)
 * @param[out] average_E_MeV_u  result (pointer to float)
 */
void AT_fluenceweighted_E_MeV_u( const long*     n,
    const float* E_MeV_u,
    const float* fluence_cm2,
    float* average_E_MeV_u);

/**
 * Computes the dose-weighted average energy of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  particle_no  particle index (pointer to array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out] doseweighted_E_MeV_u  result (pointer to float)
 */
void AT_doseweighted_E_MeV_u( const long*     n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* doseweighted_E_MeV_u);

/**
 * Computes the fluence-weighted average LET of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  particle_no  particle index (pointer to array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out] fluenceweighted_LET_MeV_cm2_g  result (pointer to float)
 */
void AT_fluenceweighted_LET_MeV_cm2_g( const long*     n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* fluenceweighted_LET_MeV_cm2_g);

/**
 * Computes the dose-weighted average LET of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters (pointer to single variable)
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  particle_no  particle index (pointer to array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (pointer to array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (pointer to array of size n)
 * @param[out] doseweighted_LET_MeV_cm2_g  result (pointer to float)
 */
void AT_doseweighted_LET_MeV_cm2_g( const long*     n,
    const float* E_MeV_u,
    const long* particle_no,
    const float* fluence_cm2,
    const long* material_no,
    float* doseweighted_LET_MeV_cm2_g);

#endif /* AT_PHYSICSROUTINES_H_ */
