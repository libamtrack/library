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
 * Structure to carry essential single, monoenergetic particle field information
 */
typedef struct {
  double   E_MeV_u;                /** energy of the particles in MeV/u */
  double   particle_no;            /** particle index */
  double   fluence_cm2;            /** fluence_cm2 */
} single_field_data;


/**
 * Structure to carry essential detector information
 */
typedef struct {
  double   material_no;            /** material index */
} detector_data;


/**
 * Structure to carry derived single, monoenergetic particle field information, e.g. needed for array building
 * in AT_GSM, AT_SPIFF, ...
 *
 * This structure replaces the f1_parameter array
 */
typedef struct {
  double   LET_MeV_cm2_g;                /** LET (in MeV*cm2/g) of the particle field */
  double   dEdx_MeV_cm2_g;               /** Energy loss dE/dx (in MeV*cm2/g) of the particle field, depending on the chosen RDD, this does not have to meet the LET (!) */
  double   d_min_Gy;                     /** lowest local dose found in the field (in Gy) */
  double   d_max_Gy;                     /** highest local dose found in the field (not considering track overlap which can cause core overlap and thus even higher doses, in Gy) */
  double   r_min_m;                      /** radius closest to the track center in the field (in m) */
  double   r_max_m;                      /** widest radius from track center in the field (in m) */
  double   single_impact_fluence_cm2;    /** fluences at which every point of the detector lies within the area ONE track only */
  double   single_impact_dose_Gy;        /** corresponding dose */
  double   normalization;                /** normalization constant for RDD (e.g. to meet LET) */
} single_field_information;


/**
 * Structure to carry derived mixed particle field information, e.g. needed for array building
 * in AT_GSM, AT_SPIFF, ...
 *
 * This structure replaces the f_parameter array
 */
typedef struct {
  double   total_fluence_cm2;                /** total fluence of all particles the field */
  double   total_dose_Gy;                    /** total dose delivered by the field */
  double   fluenceweighted_E_MeV_u;          /** fluence-weighted average energy (in MeV/u) */
  double   doseweighted_E_MeV_u;             /** dose-weighted average energy (in MeV/u) */
  double   fluenceweighted_LET_MeV_cm2_g;    /** fluence-weighted average LET (in MeV*cm2/g) */
  double   doseweighted_LET_MeV_cm2_g;       /** dose-weighted average LET (in MeV*cm2/g) */
  double   u;                                /** average number of track contributing to a detector voxel, needed by AT_SPIFF, AT_SPISS */
} mixed_field_information;


/**
 *  Returns relativistic speed for single value of energy
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @return     beta                     relative particle speed beta = v/c
 */
inline double AT_beta_from_E_single( const double  E_MeV_u );


/**
 *  Returns relativistic speed for many particles
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV]
 * @param[out] beta                     vector of relative particle speed beta = v/c
 * @return     status code
 */
int AT_beta_from_E( const long  n,
    const double  E_MeV_u[],
    double  beta[]);


/**
 *  Return energy per nucleon of particle with relative speed beta
 *
 * @param[in]  beta                     relative particle speed beta = v/c
 * @return                              energy of particle per nucleon [MeV]
 */
inline double AT_E_from_beta_single(  const double beta );


/**
 *  Return energy per nucleon of particle with relative speed beta
 *
 * @param[in]  n                        number of particles
 * @param[in]  beta                     vector of relative particle speed beta = v/c
 * @param[out] E_MeV_u                  vector of energies of particle per nucleon [MeV]
 * @return     status code
 */
int AT_E_from_beta(  const long  n,
    const double  beta[],
    double  E_MeV_u[]);


/**
 * Effective charge according to Barkas-Bethe-approximation:
 *
 * Zeff = Z * exp( -125 * beta / Z^(2/3) )
 *
 * calculated for particle with given relative speed beta
 *
 * @param[in]  beta                     relative particle speed beta = v/c
 * @param[in]  Z                        atomic number Z of ion
 * @return     effective_charge of ion
 */
inline double AT_effective_charge_from_beta_single(  const double beta,
    const long Z);


/**
 * Effective charge according to Barkas-Bethe-approximation:
 *
 * Zeff = Z * exp( -125 * beta / Z^(2/3) )
 *
 * calculated for particle with given relative speed beta
 *
 * @param[in]  n                        number of particles
 * @param[in]  beta                     vector of relative particle speed beta = v/c
 * @param[in]  Z                        atomic number Z of ion
 * @param[out] effective_charge of ion
 * @return     status code
 */
int AT_effective_charge_from_beta(  const long  n,
    const double  beta[],
    const long    Z[],
    double        effective_charge[]);


/**
 * Get Bohr's energy spread (Wilson, 1947, Phys Rev 71, 385)
 * @param[in]  n                        number of particles
 * @param[in]  material_no              index number for detector material
 * @param[out] dsE2dz                   Increase of energy spread (variance s_E^2) with thickness of material: d(s_E^2)/dz
 */
void AT_Bohr_Energy_Straggling_g_cm2(  const long*  n,
    const long*  material_no,
    double*  dsE2dz);


/**
 * Effective charge according to Barkas-Bethe-approximation:
 *
 * Zeff = Z * exp( -125 * beta / Z^(2/3) )
 *
 * calculated for particle with given energy per nucleon
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_no              type of the particles in the mixed particle field
 * @return     effective_charge         Effective charge according to Barkas-Bethe-approximation
 */
double AT_effective_charge_from_E_MeV_u_single(  const double E_MeV_u,
    const long  particle_no);


/**
 * Effective charge according to Barkas-Bethe-approximation:
 *
 * Zeff = Z * exp( -125 * beta / Z^(2/3) )
 *
 * calculated for particle with given energy per nucleon
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV]
 * @param[in]  particle_no              type of the particles in the mixed particle field (array of size n)
 * @param[out] effective_charge         Effective charge according to Barkas-Bethe-approximation
 * @return     status code
 */
int AT_effective_charge_from_E_MeV_u(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    double        effective_charge[]);


/**
 * Max relativistic energy transfer for single particle TODO
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @return max_E_transfer_MeV
 */
inline double AT_max_relativistic_E_transfer_MeV_single( const double E_MeV_u );


/**
 * Max classic energy transfer for single particle TODO
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @return max_E_transfer_MeV
 */
inline double AT_max_classic_E_transfer_MeV_single( const double E_MeV_u );


/**
 * Max energy transfer for single particle TODO
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV]
 * @param[out] max_E_transfer_MeV
 */
inline double AT_max_E_transfer_MeV_single( const double E_MeV_u);


/**
 * Max energy transfer TODO
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[out] max_E_transfer_MeV       (array of size n)
 * @return     status code
 */
int AT_max_E_transfer_MeV(  const long  n,
    const double  E_MeV_u[],
    double        max_E_transfer_MeV[]);


/**
 * Returns dose in Gy for each given particle
 * @param[in]  n            number of particle types in the mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  fluence_cm2  fluence for each particle type (array of size n)
 * @param[in]  particle_no  type of the particles in the mixed particle field (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out] D_Gy         vector of size n to be allocated by the user which will be used to return the results
 */
void AT_D_Gy(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    double        D_Gy[]);


/**
 * Returns fluence in 1/cm2 for each given particle
 * @param[in]  n            number of particle types in the mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  D_Gy         dose / Gy for each particle type (array of size n)
 * @param[in]  particle_no  type of the particles in the mixed particle field (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out] fluence_cm2         vector of size n to be allocated by the user which will be used to return the results
 */
void AT_fluence_cm2(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  D_Gy[],
    const long    material_no,
    double        fluence_cm2[]);


/**
 * Converts pair-wise physical beam parameters of a symmetric, double Gaussian shape beam (lateral), i.e.
 * central (peak) fluence / width (std.dev.)
 * and accelerator parameters, i.e.
 * total number of particle / FWHM
 *
 * The routine completes the missing data, e.g. FWHM if sigma given, fluence_cm2 (if set 0) if N given etc.
 * If both sigma_cm and FWHM or fluence_cm2 and N are given the physical parameters are taken and the acc. reevaluated
 *
 * @param[in]      n             length of vectors for parameters
 * @param[in,out]  fluence_cm2   fluence in beam center (array of size n)
 * @param[in,out]  sigma_cm      beam width stdev (array of size n)
 * @param[in,out]  N             vector of size n to be allocated by the user which will be used to return the results (absolute particle number)
 * @param[in,out]  FWHM_mm       vector of size n to be allocated by the user which will be used to return the results (in mm)
 */
void AT_convert_beam_parameters(  const long  n,
    double fluence_cm2[],
    double sigma_cm[],
    double N[],
    double FWHM_mm[]);


/**
 * Interparticle distance TODO
 * @param[in]      n               length of vectors for parameters
 * @param[in]      LET_MeV_cm2_g   LET for each particle type (array of size n)
 * @param[in]      fluence_cm2     fluence for each particle type (array of size n)
 * @param[out]     results_m       interparticle distance for each particle type (array of size n)
 */
void AT_interparticleDistance_m(       const long   n,
    const double  LET_MeV_cm2_g[],
    const double  fluence_cm2[],
    double        results_m[]
);


/**
 * Inverse interparticle distance TODO
 * @param[in]      n               length of vectors for parameters
 * @param[in]      LET_MeV_cm2_g   LET for each particle type (array of size n)
 * @param[in]      distance_m      interparticle distance for each particle type (array of size n)
 * @param[out]     results_Gy      TODO
 */
void AT_inv_interparticleDistance_Gy(  const long   n,
    const double   LET_MeV_cm2_g[],
    const double   distance_m[],
    double         results_Gy[]
);


/**
 * Computes the fluences at which (for a given material and electron-range model) every
 * point of the detector lies within the area ONE track only
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  E_MeV_u      energy of particle
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  er_model     index of electron-range model
 * @see          AT_ElectronRange.h for definition
 * @return     single_impact_fluence_cm2  results (one for each entry in the parameter vectors)
  */
double AT_single_impact_fluence_cm2_single( const double E_MeV_u,
    const long material_no,
    const long er_model);


/**
 * Computes the fluences at which (for a given material and electron-range model) every
 * point of the detector lies within the area ONE track only
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  er_model     index of electron-range model
 * @see          AT_ElectronRange.h for definition
 * @param[out] single_impact_fluence_cm2  results (one for each entry in the parameter vectors)
  */
void AT_single_impact_fluence_cm2( const long n,
    const double  E_MeV_u[],
    const long    material_no,
    const long    er_model,
    double        single_impact_fluence_cm2[]);


/**
 * TODO
 * @param LET_MeV_cm2_g
 * @param single_impact_fluence_cm2
 * @return
 */
inline double AT_single_impact_dose_Gy_single( const double LET_MeV_cm2_g,
    const double single_impact_fluence_cm2);


/**
 * Computes the total dose of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @return     total_dose_Gy  result
  */
double  AT_total_D_Gy( const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no);


/**
 * Computes the total fluence of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  D_Gy  doses of particles in the mixed particle field (array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @return       total_fluence_cm  result
  */
double AT_total_fluence_cm2( const long n,
    const double   E_MeV_u[],
    const long     particle_no[],
    const double   D_Gy[],
    const long     material_no);


/**
 * Computes the fluence-weighted average energy of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size n)
 * @return     average_E_MeV_u  result
 */
double AT_fluenceweighted_E_MeV_u( const long    n,
    const double E_MeV_u[],
    const double fluence_cm2[]);


/**
 * Computes the dose-weighted average energy of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @return     doseweighted_E_MeV_u  result
 */
double AT_doseweighted_E_MeV_u( const long   n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no);


/**
 * Computes the fluence-weighted average LET of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @return     fluenceweighted_LET_MeV_cm2_g  result
 */
double AT_fluenceweighted_LET_MeV_cm2_g( const long     n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no);


/**
 * Computes the dose-weighted average LET of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  n            length of vectors for parameters
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size n)
 * @return     doseweighted_LET_MeV_cm2_g  result
 */
double AT_doseweighted_LET_MeV_cm2_g( const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no);

#endif /* AT_PHYSICSROUTINES_H_ */
