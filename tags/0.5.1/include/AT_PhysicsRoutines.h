#ifndef AT_PHYSICSROUTINES_H_
#define AT_PHYSICSROUTINES_H_

/**
 * @brief Physics related routines
 */

/*
 *    AT_PhysicsRoutines.h
 *    ==================
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

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "AT_Constants.h"
#include "AT_DataMaterial.h"
#include "AT_DataParticle.h"
#include "AT_DataStoppingPower.h"
#include "AT_ElectronRange.h"
#include "AT_NumericalRoutines.h"


/**
 * Structure to carry essential single, monoenergetic particle field information
 */
typedef struct {
  double   E_MeV_u;                /** energy of the particles in MeV/u */
  double   particle_no;            /** particle index */
  double   fluence_cm2;            /** fluence_cm2 */
} AT_single_field_data;


/**
 * Structure to carry essential detector information
 */
typedef struct {
  double   material_no;            /** material index */
} AT_detector_data;


/**
 * Structure to carry derived single, monoenergetic particle field information, e.g. needed for array building
 * in AT_GSM, AT_CPPSC, ...
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
} AT_single_field_information;


/**
 * Structure to carry derived mixed particle field information, e.g. needed for array building
 * in AT_GSM, AT_CPPSC, ...
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
  double   u;                                /** average number of track contributing to a detector voxel, needed by AT_CPPSC, AT_SPISS */
} AT_mixed_field_information;


/**
 *  Returns relativistic speed for single value of energy
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @return     beta                     relative particle speed beta = v/c
 */
 double AT_beta_from_E_single( const double  E_MeV_u );


/**
 *  Returns relativistic speed for many particles
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[out] beta                     vector of relative particle speed beta = v/c (array of size n)
 * @return     status code
 */
int AT_beta_from_E( const long  n,
    const double  E_MeV_u[],
    double  beta[]);


/**
 *  Returns energy per nucleon of particle with relative speed beta
 *
 * @param[in]  beta                     relative particle speed beta = v/c
 * @return                              energy of particle per nucleon [MeV]
 */
 double AT_E_from_beta_single(  const double beta );


/**
 *  Returns energy per nucleon of particle with relativistic speed beta
 *
 * @param[in]  n                        number of particles
 * @param[in]  beta                     vector of relative particle speed beta = v/c (array of size n)
 * @param[out] E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @return     status code
 */
int AT_E_from_beta(  const long  n,
    const double  beta[],
    double  E_MeV_u[]);


/**
 *  Returns energy for single value of relativistic gamma
 *
 * @param[in]  gamma
 * @return     E_MeV_u                  energy of particle per nucleon [MeV]
 */
 double AT_E_from_gamma_single( const double gamma );

/**
 *  Returns energy from relativistic gamma
 *
 * @param[in]  n                        number of particles
 * @param[in]  gamma                    vector of results (array of size n)
 * @param[out] E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @return     status code
 */
int AT_E_from_gamma( const long  n,
    const double  gamma[],
    double        E_MeV_u[]);



/**
 *  Returns energy per nucleon of particle with given momentum per nucleon
 *
 * @param[in]  momentum_MeV_c_u         momentum per particle [MeV/c]
 * @return                              energy of particle per nucleon [MeV]
 */
 double AT_E_MeV_u_from_momentum_single( 	const double momentum_MeV_c_u);

/**
 *  Returns energy per nucleon for particles with given momentum per nucleon
 *
 * @param[in]  n                        number of particles
 * @param[in]  momentum_MeV_c_u         vector of particle momenta per nucleon [MeV/c], (array of size n)
 * @param[out] E_MeV_u                  vector of energies of particle per nucleon [MeV], (array of size n)
 * @return     status code
 */
int AT_E_MeV_u_from_momentum_MeV_c_u(  const long  n,
    const double  momentum_MeV_c_u[],
    double        E_MeV_u[]);

/**
 *  Returns relativistic gamma for single value of energy
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @return     gamma
 */
 double AT_gamma_from_E_single( const double E_MeV_u );

/**
 *  Returns relativistic gamma
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[out] gamma                    vector of results (array of size n)
 * @return     status code
 */
int AT_gamma_from_E( const long  n,
    const double  E_MeV_u[],
    double        gamma[]);

/**
 * Effective charge according to Barkas-Bethe-approximation:
 *
 * Zeff = Z *[1- exp( -125 * beta / Z^(2/3) )]
 *
 * calculated for particle with given relative speed beta
 *
 * @param[in]  beta                     relative particle speed beta = v/c
 * @param[in]  Z                        atomic number Z of ion
 * @return     effective_charge of ion
 */
 double AT_effective_charge_from_beta_single(  const double beta,
    const long Z);


/**
 * Effective charge according to Barkas-Bethe-approximation:
 *
 * Zeff = Z *[1-exp( -125 * beta / Z^(2/3) )]
 *
 * calculated for particle with given relative speed beta
 *
 * @param[in]  n                        number of particles
 * @param[in]  beta                     vector of relative particle speed beta = v/c (array of size n)
 * @param[in]  Z                        atomic number Z of ion (array of size n)
 * @param[out] effective_charge of ion
 * @return     status code
 */
int AT_effective_charge_from_beta(  const long  n,
    const double  beta[],
    const long    Z[],
    double        effective_charge[]);


/**
 * Get energy spread with depth according to Bohr's classical theory
 * Bohr, N. (1915), Phil. Mag. 30, 581ff, see also Evans, R.D. (1955), The atomic nucleus, McGraw Hill, New York, p. 661
 * In the literature dsE2dz is often given in units ergs2/cm. Here we report it mass-normalized MeV2*cm2/g
 * Since the effective charge of the particle enters the equation, particle types and energies have to be given
 * The equation is however limited to energies > 10 MeV/u and not too heavy ions
 * TODO: add William extension for relativistic effects (Williams, E.J. (1945), Revs. Mod. Phys. 17, 217ff)
 * @param[in]  n                        number particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_no              type of the particles in the mixed particle field (array of size n)
 * @param[in]  material_no              index number for slab material
 * @param[out] dsE2dz_MeV2_cm2_g        Increase of energy straggling variance sigma_E^2 per unit length of material (array of size n)
 */
void AT_energy_straggling_MeV2_cm2_g(  const long  n,
	const double	E_MeV_u[],
	const long	particle_no[],
    const long  material_no,
    double	dsE2dz_MeV2_cm2_g[]);

/**
 * Get energy spread of an ion beam after traversing
 * a material slab according to Bohr's classical theory.
 * Bohr, N. (1915), Phil. Mag. 30, 581ff, see also Evans, R.D. (1955), The atomic nucleus, McGraw Hill, New York, p. 661
 * Please note that the effective charge is assumed to be constant over the material slab
 * If this is not the case you should apply this routine multiple times to subslices
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_no              type of the particles in the mixed particle field (array of size n)
 * @param[in]  material_no              index number for slab material
 * @param[in]  slab_thickness_m         thickness of slab in m
 * @param[in]  initial_sigma_E_MeV_u    energy spread - 1 sigma - before traversing the slab - can be 0 (array of size n)
 * @param[out] sigma_E_MeV_u            energy spread - 1 sigma - after traversing the slab (array of size n)
 */
void AT_energy_straggling_after_slab_E_MeV_u( const long  n,
	const double	E_MeV_u[],
	const long	particle_no[],
    const long	material_no,
    const double	slab_thickness_m,
    const double	initial_sigma_E_MeV_u[],
    double	sigma_E_MeV_u[]);

/**
 * Effective charge according to Barkas-Bethe-approximation
 * for particle with given energy per nucleon
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_no              type of the particles in the mixed particle field
 * @return     effective_charge         Effective charge according to Barkas-Bethe-approximation
 */
double AT_effective_charge_from_E_MeV_u_single(  const double E_MeV_u,
    const long  particle_no);


/**
 * Effective charge according to Barkas-Bethe-approximation:
 * for particles with given kinetic energy per nucleon
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_no              type of the particles in the mixed particle field (array of size n)
 * @param[out] effective_charge         Effective charge according to Barkas-Bethe-approximation (array of size n)
 * @return     status code
 */
int AT_effective_charge_from_E_MeV_u(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    double        effective_charge[]);


/**
 * Max relativistic energy transfer for single particle
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV/u]
 * @return max_E_transfer_MeV
 */
 double AT_max_relativistic_E_transfer_MeV_single( const double E_MeV_u );


/**
 * Max classic energy transfer for single particle
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV/u]
 * @return max_E_transfer_MeV
 */
 double AT_max_classic_E_transfer_MeV_single( const double E_MeV_u );


/**
 * Max energy transfer for single particle
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV/u]
 * @return     max_E_transfer_MeV
 */
 double AT_max_E_transfer_MeV_single( const double E_MeV_u);


/**
 * Kinetic energy maximally transferred from an ion to an electron
 * in a collision - relativistic or non-relativistic
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  energies of particle per nucleon [MeV/u]; if positive, the computation will be relativistic; if negative, the classic formular will be used (array of size n)
 * @param[out] max_E_transfer_MeV       maximal energies transferred (array of size n)
 * @return     status code
 */
int AT_max_E_transfer_MeV(  const long  n,
    const double  E_MeV_u[],
    double        max_E_transfer_MeV[]);

/**
 *  Returns relativistic momentum (per nucleon) of particle
 *
 * @param	  	E_MeV_u                 kinetic Energy per nucleon
 * @return                              momentum [MeV/c]
 */
 double AT_momentum_from_E_MeV_c_u_single( const double E_MeV_u);

/**
 *  Returns relativistic momenta per nucleon for particles with given kinetic energy
 *
 * @param[in]	n						number of particles
 * @param[in]  	E_MeV_u                 kinetic energy per nucleon (array of size n)
 * @param[out]	momentum_MeV_c  		momentum per nucleon (array of size n)
 * @return                              return code
 */
int AT_momentum_MeV_c_u_from_E_MeV_u( const long  n,
    const double  E_MeV_u[],
    double        momentum_MeV_c[]);

/**
 * Returns dose in Gy for particle with given fluence and energy
 * @param[in]  E_MeV_u      energy per unit mass
 * @param[in]  fluence_cm2  fluence in 1/cm2
 * @param[in]  particle_no  type of the particle
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no  stopping power source index
 * @return     D_Gy         dose in Gy
 */
double AT_dose_Gy_from_fluence_cm2_single(  const double  E_MeV_u,
    const long    particle_no,
    const double  fluence_cm2,
    const long    material_no,
    const long    stopping_power_source_no);


/**
 * Returns dose in Gy for each given particle
 * @param[in]  n            number of particle types in the mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  fluence_cm2  fluence for each particle type (array of size n)
 * @param[in]  particle_no  type of the particles in the mixed particle field (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no  stopping power source index
 * @param[out] dose_Gy          be allocated by the user which will be used to return the results (array of size n)
 */
void AT_dose_Gy_from_fluence_cm2(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    const long    stopping_power_source_no,
    double        dose_Gy[]);


/**
 * Returns fluence in 1/cm2 for particles with given dose and energy
 * @param[in]  E_MeV_u      energy of particle
 * @param[in]  particle_no  type of the particles
 * @see          AT_DataParticle.h for definition
 * @param[in]  D_Gy         dose in Gy
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no  TODO
 * @return fluence in 1/cm2
 */
double AT_fluence_cm2_from_dose_Gy_single( const double  E_MeV_u,
    const long    particle_no,
    const double  D_Gy,
    const long    material_no,
    const long    stopping_power_source_no);


/**
 * Returns fluence in 1/cm2 for each given particle
 * @param[in]  n            number of particle types in the mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no  type of the particles in the mixed particle field (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  D_Gy         dose / Gy for each particle type (array of size n)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no  TODO
 * @param[out] fluence_cm2         to be allocated by the user which will be used to return the results (array of size n)
 */
void AT_fluence_cm2_from_dose_Gy(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  D_Gy[],
    const long    material_no,
    const long    stopping_power_source_no,
    double        fluence_cm2[]);


/**
 * Converts physical beam parameters of a symmetric, double lateral Gaussian shape beam, i.e.
 * central (=peak) fluence and width (= 1 standard deviation)
 * to technical, accelerator parameters, i.e.
 * total number of particles and FWHM
 *
 * @param[in]      n             length of vectors for parameters
 * @param[in]      fluence_cm2   fluence in beam center (array of size n)
 * @param[in]      sigma_cm      beam width stdev (array of size n)
 * @param[out]     N             resulting absolute particle numbers (array of size n)
 * @param[out]     FWHM_mm       resulting FWHMs (in mm) (array of size n)
 */
void AT_beam_par_physical_to_technical(  const long  n,
    const double fluence_cm2[],
    const double sigma_cm[],
    double N[],
    double FWHM_mm[]);

/**
 * Converts technical, accelerator parameters of a symmetric, double lateral Gaussian shape beam, i.e.
 * total number of particles and FWHM to
 * physical beam parameters, i.e.
 * central (=peak) fluence and width (= 1 standard deviation)
 *
 * @param[in]      n             length of vectors for parameters
 * @param[in]      N             absolute particle numbers (array of size n)
 * @param[in]      FWHM_mm       FWHMs (in mm) (array of size n)
 * @param[out]     fluence_cm2   resulting fluence in beam center (array of size n)
 * @param[out]     sigma_cm      resulting beam width stdev (array of size n)
 */
void AT_beam_par_technical_to_physical(  const long  n,
    const double N[],
    const double FWHM_mm[],
    double fluence_cm2[],
    double sigma_cm[]);

/**
 * Interparticle distance
 *
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
 * Inverse interparticle distance
 *
 * @param[in]      n               length of vectors for parameters
 * @param[in]      LET_MeV_cm2_g   LET for each particle type (array of size n)
 * @param[in]      distance_m      interparticle distance for each particle type (array of size n)
 * @param[out]     results_Gy
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
 * @param[out] single_impact_fluence_cm2  results (one for each entry in the parameter vectors) (array of size n)
  */
void AT_single_impact_fluence_cm2( const long n,
    const double  E_MeV_u[],
    const long    material_no,
    const long    er_model,
    double        single_impact_fluence_cm2[]);


/**
 * Dose for the fluence at a single impact
 *
 * @param[in] LET_MeV_cm2_g              LET of particle
 * @param[in] single_impact_fluence_cm2  the fluence corresponing to a single impact
 * @return    single impact dose
 */
double AT_single_impact_dose_Gy_single( const double LET_MeV_cm2_g,
    const double single_impact_fluence_cm2);

/**
 * Doses for the fluences at a single impact
 *
 * @param[in]  n                          number of particles
 * @param[in]  E_MeV_i                    Energy
 * @param[in]  particle_no                particle type
 * @param[in]  material_no                material
 * @param[in]  er_model                   electron-range model
 * @param[in]  stopping_power_source_no   TODO
 * @param[out] single_impact_dose_Gy      resulting single impact doses
 */
void AT_single_impact_dose_Gy( const long n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    const long    er_model,
    const long    stopping_power_source_no,
    double        single_impact_dose_Gy[]);

/**
 * Computes the total dose of a mixed particle field
 *
 * @param[in]  number_of_field_components            number of components in the mixed field
 * @param[in]  E_MeV_u                               energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no                           particle index (array of size number_of_field_components)
 * @see AT_DataParticle.h for definition
 * @param[in]  fluence_cm2                           fluences of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  material_no                           material index
 * @see AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no              TODO
 * @return     total_dose_Gy                         result
  */
double AT_total_D_Gy( const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    const long    stopping_power_source_no);


/**
 * Computes the total fluence of a mixed particle field
 *
 * @param[in]  number_of_field_components            number of components in the mixed field
 * @param[in]  E_MeV_u                               energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no                           particle index (array of size number_of_field_components)
 * @see AT_DataParticle.h for definition
 * @param[in]  D_Gy                                  doses of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  material_no                           material index
 * @param[in]  stopping_power_source_no              TODO
 * @see AT_DataMaterial.h for definition
 * @return     total_fluence_cm                      result
  */
double AT_total_fluence_cm2( const long number_of_field_components,
    const double   E_MeV_u[],
    const long     particle_no[],
    const double   D_Gy[],
    const long     material_no,
    const long     stopping_power_source_no);


/**
 * Computes the fluence-weighted average energy of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  number_of_field_components            number of components in mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size number_of_field_components)
 * @return     average_E_MeV_u  fluence-weighted mean energy
 */
double AT_fluence_weighted_E_MeV_u( const long    number_of_field_components,
    const double E_MeV_u[],
    const double fluence_cm2[]);


/**
 * Computes the dose-weighted average energy of a particle field
 *
 * Needed by SuccessiveConvolutions
 *
 * @param[in]  number_of_field_components            number of components in mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no  particle index (array of size number_of_field_components)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no  TODO
 * @return     dose-weighted mean energy
 */
double AT_dose_weighted_E_MeV_u( const long   number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    const long    stopping_power_source_no);


/**
 * Computes the fluence-weighted average LET of a particle field
 *
 * @param[in]  number_of_field_components            number of components in mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no  particle index (array of size number_of_field_components)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no  TODO
 * @return     fluence-weighted LET
 */
double AT_fluence_weighted_LET_MeV_cm2_g( const long     number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    const long    stopping_power_source_no);


/**
 * Computes the dose-weighted average LET of a particle field
 *
 * @param[in]  number_of_field_components            number of components in mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no  particle index (array of size number_of_field_components)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no  TODO
 * @return     dose-weighted LET
 */
double AT_dose_weighted_LET_MeV_cm2_g( const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    const long    stopping_power_source_no);


/**
 * Computes the stopping power ratio for a material and a reference material.
 *
 * In case of mixed particle fields, the stopping power ratios of individual components are
 * weighted by their respective fluences. Thus, this routines computes the ration of fluence-weighted
 * stopping powers, NOT of dose-weighted stopping powers.
 *
 * @param[in]  number_of_field_components            number of components in mixed field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no  particle index (array of size number_of_field_components)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  reference_material_no  material index of reference material
 * @param[in]  stopping_power_source_no  TODO
 * @return     stopping power ratio
 */
double AT_stopping_power_ratio( const long     number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    const long	  reference_material_no,
    const long    stopping_power_source_no);

/**
 * Computes the number of track contributing to a representative point in a mixed field
 *
 * @param[in]  number_of_field_components            number of components in mixed particle field
 * @param[in]  E_MeV_u      energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no  particle index (array of size number_of_field_components)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2  fluences of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  er_model     chosen electron-range-model
 * @param[in]  stopping_power_source_no     TODO
 * @return     resulting mean number of tracks contributing
 */
double AT_mean_number_of_tracks_contrib(    const long number_of_field_components,
                const double E_MeV_u[],
                const long particle_no[],
                const double fluence_cm2[],
                const long material_no,
                const long er_model,
                const long stopping_power_source_no);

/**
 * Computes the kinetic variable needed for computation of
 * density effect in Bethe formula for stopping power
 * following the Sternheimer (1971) approach
 *
 * @param[in]  E_MeV_u      energy of particle
 * @return     				kinetic variable
 */
double AT_kinetic_variable_single(double E_MeV_u);

/**
 * Computes the cross section (in 1/m2) the a particle is scattered
 * in the solid angle O = 2 * pi * theta * d_theta given the
 * scatter angle theta
 *
 * @param[in]  E_MeV_u      			energy of incoming particle
 * @param[in]  particle_no  			particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no  			material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]  n						number of scattering angles given
 * @param[in]  scattering_angle			scattering angles theta (array of size n)
 * @param[out] scatter_cross_section	scatter cross section (array of size n)
 * @return     status code
 */
long AT_Rutherford_scatter_cross_section( const double E_Mev_u,
		const long particle_no,
		const long material_no,
		const long n,
		const double scattering_angle[],
		double scatter_cross_section[]);

/**
 * Computes the gyroradius for a particle in a magnetic field. For this,
 * the effective charge of the particle (as a function of the kinetic energy)
 * is used.
 *
 * @param[in]  E_MeV_u      			energy of incoming particle
 * @param[in]  particle_no  			particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]  B_T						magnetic B-field strength in Tesla
 * @return     gyroradius in m
 */
double AT_gyroradius_m( const double E_MeV_u,
		const long particle_no,
		const double B_T);


#endif /* AT_PHYSICSROUTINES_H_ */
