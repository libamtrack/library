#ifndef AT_ENERGYLOSSSTRAGGLING_H_
#define AT_ENERGYLOSSSTRAGGLING_H_

/**
 * @brief Stopping power
 */

/*
 *    AT_EnergyLoss.h
 *    ===============
 *
 *    Copyright 2006, 2014 The libamtrack team
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

#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "AT_Error.h"
#include "AT_DataMaterial.h"
#include "AT_PhysicsRoutines.h"


#define BETHE_LOWER_LIMIT_E_MEV_U 1.0

#ifdef HAVE_CERNLIB

/**
 * Compute Vavilov distribution using CERNLIB (G116)
 *
 * @param[in]  n             array size
 * @param[in]  lambda_V      lambda (array of size n)
 * @param[in]  kappa         straggling parameter
 * @param[in]  beta          relativistic speed, between 0 and 1
 * @param[out] density       resulting density (array of size n)
 */
void AT_Vavilov_PDF( const long n, const double lambda_V[], const double kappa, const double beta,
		double density[]);

/**
 * Compute Landau distribution using CERNLIB (G115)
 *
 * @param[in]  n             array size
 * @param[in]  lambda        lambda (array of size n)
 * @param[out] density       resulting density (array of size n)
 */
void AT_Landau_PDF( const long n, const double lambda[], double density[]);

#endif /* HAVE_CERNLIB */

/**
 * Computes the Rutherford single differential cross section
 * for the energy spectrum of secondary electrons produced by
 * an HCP
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @param[in]      material_no  material index
 * @param[in]  	   n      		number of secondary electron energies
 * @param[in]      T_MeV 	    electron energies (array of size n)
 * @param[out]     dsdT_m2_MeV  Rutherford SDCS for given electron energies (array of size n)
 * @return         status code
 */
int AT_Rutherford_SDCS( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const long n,
		const double T_MeV[],
		double dsdT_m2_MeV[]);

/**
 * Computes leading term of the Bethe formula
 * for many particles according to ICRU49, p.6,
 * after Cohen and Taylor (1986)
 * @param[in]  	    E_MeV_u      energies of particle per nucleon
 * @param[in]  	    particle_no  particle indices
 * @see             AT_DataParticle.h for definition
 * @param[in]       material_no  material index
 * @see             AT_DataMaterial.h for definition
 * @param[in]       use_effective_charge 	if true the effective projectile charge (using the Barkas parametrization) will be used instead of the atomic number
 * @return			result
 */
double AT_el_energy_loss_leading_term_MeV_cm2_g(	const double 	E_MeV_u,
						const long 		particle_no,
						const long 		material_no,
						const bool		use_effective_charge);

/**
 * Computes the electronic energy loss using the Bethe formula
 * according to ICRU49, p.6, Eq. 2.1
 * BUT WITHOUT shell, Bloch or Barkas correction!
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      E_restricted_keV 	if positive and smaller than maximally transferable energy, the restricted stopping power will be computed
 * @return     result
 */
double AT_Bethe_energy_loss_MeV_cm2_g_single(	const double 	E_MeV_u,
												const long 	    particle_no,
												const long 		material_no,
												const double	E_restricted_keV,
												const bool      use_effective_charge);


/**
 * Computes the mass stopping power using the Bethe formula
 * for many particles according to ICRU49, p.6, Eq. 2.1
 * BUT WITHOUT shell, Bloch or Barkas correction!
 * @param[in]  	   n      		number of particles
 * @param[in]  	   E_MeV_u      energies of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  particle indices (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index (single value)
 * @see          AT_DataMaterial.h for definition
 * @param[in]      E_restricted_keV 	if positive and smaller than maximally transferable energy, the restricted stopping power will be computed (single value)
 * @param[out]     Mass_Stopping_Power_MeV_cm2_g (array of size n)
 */
void AT_Bethe_energy_loss_MeV_cm2_g(	const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const double E_restricted_keV,
		const bool  use_effective_charge,
		double Mass_Stopping_Power_MeV_cm2_g[]);

/**
 * Computes the mean energy loss in a slab of
 * material using the Bethe formula
 * for many particles according to ICRU49
 * BUT WITHOUT shell, Bloch or Barkas correction!
 * No effective projectile charge is considered!
 * @param[in]  	   E_MeV_u      energies of particle per nucleon
 * @param[in]  	   particle_no  particle indices
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @return     result
 */
double AT_Bethe_mean_energy_loss_MeV( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);

/**
 * Computes the stopping number to be used with the Bethe formula
 * according to ICRU49, p.6, Eq. 2.3
 * BUT WITHOUT shell correction!
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      E_restricted_keV 	if positive and smaller than maximally transferable energy, the restricted stopping number will be computed
 * @return     result
 */
double AT_Bethe_Stopping_Number(	const double 	E_MeV_u,
									const long      particle_no,
									const long 		material_no,
									const double	E_restricted_keV);

/**
 * Computes the kappa criterium for the
 * energy loss distribution according to
 * Seltzer & Berger
 * No effective projectile charge is considered!
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @return     result
 */
double AT_kappa( const double 	E_MeV_u,
		const long      particle_no,
		const long 		material_no,
		const double    slab_thickness_um);


#ifdef HAVE_CERNLIB
/**
 * Computes the lambda parameter for the
 * Vavilov distribution acc. to Seltzer & Berger
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   energy_loss_keV      energy loss (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon
 * @param[in]  	   particle_no  		particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @param[out]     lambda_V (array of size n)
 */
void AT_lambda_from_energy_loss( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double lambda_V[]);

/**
 * Computes the energy loss from the lambda parameter of the
 * Vavilov distribution acc. to Seltzer & Berger
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   lambda_V      energy loss (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon
 * @param[in]  	   particle_no  		particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @param[out]     energy_loss_keV (array of size n)
 */
void AT_energy_loss_from_lambda( const long n,
		const double lambda_V[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double energy_loss_keV[]);

/**
 * Computes the energy loss
 * Vavilov distribution acc. to Seltzer & Berger
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   energy_loss_keV      energy loss (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon
 * @param[in]  	   particle_no  		particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @param[out]     fDdD (array of size n)
 */
void AT_Vavilov_energy_loss_distribution( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double fDdD[]);

/**
 * Computes the energy loss distribution
 * Uses Landau, Vavilov or Gauss depending on
 * kappa parameter
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   energy_loss_keV      energy loss (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon
 * @param[in]  	   particle_no  		particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @param[out]     fDdD (array of size n)
 */
void AT_energy_loss_distribution( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double fDdD[]);

/**
 * Computes the most probable energy loss
 * Uses Landau, Vavilov or Gauss depending on
 * kappa parameter
 * No effective projectile charge is considered!
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @return     result
 */
double AT_energy_loss_mode( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);

/**
 * Computes the width of the energy loss distribution
 * Uses Landau, Vavilov or Gauss depending on
 * kappa parameter
 * No effective projectile charge is considered!
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @return     result
 */
double AT_energy_loss_FWHM( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);

double AT_lambda_Vavilov_Mode(const double kappa, const double beta);
double AT_lambda_Vavilov_Mean(const double kappa, const double beta);
double AT_lambda_Vavilov_Variance(const double kappa, const double beta);
double AT_lambda_Vavilov_Skewness(const double kappa, const double beta);

double AT_lambda_Vavilov_FWHM_left(const double kappa, const double beta);
double AT_lambda_Vavilov_FWHM_right(const double kappa, const double beta);
double AT_lambda_Vavilov_FWHM(const double kappa, const double beta);
double AT_energy_loss_keV_Vavilov_FWHM(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);

double AT_lambda_Landau_Mode();
double AT_lambda_Landau_Mean(const double kappa, const double beta);
double AT_lambda_Landau_FWHM_left();
double AT_lambda_Landau_FWHM_right();
double AT_lambda_Landau_FWHM();
double AT_energy_loss_keV_Landau_FWHM(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);
double AT_energy_loss_keV_Landau_Mode(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);
#endif /* HAVE_CERNLIB */

double AT_Gauss_Mode();
double AT_Gauss_Mean();
double AT_Gauss_FWHM();

#endif /* AT_ENERGYLOSSSTRAGGLING_H_ */
