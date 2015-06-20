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

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "AT_CernlibFuns.h"

#include "AT_Error.h"
#include "AT_DataMaterial.h"
#include "AT_PhysicsRoutines.h"


/////////////////////////////////////////////////////////////////////////////
// ALL DEFINITIONS AND ROUTINES CONCERNING THE ENERGY LOSS STRAGGLING FOLLOW:
// CERN Program Library Long Writeup W5013, GEANT, 1994.
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
// MEAN ENERGY LOSS
/////////////////////////////////////////////////////////////////////////////


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
double AT_mean_energy_loss_keV( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);


/**
 * Parameter xi - reduced mean energy loss
 * @param[in]  	   E_MeV_u      energies of particle per nucleon
 * @param[in]  	   particle_no  particle indices
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @return			xi
 */
double AT_xi_keV(	const double 	E_MeV_u,
				    const long 		particle_no,
				    const long 		material_no,
				    const double    slab_thickness_um);

/////////////////////////////////////////////////////////////////////////////
// KAPPA PARAMETER
/////////////////////////////////////////////////////////////////////////////

/**
 * Computes the kappa criterium for the
 * energy loss distribution according to
 * Seltzer and Berger, and CERN W5013
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of particles
 * @param[in]  	   E_MeV_u      		energy of particle per amu (array of size n)
 * @param[in]  	   particle_no  		particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um (array of size n)
 * @param[out]	   kappa				kappa parameter (array of size n)
 */
void AT_kappa_multi( const long n,
		const double 	E_MeV_u[],
		const long      particle_no[],
		const long 		material_no,
		const double    slab_thickness_um[],
		double          kappa[]);

/**
 * TODO
 *
 * @param E_MeV_u
 * @param particle_no
 * @param material_no
 * @param slab_thickness_um
 * @return
 */
double AT_kappa_single( const double 	E_MeV_u,
		const long      particle_no,
		const long 		material_no,
		const double    slab_thickness_um);


/////////////////////////////////////////////////////////////////////////////
// ENERGY LOSS STRAGGLING: LANDAU THEORY
/////////////////////////////////////////////////////////////////////////////

/**
 * Computes the Landau probability density function using CERNLIB (G115)
 *
 * @param[in]  n                    array size
 * @param[in]  lambda_landau        Landau lambda (array of size n)
 * @param[out] density              resulting density (array of size n)
 */
void AT_Landau_PDF( const long n, 
        const double lambda_landau[], 
        double density[]);

/**
 * Computes the Landau inverse distribution function using CERNLIB (G115)
 *
 * @param[in]  n                    array size
 * @param[in]  rnd                  random number from uniform distribution between 0 and 1 (array of size n)
 * @param[out] lambda_landau        resulting Landau lambda (array of size n)
 */
void AT_Landau_IDF( const long n, 
        const double rnd[], 
        double lambda_landau[]);

/**
 * Computes the lambda parameter for the
 * Landau distribution acc. to CERN W5013
 *
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   energy_loss_keV      energy loss (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon
 * @param[in]  	   particle_no  		particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @param[out]     lambda_landau (array of size n)
 */
void AT_lambda_landau_from_energy_loss_multi( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double lambda_landau[]);

/**
 *
 * @return
 */
double AT_lambda_landau_from_energy_loss_single( const double energy_loss_keV,
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);


/**
 * Computes the mean lambda, introduced to enable
 * average value for Landau distribution. See Geant3 W5013, p.254
 *
 * @param[in]  	   n      				number of particles
 * @param[in]  	   E_MeV_u      		energy of particle per amu (array of size n)
 * @param[in]  	   particle_no  		particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um (array of size n)
 * @param[out]	   lambda_mean			mean lambda for given particle (array of size n)
 */
void AT_lambda_mean_multi( const long n,
		const double	E_MeV_u[],
		const long      particle_no[],
		const long 		material_no,
		const double    slab_thickness_um[],
		double 			lambda_mean[]);

double AT_lambda_mean_single( const double	E_MeV_u,
		const long      particle_no,
		const long 		material_no,
		const double    slab_thickness_um);


/**
 * Computes the mean lambda, introduced to enable
 * average value for Landau distribution. See Geant3 W5013, p.254
 *
 * @param[in]  	   n      				number of particles
 * @param[in]  	   E_MeV_u      		energy of particle per amu (array of size n)
 * @param[in]  	   particle_no  		particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um (array of size n)
 * @param[out]	   lambda_max			maximum lambda for given particle (array of size n)
 */
void AT_lambda_max_multi( const long n,
		const double	E_MeV_u[],
		const long      particle_no[],
		const long 		material_no,
		const double    slab_thickness_um[],
		double 			lambda_max[]);

double AT_lambda_max_single( double lambda_mean );

double AT_lambda_Landau_Mode();
double AT_lambda_Landau_Mean(const double kappa, const double beta );
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


/**
 * Computes the energy loss from the lambda parameter
 * of the Landau distribution acc. to CERN W5013
 *
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   lambda_landau      Landau lambda (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  		particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um (array of size n)
 * @param[out]     energy_loss_keV (array of size n)
 */
void AT_energy_loss_from_lambda_landau_multi( const long n,
		const double lambda_landau[],
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const double slab_thickness_um[],
		double energy_loss_keV[]);

double AT_energy_loss_from_lambda_landau_single( const double lambda_landau,
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);

/**
 * Computes the energy loss
 * Landau distribution acc. to CERN W5013
 *
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
void AT_Landau_energy_loss_distribution( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double fDdD[]);


/////////////////////////////////////////////////////////////////////////////
// ENERGY LOSS STRAGGLING: VAVILOV THEORY
/////////////////////////////////////////////////////////////////////////////

/**
 * Computes the Vavilov probability density function using CERNLIB (G116)
 *
 * @param[in]  n                   array size
 * @param[in]  lambda_vavilov      Vavilov lambda (array of size n)
 * @param[in]  kappa               straggling parameter
 * @param[in]  beta                relativistic speed, between 0 and 1
 * @param[out] density             resulting density (array of size n)
 */
void AT_Vavilov_PDF( const long n, const double lambda_vavilov[], const double kappa, const double beta,
		double density[]);

/**
 * Computes the Vavilov probability density function using CERNLIB (G116)
 *
 * @param[in]  n                   array size
 * @param[in]  rnd                 random number from uniform distribution between 0 and 1 (array of size n)
 * @param[in]  kappa               straggling parameter  (array of size n)
 * @param[in]  beta                relativistic speed, between 0 and 1 (array of size n)
 * @param[out] lambda_vavilov      resulting Vavilov lambda (array of size n)
 */
void AT_Vavilov_IDF( const long n, const double rnd[], const double kappa[], const double beta[],
		double lambda_vavilov[]);

/**
 * Computes the lambda parameter for the
 * Vavilov distribution acc. to CERN W5013
 *
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   energy_loss_keV      energy loss (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon
 * @param[in]  	   particle_no  		particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um
 * @param[out]     lambda_vavilov (array of size n)
 */
void AT_lambda_vavilov_from_energy_loss_multi( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double lambda_vavilov[]);

/**
 *
 * @return
 */
double AT_lambda_vavilov_from_energy_loss_single( const double energy_loss_keV,
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);


double AT_lambda_Vavilov_Mode(const double kappa, const double beta );
double AT_lambda_Vavilov_Mean(const double kappa, const double beta );
double AT_lambda_Vavilov_Variance(const double kappa, const double beta );
double AT_lambda_Vavilov_Skewness(const double kappa, const double beta );
double AT_lambda_Vavilov_FWHM_left(const double kappa, const double beta );
double AT_lambda_Vavilov_FWHM_right(const double kappa, const double beta );
double AT_lambda_Vavilov_FWHM(const double kappa, const double beta );
double AT_energy_loss_keV_Vavilov_FWHM(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um);

/**
 * Computes the energy loss from the lambda parameter of the
 * Vavilov distribution acc. to CERN W5013
 *
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   lambda_vavilov      Vavilov lambda (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  		particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um (array of size n)
 * @param[out]     energy_loss_keV (array of size n)
 */
void AT_energy_loss_from_lambda_vavilov_multi( const long n,
		const double lambda_vavilov[],
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const double slab_thickness_um[],
		double energy_loss_keV[]);

/**
 * Computes the energy loss
 * Vavilov distribution acc. to CERN W5013
 *
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




/////////////////////////////////////////////////////////////////////////////
// ENERGY LOSS STRAGGLING: GAUSS THEORY
/////////////////////////////////////////////////////////////////////////////

/**
 * Computes Gauss probability density function (for compatibility)
 *
 * @param[in]  n             array size
 * @param[in]  lambda_gauss  Gauss lambda (array of size n)
 * @param[out] density       resulting density (array of size n)
 */
void AT_Gauss_PDF( const long n,
		const double lambda_gauss[],
		double density[]);

/**
 * Compute Gauss inverse distribution function (for compatibility)
 *
 * @param[in]  n             array size
 * @param[in]  rnd           random number from uniform distribution between 0 and 1 (array of size n)
 * @param[out] lambda_gauss  resulting Gauss lambda (array of size n)
 */
void AT_Gauss_IDF( const long n,
		const double rnd[],
		double lambda_gauss[]);


/**
 * Computes the energy loss from the lambda parameter of the
 * Gauss distribution for compatibility with CERN W5013
 *
 * No effective projectile charge is considered!
 * @param[in]  	   n      				number of energy loss data
 * @param[in]  	   lambda_gauss      Gauss lambda (array of size n)
 * @param[in]  	   E_MeV_u      		energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  		particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  		material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_um	slab thickness in um (array of size n)
 * @param[out]     energy_loss_keV (array of size n)
 */
void AT_energy_loss_from_lambda_gauss_multi( const long n,
		const double lambda_gauss[],
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const double slab_thickness_um[],
		double energy_loss_keV[]);


double AT_Gauss_Mode();
double AT_Gauss_Mean();
double AT_Gauss_FWHM();


/////////////////////////////////////////////////////////////////////////////
// ENERGY LOSS STRAGGLING: GENERAL THEORY
/////////////////////////////////////////////////////////////////////////////

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



#endif /* AT_ENERGYLOSSSTRAGGLING_H_ */
