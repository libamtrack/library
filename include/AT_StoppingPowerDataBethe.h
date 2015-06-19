#ifndef AT_STOPPINGPOWERDATABETHE_H_
#define AT_STOPPINGPOWERDATABETHE_H_

/**
 * @brief Stopping power
 */

/*
 *    AT_DataStoppingPower.h
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

#include <stdbool.h>

#include "AT_EnergyLoss.h"


#define BETHE_LOWER_LIMIT_E_MEV_U 1.0

/**
 * Returns the electronic mass stopping power in MeV*cm2/g
 * calculated by the Bethe formula, including density correction.
 *
 * @param[in]   n							    number of energies / particles
 * @param[in]   E_MeV_u					    	kinetic energies in MeV per amu (array of size n)
 * @param[in]   particle_no                     particle numbers (array of size n)
 * @param[in]   material_no                 	material number
 * @param[in]   info							not used
 * @param[out]  mass_stopping_power_MeV_cm2_g   array to return stopping powers (array of size n)
 * @return
 */
int AT_Bethe_wrapper( const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const char info[],
		double mass_stopping_power_MeV_cm2_g[]);

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


#endif /* AT_STOPPINGPOWERDATABETHE_H_ */
