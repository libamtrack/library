#ifndef AT_DATARANGE_H_
#define AT_DATARANGE_H_

/**
 * @brief Range
 */

/*
 *    AT_DataRange.h
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

#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "AT_Error.h"
#include "AT_NumericalRoutines.h"
#include "AT_PhysicsRoutines.h"
#include "AT_StoppingPower.h"



#endif /* AT_DATARANGE_H_ */

/**
 * Structure to hold parameters other then parameter to integrate over
 * for GSL integration routine used in AT_CSDA_range_Bethe_g_cm2_single
 */
typedef struct {
  long    particle_no;
  long    material_no;
} AT_CSDA_range_Bethe_parameters;

/**
 * Integrand function for CSDA computation:
 * inverse stopping power as function of energy (all other dependencies [particle, material] in params structre)
 * @param[in] E_MeV_u
 * @param[in] params
 * @return TODO
 */
double AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_int( double  E_MeV_u,
    void*   params);

/**
 * Structure to hold parameters
 */
typedef struct {
  double  E_initial_MeV_u;
  long    particle_no;
  long    material_no;
  double  range_g_cm2;
} AT_CSDA_range_difference_parameters;


/**
 * Integrand function for CSDA computation:
 * inverse stopping power as function of energy (all other dependencies [particle, material] in params structre)
 * @param[in] E_MeV_u
 * @param[in] params
 * @return TODO
 */
double AT_CSDA_range_g_cm_int( double  E_MeV_u,
    void*   params);

/**
 * Computes the CSDA range using the Bethe formula (AT_Stopping_Power_Bethe_Number)
 * according to ICRU49, p.6, Eq. 2.1
 * BUT WITHOUT shell or density, Bloch or Barkas correction!
 * @param[in]  	   n                    number of particles
 * @param[in]  	   E_initial_MeV_u      initial energy of particle per nucleon (array of size n)
 * @param[in]  	   E_final_MeV_u        final energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out]    CSDA_range_cm2_g resulting range (array of size n)
 */
void AT_CSDA_range_Bethe_g_cm2_multi(	const long    n,
		const double 	E_initial_MeV_u[],
		const double 	E_final_MeV_u[],
		const long 		particle_no[],
		const long 		material_no,
		double          CSDA_range_cm2_g[]);

/**
 * Computes the CSDA range using the Bethe formula (AT_Stopping_Power_Bethe_Number)
 * according to ICRU49, p.6, Eq. 2.1
 * BUT WITHOUT shell or density, Bloch or Barkas correction!
 * @param[in]  	   E_initial_MeV_u      initial energy of particle per nucleon
 * @param[in]  	   E_final_MeV_u       final energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @return     result
 */
double AT_CSDA_range_Bethe_g_cm2_single(	const double 	E_initial_MeV_u,
		const double 	E_final_MeV_u,
		const long 		particle_no,
		const long 		material_no);

/**
 * Solver function for CSDA energy after slab
 * @param[in] E_final_MeV_u
 * @param[in] params
 * @return TODO
 */
double AT_CSDA_range_difference_solver( double  E_final_MeV_u,
	    void*   params);

/**
 * Computes the ion energy after transversing a slab of material using Bethe stopping power
 * and CSDA approach
 * @param[in]  	   E_initial_MeV_u      initial energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no          particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no          material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_m     thickness of slab to transversed
 */
double AT_CSDA_energy_after_slab_E_MeV_u_single( const double E_initial_MeV_u,
		const long   particle_no,
		const long   material_no,
		const double slab_thickness_m);

/**
 * Computes the ion energy after transversing a slab of material using Bethe stopping power
 * and CSDA approach for many energies / particles
 *
 * @param[in]  	   n                    number of particles
 * @param[in]  	   E_initial_MeV_u      initial energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no          particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no          material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_m     thickness of slab to transversed
 * @param[out]     E_final_MeV_u        final energy after slab (array of size n)
 */
void AT_CSDA_energy_after_slab_E_MeV_u_multi( const long n,
		const double E_initial_MeV_u[],
		const long   particle_no[],
		const long   material_no,
		const double slab_thickness_m,
		double E_final_MeV_u[]);

/**
 * Computes the water equivalent path length (WEPL) using the Bethe formula
 * @param[in]  	   n            number of particles
 * @param[in]  	   E_MeV_u      energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_m  thickness of slab of material different than water, in meter
 * @param[out]    WEPL resulting water equivalent path length (array of size n)
 */
void AT_WEPL_Bethe_multi(	const long    n,
		const double    E_MeV_u[],
		const long 		particle_no[],
		const long 		material_no,
		const double    slab_thickness_m,
		double          WEPL[]);

/**
 * Computes the water equivalent path length (WEPL) using the Bethe formula
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      slab_thickness_m  thickness of slab of material different than water, in meter
 * @return     result
 */
double AT_WEPL_Bethe_single(	const double 	E_MeV_u,
		const long 		particle_no,
		const long 		material_no,
		const double    slab_thickness_m);

