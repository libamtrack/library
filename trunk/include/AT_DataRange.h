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
#include "AT_DataStoppingPower.h"



#endif /* AT_DATARANGE_H_ */

/**
 * Structure to hold parameters other then parameter to integrate over
 * for GSL integration routine used in AT_CSDA_range_Bethe_cm2_g_single
 */
typedef struct {
  long    particle_no;
  long    material_no;
} AT_CSDA_range_Bethe_parameters;

double AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_int( double  r_m,
    void*   params);

/**
 * Computes the CSDA range using the Bethe formula (AT_Stopping_Power_Bethe_Number)
 * according to ICRU49, p.6, Eq. 2.1
 * BUT WITHOUT shell or density, Bloch or Barkas correction!
 * @param[in]  	   n            number of particles
 * @param[in]  	   E_MeV_u      energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out]    CSDA_range_cm2_g resulting range (array of size n)
 */
void AT_CSDA_range_Bethe_cm2_g_multi(	const long    n,
		                    const double 	E_MeV_u[],
							const long 		particle_no[],
							const long 		material_no,
							double          CSDA_range_cm2_g[]);

/**
 * Computes the CSDA range using the Bethe formula (AT_Stopping_Power_Bethe_Number)
 * according to ICRU49, p.6, Eq. 2.1
 * BUT WITHOUT shell or density, Bloch or Barkas correction!
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @return     result
 */
double AT_CSDA_range_Bethe_cm2_g_single(	const double 	E_MeV_u,
							const long 		particle_no,
							const long 		material_no);

/**
 * Computes the water equivalent path length (WEPL) using the Bethe formula
 * @param[in]  	   n            number of particles
 * @param[in]  	   E_MeV_u      energy of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  particle index (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[out]    WEPL resulting water equivalent path length (array of size n)
 */
void AT_WEPL_Bethe_multi(	const long    n,
		const double 	E_MeV_u[],
		const long 		particle_no[],
		const long 		material_no,
		double          WEPL[]);

/**
 * Computes the water equivalent path length (WEPL) using the Bethe formula
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @return     result
 */
double AT_WEPL_Bethe_single(	const double 	E_MeV_u,
		const long 		particle_no,
		const long 		material_no);

