#ifndef AT_STOPPINGPOWER_H_
#define AT_STOPPINGPOWER_H_

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

#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "AT_Error.h"
#include "AT_NumericalRoutines.h"
#include "AT_PhysicsRoutines.h"
#include "AT_EnergyLoss.h"

#include "AT_StoppingPowerData.h"

/**
 * Retrieves the mass stopping power in MeV*cm2/g for
 * the requested energies and particles for a specified
 * material and data source
 * @param[in] stopping_power_source_no
 * @param[in] n
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @param[out] stopping_power_MeV_cm2_g
 * @return
 */
int AT_Mass_Stopping_Power( const long stopping_power_source_no,
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_MeV_cm2_g[]);

/**
 * Retrieves the stopping power in keV/um for
 * the requested energies and particles for a specified
 * material and data source
 * @param[in] stopping_power_source_no
 * @param[in] n
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @param[out] stopping_power_keV_um
 * @return
 */
int AT_Stopping_Power( const long stopping_power_source_no,
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_keV_um[]);
/**
 * TODO
 * @param[in] stopping_power_source_no
 * @param[in] Stopping_Power_MeV_cm2_g
 * @param[in] particle_no
 * @param[in] material_no
 * @return range [m]
 */
double AT_Energy_MeV_u_from_Stopping_Power_single( const long stopping_power_source_no,
		const double Stopping_Power_MeV_cm2_g, const long particle_no,
		const long material_no);



#endif /* AT_STOPPINGPOWER_H_ */
