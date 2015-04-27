#ifndef AT_STOPPINGPOWER_H_
#define AT_STOPPINGPOWER_H_

/**
 * @brief Main routines for stopping power data
 */

/*
 *    AT_StoppingPower.h
 *    ==================
 *
 *    Copyright 2006, 2015 The libamtrack team
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
 * Retrieves the electronic mass stopping power in MeV*cm2/g
 * for the requested energies and particles for a specified
 * material and data source. The data source is thereby
 * given via its name (s. AT_StoppingPowerData.h from
 * details), except for data that should be read for
 * a file, in this case the (path and) filename has to be
 * provided. In this case, the user has to make sure that
 * energy and stopping power units are correct and that
 * the data match the given material (use material.no = 0
 * for custom-defined material).
 *
 * The file has to be plain
 * ASCII with three columns (separated by space)
 *    charge, energy, and stopping power
 * and sorted in ascending order by first charge than energy
 * any alphanumeric comment can be inserted (in separate
 * lines)
 *
 * @param[in]      stopping_power_source       name of the data source
 * @param[in]      n		               number of energies / particles
 * @param[in]      E_MeV_u                     kinetic energies in MeV per amu (array of size n)
 * @param[in]      particle_no                 particle numbers (array of size n)
 * @param[in]      material_no                 material number
 * @param[out]     stopping_power_MeV_cm2_g    array to return stopping powers (array of size n)
 * @return         status
 */
int AT_Mass_Stopping_Power( const char stopping_power_source[],
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_MeV_cm2_g[]);


/**
 * Retrieves the electronic stopping power in keV/um for
 * the requested energies and particles for a specified
 * material and data source. The data source is thereby
 * given via its name (s. AT_StoppingPowerData.h for
 * details), except for data that should be read from
 * a file, in this case the (path and) filename has to be
 * provided. In this case, the user has to make sure that
 * energy and stopping power units are correct and that
 * the data match the given material (use material.no = 0
 * for custom-defined material) for density scaling.
 * 
 * The file has to be plain
 * ASCII with three columns (separated by space)
 *    charge, energy, and stopping power
 * and sorted in ascending order by first charge than energy
 * any alphanumeric comment can be inserted (in separate
 * lines)
 *
 * @param[in]   stopping_power_source		name of the data source
 * @param[in]   n							number of energies / particles
 * @param[in]   E_MeV_u						kinetic energies in MeV per amu (array of size n)
 * @param[in]   particle_no                 particle numbers (array of size n)
 * @param[in]   material_no                 material number
 * @param[out]  stopping_power_keV_um       array to return stopping powers (array of size n)
 * @return		status
 */
int AT_Stopping_Power( const char stopping_power_source[],
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_keV_um[]);


/**
 * Retrieves the electronic mass stopping power in MeV*cm2/g
 * for the requested energies and particles for a specified
 * material and data source. The data source is thereby
 * given via its integer id (s. AT_StoppingPowerData.h for
 * details). Data that should be read from a file
 * cannot be used with this method.
 *
 * @param[in]   stopping_power_source_no	id of the data source
 * @param[in]   n							number of energies / particles
 * @param[in]   E_MeV_u						kinetic energies in MeV per amu (array of size n)
 * @param[in]   particle_no                 particle numbers (array of size n)
 * @param[in]   material_no                 material number
 * @param[out]  stopping_power_MeV_cm2_g    array to return stopping powers (array of size n)
 * @return		status
 */
int AT_Mass_Stopping_Power_with_no( const long stopping_power_source_no,
		const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double stopping_power_MeV_cm2_g[]);

/**
 * Retrieves the electronic stopping power in keV/um for
 * the requested energies and particles for a specified
 * material and data source. The data source is thereby
 * given via its integer id (s. AT_StoppingPowerData.h for
 * details). Data that should be read from a file
 * cannot be used with this method.
 *
 * @param[in]   stopping_power_source_no	id of the data source
 * @param[in]   n							number of energies / particles
 * @param[in]   E_MeV_u						kinetic energies in MeV per amu (array of size n)
 * @param[in]   particle_no                 particle numbers (array of size n)
 * @param[in]   material_no                 material number
 * @param[out]  stopping_power_keV_um       array to return stopping powers (array of size n)
 * @return		status
 */
int AT_Stopping_Power_with_no( const long stopping_power_source_no,
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
