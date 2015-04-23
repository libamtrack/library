#ifndef AT_STOPPINGPOWERDATAPSTAR_H_
#define AT_STOPPINGPOWERDATAPSTAR_H_

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


#include "AT_DataParticle.h"
#include "AT_DataMaterial.h"
#include "AT_PhysicsRoutines.h"

int AT_PSTAR_wrapper( const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const char info[],
		double mass_stopping_power_MeV_cm2_g[]);

/** PSTAR data downloaded from NIST database: http://www.nist.gov/pml/data/star/index.cfm
 * Stopping-Power and Range Tables
 * for Electrons, Protons, and Helium Ions
 *
 * M.J. Berger,1) J.S. Coursey,2) M.A. Zucker2) and J. Chang2)
 *
 * 1) NIST, Physics Laboratory, Ionizing Radiation Division
 * 2) NIST, Physics Laboratory, ECSED
*/

/**
 * Basic structure to hold stopping power data for
 * a specified data source, material, and particle
 * The struct makes no assumption on the nature and the units
 * of energy and stopping power
 *
 * @struct AT_stopping_power_tabulated_source
 */
typedef struct {
	const long   number_of_data_points; /**< number of data points for given material and source */
	const long   material_no;
	const long   particle_no;
	const double energy_and_stopping_power[][2];
} PSTAR_data_for_material_struct;

/**
 * Structure to hold the single stopping power data structs for a given material
 *
 * @struct 	AT_stopping_power_tabulated_source_group_for_all_materials_struct
 * @see 	AT_stopping_power_tabulated_source
 */
typedef struct {
	const long number_of_materials; /**< number of data points for given source */
	const long material_no[MATERIAL_DATA_N];
	const PSTAR_data_for_material_struct * stopping_power_source_data[];
} PSTAR_data_struct;

#endif /* AT_STOPPINGPOWERDATAPSTAR_H_ */
