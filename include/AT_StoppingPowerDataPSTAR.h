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

/**
 * Returns the electronic mass stopping power in MeV*cm2/g
 * as given by the NIST PSTAR tables for a number of materials.
 * The data are scaled by the effective charge if projectiles
 * other than protons are requested.
 *
 * Data were downloaded from the NIST website:
 * http://www.nist.gov/pml/data/star/index.cfm
 *
 * See: Stopping-Power and Range Tables
 * for Electrons, Protons, and Helium Ions
 *
 * M.J. Berger,1) J.S. Coursey,2) M.A. Zucker2) and J. Chang2)
 *
 * 1) NIST, Physics Laboratory, Ionizing Radiation Division
 * 2) NIST, Physics Laboratory, ECSED
 *
 *
 * @param[in]   n							    number of energies / particles
 * @param[in]   E_MeV_u					    	kinetic energies in MeV per amu (array of size n)
 * @param[in]   particle_no                     particle numbers (array of size n)
 * @param[in]   material_no                 	material number
 * @param[in]   info							not used
 * @param[out]  mass_stopping_power_MeV_cm2_g   array to return stopping powers (array of size n)
 * @return
 */
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

#define N_PSTAR_DATAPOINTS   132
#define N_PSTAR_MATERIALS    10

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
	const long   particle_no;
	const double energy_and_stopping_power[N_PSTAR_DATAPOINTS][2];
} PSTAR_data_for_material_struct;

/**
 * Structure to hold the single stopping power data structs for a given material
 *
 * @struct 	AT_stopping_power_tabulated_source_group_for_all_materials_struct
 * @see 	AT_stopping_power_tabulated_source
 */
typedef struct {
	const long material_no[N_PSTAR_MATERIALS];
	const PSTAR_data_for_material_struct * stopping_power_source_data[N_PSTAR_MATERIALS];
} PSTAR_data_struct;

#endif /* AT_STOPPINGPOWERDATAPSTAR_H_ */
