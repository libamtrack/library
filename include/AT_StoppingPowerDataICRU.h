#ifndef AT_STOPPINGPOWERDATAICRU_H_
#define AT_STOPPINGPOWERDATAICRU_H_

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

#include "AT_DataMaterial.h"

/**
 * Returns the electronic mass stopping power in MeV*cm2/g
 * as given by ICRU49 (H+He) and ICRU 73 (Li...Ar). The data
 * have be transformed to fit the units (E/u and MeV*cm2/g).
 *
 * @param[in]   n							    number of energies / particles
 * @param[in]   E_MeV_u					    	kinetic energies in MeV per amu (array of size n)
 * @param[in]   particle_no                     particle numbers (array of size n)
 * @param[in]   material_no                 	material number
 * @param[in]   info							not used
 * @param[out]  mass_stopping_power_MeV_cm2_g   array to return stopping powers (array of size n)
 * @return
 */
int AT_ICRU_wrapper( const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const char info[],
		double mass_stopping_power_MeV_cm2_g[]);


/**
 * @struct AT_stopping_power_ICRU_table
 * TODO
 */
typedef struct {
	const long material_no;
	const long number_of_data_points;
	const double energy_and_stopping_power[19][53];
} AT_stopping_power_ICRU_table_struct;

#endif /* AT_STOPPINGPOWERDATAICRU_H_ */
