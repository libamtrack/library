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


int AT_Bethe_wrapper( const int n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double mass_stopping_power_MeV_cm2_g[]);



#endif /* AT_STOPPINGPOWERDATABETHE_H_ */
