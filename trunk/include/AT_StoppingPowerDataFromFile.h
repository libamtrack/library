#ifndef AT_STOPPINGPOWERDATAFROMFILE_H_
#define AT_STOPPINGPOWERDATAFROMFILE_H_

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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/**
 * Wrapper for stopping power data read from file
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return stopping power (MeV cm2 per g)
 */
double AT_FromFile_wrapper( const double E_MeV_u, const long particle_no,
		const long material_no);



#endif /* AT_STOPPINGPOWERDATAFROMFILE_H_ */
