#ifndef AT_BortfeldModel_H_
#define AT_BortfeldModel_H_

/**
 * @brief Numerical Routines
 */


/*
 *    AT_BortfeldModel.h
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

#include "AT_NumericalRoutines.h"

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
#include <sys/malloc.h>
#else

#include <malloc.h>

#endif


/**
 * TODO
 * @param[in] z_cm
 * @param[in] E_MeV_u
 * @param[in] fluence_cm2
 * @param[in] sigma_E_MeV_u
 * @param[in] material_no
 * @param[in] eps
 * @return
 */
double AT_dose_Bortfeld_Gy_single(const double z_cm,
                                  const double E_MeV_u,
                                  const double fluence_cm2,
                                  const double sigma_E_MeV_u,
                                  const long material_no,
                                  const double eps);

/**
 * TODO
 * @param[in] n
 * @param[in] z_cm TODO (array of size n)
 * @param[in] E_MeV_u
 * @param[in] fluence_cm2
 * @param[in] sigma_E_MeV_u
 * @param[in] material_no
 * @param[in] eps
 * @param[out] dose_Gy TODO (array of size n)
 */
void AT_dose_Bortfeld_Gy_multi(const long n,
                               const double z_cm[],
                               const double E_MeV_u,
                               const double fluence_cm2,
                               const double sigma_E_MeV_u,
                               const long material_no,
                               const double eps,
                               double dose_Gy[]);

#endif /* AT_BortfeldModel_H_ */