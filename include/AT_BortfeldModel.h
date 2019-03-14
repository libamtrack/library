#ifndef AT_BortfeldModel_H_
#define AT_BortfeldModel_H_

/**
 * @brief Bortfeld model of Bragg curve
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

#include "AT_BortfeldModel.h"

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
#include <sys/malloc.h>
#else

#include <malloc.h>

#endif


/**
 * Computes dose at given depth for proton beams according to analytical model of T. Bortfeld
 * Bortfeld, 1997, An analytical approximation of the Bragg curve for therapeutic
 * proton beams, Med. Phys. 24(12), 2024ff.
 * @param[in] z_cm            depth in medium [cm]
 * @param[in] fluence_cm2     proton fluence [1/cm2]
 * @param[in] E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[in] eps             fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
 * @return                    dose at given depth [Gy]
 */
double AT_dose_Bortfeld_Gy_single(const double z_cm,
                                  const double fluence_cm2,
                                  const double E_MeV_u,
                                  const double sigma_E_MeV_u,
                                  const long material_no,
                                  const double eps);

/**
 * Computes dose at given depth for proton beams according to analytical model of T. Bortfeld
 * Bortfeld, 1997, An analytical approximation of the Bragg curve for therapeutic
 * proton beams, Med. Phys. 24(12), 2024ff.
 * @param[in]  n               number of depth steps
 * @param[in]  z_cm            depths in medium [cm] (array of size n)
 * @param[in]  fluence_cm2     proton fluence [1/cm2]
 * @param[in]  E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[in] eps             fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
 * @param[out] dose_Gy         doses at given depth [Gy] (array of size n)
 */
void AT_dose_Bortfeld_Gy_multi(const long n,
                               const double z_cm[],
                               const double fluence_cm2,
                               const double E_MeV_u,
                               const double sigma_E_MeV_u,
                               const long material_no,
                               const double eps,
                               double dose_Gy[]);


/**
 * Computes track averaged LET according to Wilkens model
 * @param[in] z_cm            depth in medium [cm]
 * @param[in] E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @return                    track averaged LET at given depth [keV/um]
 */
double AT_LET_t_Wilkens_keV_um_single(const double z_cm,
        const double E_MeV_u,
        const double sigma_E_MeV_u,
        const long material_no);


/**
 * Computes track averaged LET according to Wilkens model
 * @param[in]  n               number of depth steps
 * @param[in]  z_cm            depths in medium [cm] (array of size n)
 * @param[in]  E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in]  sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in]  material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[out] LET_keV_um      track averaged LET at given depth [keV/um] (array of size n)
 */
void AT_LET_t_Wilkens_keV_um_multi(const long n,
                               const double z_cm[],
                               const double E_MeV_u,
                               const double sigma_E_MeV_u,
                               const long material_no,
                               double LET_keV_um[]);

#endif /* AT_BortfeldModel_H_ */