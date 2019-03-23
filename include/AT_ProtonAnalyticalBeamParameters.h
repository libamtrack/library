#ifndef AT_ProtonAnalyticalBeamParameters_H_
#define AT_ProtonAnalyticalBeamParameters_H_

/**
 * @brief Proton analytical models of dose, LET and RBE
 */


/*
 *    AT_ProtonAnalyticalBeamParameters.h
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
#include <sys/malloc.h>
#else

#include <malloc.h>

#endif


#include "AT_ProtonAnalyticalModels.h"


/**
 * Location of maximum dose according to Bortfeld dose model
 * @see AT_dose_Bortfeld_Gy_single
 * @param[in] E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[in] eps             fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
 * @return                    depth at which dose reaches maximum value [cm]
 */
double AT_max_location_Bortfeld_cm(const double E_MeV_u,
                                   const double sigma_E_MeV_u,
                                   const long material_no,
                                   const double eps);

/**
 * Range defined as the location on distal part where dose drops to certain fraction of maximum value.
 * @see AT_dose_Bortfeld_Gy_single for details of Bortfeld dose models
 * @param[in] E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[in] eps             fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
 * @param[in] dose_drop       fraction of max dose at which range is calculated
 * if negative a default value of 0.8 is assumed
 * @param[in] search_direction   is positive a search is done in distal part (behind Bragg peak maximum), otherwise
 * search is done in proximal part
 * if negative a default value of 0.8 is assumed
 * @return                    range [cm]
 */
double AT_range_Bortfeld_cm(const double E_MeV_u,
                            const double sigma_E_MeV_u,
                            const long material_no,
                            const double eps,
                            const double dose_drop,
                            const short search_direction);


/**
 * Full width at half-maximum (FWHM) defined as the width of the dose profile at 50% of maximum dose value.
 * @see AT_dose_Bortfeld_Gy_single for details of Bortfeld dose models
 * @param[in] E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[in] eps             fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
 * @return                    FWHM [cm]
 */
double AT_fwhm_Bortfeld_cm(const double E_MeV_u,
                           const double sigma_E_MeV_u,
                           const long material_no,
                           const double eps);

/**
 * Ratio between maximum dose value and the entrance dose (so called "max_to_plateau" ratio).
 * @see AT_dose_Bortfeld_Gy_single for details of Bortfeld dose models
 * @param[in] E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[in] eps             fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
 * @return                    "max_to_plateau" ratio
 */
double AT_max_plateau_Bortfeld(const double E_MeV_u,
                               const double sigma_E_MeV_u,
                               const long material_no,
                               const double eps);

/**
 * Energy of the beam which has given range.
 * @see AT_dose_Bortfeld_Gy_single for details of Bortfeld dose models
 * @param[in] range_cm        range [cm]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[in] eps             fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
 * @param[in] dose_drop       fraction of max dose at which range is calculated
 * if negative a default value of 0.8 is assumed
 * @return                    energy [MeV/u]
 */
double AT_energy_Bortfeld_MeV_u(const double range_cm,
                                const double sigma_E_MeV_u,
                                const long material_no,
                                const double eps,
                                const double dose_drop);

/**
 * @see AT_dose_Bortfeld_Gy_single for details of Bortfeld dose models
 * @param[in] range_cm        range [cm]
 * @param[in] fwhm_cm         FWHM [cm]
 * @param[in] max_to_plateau  "max_to_plateau" ratio
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[in] dose_drop       fraction of max dose at which range is calculated
 * if negative a default value of 0.8 is assumed
 * @param[out] E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[out] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * @param[out] eps             fraction of primary fluence contributing to the tail of energy spectrum
 */
void AT_fit_Bortfeld(const double range_cm,
                     const double fwhm_cm,
                     const double max_to_plateau,
                     const long material_no,
                     const double dose_drop,
                     double * E_MeV_u,
                     double * sigma_E_MeV_u,
                     double * eps);


#endif /* AT_ProtonAnalyticalBeamParameters_H_ */
