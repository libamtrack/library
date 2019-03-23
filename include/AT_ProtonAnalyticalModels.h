#ifndef AT_ProtonAnalyticalModels_H_
#define AT_ProtonAnalyticalModels_H_

/**
 * @brief Proton analytical models of dose, LET and RBE
 */


/*
 *    AT_ProtonAnalyticalModels.h
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
#include <sys/malloc.h>
#else

#include <malloc.h>

#endif

#include "AT_DataMaterial.h"



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

/**
 * Computes dose averaged LET according to Wilkens model
 * @param[in] z_cm            depth in medium [cm]
 * @param[in] E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @return                    track averaged LET at given depth [keV/um]
 */
double AT_LET_d_Wilkens_keV_um_single(const double z_cm,
                                      const double E_MeV_u,
                                      const double sigma_E_MeV_u,
                                      const long material_no);


/**
 * Computes dose averaged LET according to Wilkens model
 * @param[in]  n               number of depth steps
 * @param[in]  z_cm            depths in medium [cm] (array of size n)
 * @param[in]  E_MeV_u         initial kinetic energy of proton beam [MeV/u]
 * @param[in]  sigma_E_MeV_u   kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in]  material_no     material code number
 * @see          AT_DataMaterial.h for definition
 * @param[out] LET_keV_um      track averaged LET at given depth [keV/um] (array of size n)
 */
void AT_LET_d_Wilkens_keV_um_multi(const long n,
                                   const double z_cm[],
                                   const double E_MeV_u,
                                   const double sigma_E_MeV_u,
                                   const long material_no,
                                   double LET_keV_um[]);

/**
 *Proton RBE models code numbers
 */
enum AT_RBEModels {
    RBE_One = 1,     /**< RBE equals 1.0 */
    RBE_OnePointOne = 2,     /**< RBE equals 1.1 */
    RBE_Carabe = 3,     /**< Carabe model (Phys. Med. Biol. 57(5), 1159-1172 (2012)) */
    RBE_Wedenberg = 4,     /**< Wedenberg model (Acta Oncol. 52(3), 580-588 (2013)) */
    RBE_McNamara = 5,     /**< McNamara model (Phys. Med. Biol. 60(21), 8399-8416)*/
};


/**
 * Computes proton RBE according to one of several analytical models
 * @param[in] z_cm              depth in medium [cm]
 * @param[in] entrance_dose_Gy  entrance dose [Gy]
 * @param[in] E_MeV_u           initial kinetic energy of proton beam [MeV/u]
 * @param[in] sigma_E_MeV_u     kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
  * @param[in] eps              fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
* @param[in] ref_alpha_beta_ratio   ratio of alpha to beta (parameters in linear-quadratic model) for reference radiation
 * @return                    proton RBE
 */
double AT_proton_RBE_single(const double z_cm,
                            const double entrance_dose_Gy,
                            const double E_MeV_u,
                            const double sigma_E_MeV_u,
                            const double eps,
                            const double ref_alpha_beta_ratio,
                            const int rbe_model_no);


/**
 * Computes proton RBE according to one of several analytical models
 * @param[in]  n                 number of depth steps
 * @param[in]  z_cm              depth in medium [cm] (array of size n)
 * @param[in]  entrance_dose_Gy  entrance dose [Gy]
 * @param[in]  E_MeV_u           initial kinetic energy of proton beam [MeV/u]
 * @param[in]  sigma_E_MeV_u     kinetic energy spread (standard deviation) [MeV/u]
 * if negative a default value of 0.01 * E_MeV_u is assumed
 * @param[in] eps               fraction of primary fluence contributing to the tail of energy spectrum
 * if negative a default value of 0.03 is assumed
 * @param[in]  ref_alpha_beta_ratio   ratio of alpha to beta (parameters in linear-quadratic model) for reference radiation
 * @param[out] rbe             proton RBE (array of size n)
 */
void AT_proton_RBE_multi(const long n,
                         const double z_cm[],
                         const double entrance_dose_Gy,
                         const double E_MeV_u,
                         const double sigma_E_MeV_u,
                         const double eps,
                         const double ref_alpha_beta_ratio,
                         const int rbe_model_no,
                         double rbe[]);



#endif /* AT_ProtonAnalyticalModels_H_ */