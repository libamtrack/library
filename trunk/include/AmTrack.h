#ifndef AmTrack_H_
#define AmTrack_H_

/**
 * @file
 * @brief libamtrack main file holding the amorphous track routines for RE/RBE calculation
 */

/*
 *    AmTrack.h
 *    =========
 *
 *    Created on: 28.07.2009
 *    Author: greilich
 *
 *
 *    Copyright 2006, 2009 Steffen Greilich / the libamtrack team
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
#include <time.h>

#include "AT_Constants.h"
#include "AT_RDD.h"
#include "AT_SuccessiveConvolutions.h"
#include "AT_GammaResponse.h"
#include "AT_PhysicsRoutines.h"
#include "AT_NumericalRoutines.h"

#include "AT_Wrapper_Java.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>


/**
* Prints name for debugging
*
* @return  none
*/
void AT_PrintName(  void);


/**
* Computes HCP response and RE/RBE using compound Poison process and
* successive convolutions (CPP_SC, the 'SPIFF' algorithm)
*
* @param[in]  n      number of particle types in the mixed particle field (pointer to single variable)
* @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
* @param[in]  particle_no    type of the particles in the mixed particle field (pointer to array of size n)
* @see          AT_DataParticle.h for definition
* @param[in]  fluence_cm2    fluences for the given particles, doses in Gy if negative (pointer to array of size n)
* @param[in]  material_no    index number for detector material (pointer to single variable)
* @see          AT_DataMaterial.h for definition
* @param[in]  RDD_model    index number for chosen radial dose distribution (pointer to single variable)
* @param[in]  RDD_parameters    parameters for chosen radial dose distribution (pointer to array of size depending on chosen model)
* @see          AT_RDD.h for definition
* @param[in]  ER_model    index number for chosen electron-range model (pointer to single variable)
* @param[in]  ER_parameters    parameters for chosen electron-range model (pointer to array of size depending on chosen model)
* @see          AT_ElectronRange.h for definition
* @param[in]  gamma_model    index number for chosen gamma response (pointer to single variable)
* @param[in]  gamma_parameters  parameters for chosen gamma response (pointer to array of size depending on chosen model)
* @see          AT_GammaResponse.h for definition
* @param[in,out]  N2      (algorithm specific) number of bins per factor of two in local dose array (pointer to single variable)
* @param[in]  fluence_factor    factor to scale the fluences given as "fluence_cm2" with (pointer to single variable)
* @param[in]  write_output    if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory (pointer to single variable)
* @param[in]  shrink_tails    (algorithm specific) if true, tails of the local dose distribution, contributing less than "shrink_tails_under" are cut (pointer to single variable)
* @param[in]  shrink_tails_under  (algorithm specific) limit for tail cutting in local dose distribution (pointer to single variable)
* @param[in]  adjust_N2    (algorithm specific) if true, "N2" will be increase if necessary at high fluence to ensure sufficient binning resolution
* @param[in]  lethal_events_mode (algorithm specific) if true, allows to do calculations for cell survival
* @param[out]  results      pointer to array of size 10 to be allocated by the user which will be used to return the results
*    results[0]    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D\n
*    results[1]    d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm\n
*    results[2]    S_HCP           (algorithm independent)  absolute particle response\n
*    results[3]    S_gamma         (algorithm independent)  absolute gamma response\n
*    results[4]    not used        (algorithm independent)\n
*    results[5]    u               (algorithm specific)     mean number of tracks contributing to representative point\n
*    results[6]    u_start         (algorithm specific)     low starting value for mean number of tracks, where linearisation is applied\n
*    results[7]    n_convolutions  (algorithm specific)     number of convolutions performed\n
*    results[8]    not used        (algorithm specific)\n
*    results[9]    not used        (algorithm specific)
* @return  none
*/
void AT_SPIFF(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const float*  fluence_cm2,
    const long*  material_no,
    const long*  RDD_model,
    const float*  RDD_parameters,
    const long*  ER_model,
    const float*  ER_parameters,
    const long*  gamma_model,
    const float*  gamma_parameters,
    long*  N2, // TODO investigate if this can be changed inside
    const float*  fluence_factor,
    const bool*  write_output,
    const bool*  shrink_tails,
    const float*  shrink_tails_under,
    const bool*  adjust_N2,
    const bool*   lethal_events_mode,
    float*  results);

/**
* Computes HCP response and RE/RBE using summation of tracks
* an a Cartesian grid (the 'GSM' algorithm)
*
* @param[in]  n      number of particle types in the mixed particle field (pointer to single variable)
* @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
* @param[in]  particle_no    type of the particles in the mixed particle field (pointer to array of size n)
* @see          AT_DataParticle.h for definition
* @param[in]  fluence_cm2    fluences for the given particles, doses in Gy if negative (pointer to array of size n)
* @param[in]  material_no    index number for detector material (pointer to single variable)
* @see          AT_DataMaterial.h for definition
* @param[in]  RDD_model    index number for chosen radial dose distribution (pointer to single variable)
* @param[in]  RDD_parameters    parameters for chosen radial dose distribution (pointer to array of size depending on chosen model)
* @see          AT_RDD.h for definition
* @param[in]  ER_model    index number for chosen electron-range model (pointer to single variable)
* @param[in]  ER_parameters    parameters for chosen electron-range model (pointer to array of size depending on chosen model)
* @see          AT_ElectronRange.h for definition
* @param[in]  gamma_model    index number for chosen gamma response (pointer to single variable)
* @param[in]  gamma_parameters  parameters for chosen gamma response (pointer to array of size depending on chosen model)
* @see          AT_GammaResponse.h for definition
* @param[in]  N_runs       (algorithm specific) number of runs within which track positions will be resampled
* @param[in]  N2      (algorithm specific) number of bins per factor of two in local dose array (pointer to single variable)
* @param[in]  fluence_factor    factor to scale the fluences given as "fluence_cm2" with (pointer to single variable)
* @param[in]  write_output    if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory (pointer to single variable)
* @param[in]  shrink_tails    (algorithm specific) if true, tails of the local dose distribution, contributing less than "shrink_tails_under" are cut (pointer to single variable)
* @param[in]  shrink_tails_under  (algorithm specific) limit for tail cutting in local dose distribution (pointer to single variable)
* @param[in]  adjust_N2    (algorithm specific) if true, "N2" will be increase if necessary at high fluence to ensure sufficient binning resolution
* @param[in]  lethal_events_mode (algorithm specific) if true, allows to do calculations for cell survival
* @param[in]  nX    (algorithm specific) number of voxels of the grid in x (and y as the grid is quadratic)
* @param[in]  voxel_size_m     side length of a voxel in m
* @param[out]  results      pointer to array of size 10 to be allocated by the user which will be used to return the results
*    results[0]    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D\n
*    results[1]    d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm\n
*    results[2]    S_HCP           (algorithm independent)  absolute particle response\n
*    results[3]    S_gamma         (algorithm independent)  absolute gamma response\n
*    results[4]    n_particles     (algorithm independent)  average number of particle tracks on the detector grid\n
*    results[5]    sd_efficiency   (algorithm specific)     standard deviation for results[0]\n
*    results[6]    sd_d_check      (algorithm specific)     standard deviation for results[1]\n
*    results[7]    sd_S_HCP        (algorithm specific)     standard deviation for results[2]\n
*    results[8]    sd_S_gamma      (algorithm specific)     standard deviation for results[3]\n
*    results[9]    sd_n_particles  (algorithm specific)     standard deviation for results[4]\n
* @return  none
*/

void AT_GSM(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const float*  fluence_cm2,
    const long*  material_no,
    const long*  RDD_model,
    const float*  RDD_parameters,
    const long*  ER_model,
    const float*  ER_parameters,
    const long*  gamma_model,
    const float*  gamma_parameters,
    const long*  N_runs,
    const long*   N2,
    const float*  fluence_factor,
    const bool*   write_output,
    const long*   nX,
    const float*  voxel_size_m,
    const bool*   lethal_events_mode,
    float*  results);

/**
* Computes HCP response and RE/RBE using Katz' Ion-Gamma-Kill approach
* according to Waligorski, 1988
*
* @param[in]  n      number of particle types in the mixed particle field (pointer to single variable)
* @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
* @param[in]  particle_no    type of the particles in the mixed particle field (pointer to array of size n)
* @see          AT_DataParticle.h for definition
* @param[in]  fluence_cm2    fluences for the given particles, doses in Gy if negative (pointer to array of size n)
* @param[in]  material_no    index number for detector material (pointer to single variable)
* @see          AT_DataMaterial.h for definition
* @param[in]  RDD_model    index number for chosen radial dose distribution (pointer to single variable)
* @param[in]  RDD_parameters    parameters for chosen radial dose distribution (pointer to array of size depending on chosen model)
* @see          AT_RDD.h for definition
* @param[in]  ER_model    index number for chosen electron-range model (pointer to single variable)
* @param[in]  ER_parameters    parameters for chosen electron-range model (pointer to array of size depending on chosen model)
* @see          AT_ElectronRange.h for definition
* @param[in]  gamma_model    index number for chosen gamma response (pointer to single variable)
* @param[in]  gamma_parameters  parameters for chosen gamma response (pointer to array of size depending on chosen model)
* @see          AT_GammaResponse.h for definition
* @param[in]  saturation_cross_section_factor  (algorithm specific)  scaling factor for the saturation cross section
* @see          Waligorski, 1988
* @param[out]  results      pointer to array of size 10 to be allocated by the user which will be used to return the results
*    results[0]    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D\n
*    results[1]    d_check         (algorithm independent)  not available with IGK\n
*    results[2]    S_HCP           (algorithm independent)  absolute particle response\n
*    results[3]    S_gamma         (algorithm independent)  absolute gamma response\n
*    results[4]    n_particles     (algorithm independent)  not available with IGK\n
*    results[5]    sI_cm2          (algorithm specific)     resulting ion saturation cross section in cm2\n
*    results[6]    gamma_dose_Gy   (algorithm specific)     dose contribution from gamma kills\n
*    results[7]    P_I             (algorithm specific)     ion kill probability\n
*    results[8]    P_G             (algorithm specific)     gamma kill probability\n
*    results[9]    not used        (algorithm specific)\n
* @return  none
*/
void AT_IGK(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const float*  fluence_cm2,
    const long*  material_no,
    const long*  RDD_model,
    const float*  RDD_parameters,
    const long*  ER_model,
    const float*  ER_parameters,
    const long*  gamma_model,
    const float*  gamma_parameters,
    const float*  saturation_cross_section_factor,
    float*  results);

/**
* Computes HCP response and RE/RBE using compound Poison process and
* statistical sampling (CPP_SS, the 'SPISS' algorithm)
*
* @param[in]  n      number of particle types in the mixed particle field (pointer to single variable)
* @param[in]  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
* @param[in]  particle_no    type of the particles in the mixed particle field (pointer to array of size n)
* @see          AT_DataParticle.h for definition
* @param[in]  fluence_cm2    fluences for the given particles, doses in Gy if negative (pointer to array of size n)
* @param[in]  material_no    index number for detector material (pointer to single variable)
* @see          AT_DataMaterial.h for definition
* @param[in]  RDD_model    index number for chosen radial dose distribution (pointer to single variable)
* @param[in]  RDD_parameters    parameters for chosen radial dose distribution (pointer to array of size depending on chosen model)
* @see          AT_RDD.h for definition
* @param[in]  ER_model    index number for chosen electron-range model (pointer to single variable)
* @param[in]  ER_parameters    parameters for chosen electron-range model (pointer to array of size depending on chosen model)
* @see          AT_ElectronRange.h for definition
* @param[in]  gamma_model    index number for chosen gamma response (pointer to single variable)
* @param[in]  gamma_parameters  parameters for chosen gamma response (pointer to array of size depending on chosen model)
* @see          AT_GammaResponse.h for definition
* @param[in] n_runs  (algorithm specific) number of points sampled for local dose distribution
* @param[in]  N2      (algorithm specific) number of bins per factor of two in local dose array (pointer to single variable)
* @param[in]  fluence_factor    factor to scale the fluences given as "fluence_cm2" with (pointer to single variable)
* @param[in]  write_output    if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory (pointer to single variable)
* @param[in]  importance_sampling  if unequal zero importance sampling will be applied to the single impact local dose distribution
* @param[in]  results      pointer to array of size 10 to be allocated by the user which will be used to return the results
*    results[0]    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D\n
*    results[1]    d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm\n
*    results[2]    S_HCP           (algorithm independent)  absolute particle response\n
*    results[3]    S_gamma         (algorithm independent)  absolute gamma response\n
*    results[4]    not used        (algorithm independent)\n
*    results[5]    not_used        (algorithm specific)\n
*    results[6]    not used        (algorithm specific)\n
*    results[7]    not_used        (algorithm specific)\n
*    results[8]    not used        (algorithm specific)\n
*    results[9]    not used        (algorithm specific)\n
*/
void AT_SPISS(  const long* n,
    const float* E_MeV_u,
    const long*  particle_no,
    const float* fluence_cm2,
    const long*  material_no,
    const long*  RDD_model,
    const float* RDD_parameters,
    const long*  ER_model,
    const float* ER_parameters,
    const long*  gamma_model,
    const float* gamma_parameters,
    const long*  n_runs, // TODO check sequence
    const long*  N2, // TODO check sequence
    const float* fluence_factor,
    const int*   write_output,
    const long*  importance_sampling,
    float*  results);

#endif /* AmTrack_H_ */
