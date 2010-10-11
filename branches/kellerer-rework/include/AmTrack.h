#ifndef AmTrack_H_
#define AmTrack_H_

/**
 * @brief libamtrack main file holding the amorphous track routines for RE/RBE calculation
 */

/*
 *    AmTrack.h
 *    =========
 *
 *    Created on: 28.07.2009
 *    Creator: greilich
 *
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
#include <time.h>

#include "AT_Constants.h"
#include "AT_RDD.h"
#include "AT_RDD_ExtendedTarget.h"
#include "AT_SuccessiveConvolutions.h"
#include "AT_GammaResponse.h"
#include "AT_Histograms.h"
#include "AT_PhysicsRoutines.h"
#include "AT_NumericalRoutines.h"
#include "AT_KatzModel.h"
#include "AT_Algorithms_GSM.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>


/**
 * Computes HCP response and RE/RBE using compound Poison process and
 * successive convolutions (CPP_SC, the 'SPIFF' algorithm)
 *
 * @param[in]  n                   number of particle types in the mixed particle field
 * @param[in]  E_MeV_u             energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no         type of the particles in the mixed particle field (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2_or_dose_Gy         fluences for the given particles, doses in Gy if negative (array of size n)
 * @param[in]  material_no         index number for detector material
 * @see          AT_DataMaterial.h for definition
 * @param[in]  rdd_model           index number for chosen radial dose distribution
 * @param[in]  rdd_parameters      parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @see          AT_RDD.h for definition
 * @param[in]  er_model            index number for chosen electron-range model
 * @see          AT_ElectronRange.h for definition
 * @param[in]  gamma_model         index number for chosen gamma response
 * @param[in]  gamma_parameters    parameters for chosen gamma response (array of size depending on chosen model)
 * @see          AT_GammaResponse.h for definition
 * @param[in,out]  N2      (algorithm specific) number of bins per factor of two in local dose array
 * @param[in]  fluence_factor      factor to scale the fluences given as "fluence_cm2" with
 * @param[in]  write_output        if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory
 * @param[in]  shrink_tails        (algorithm specific) if true, tails of the local dose distribution, contributing less than "shrink_tails_under" are cut
 * @param[in]  shrink_tails_under  (algorithm specific) limit for tail cutting in local dose distribution
 * @param[in]  adjust_N2           (algorithm specific) if true, "N2" will be increase if necessary at high fluence to ensure sufficient binning resolution
 * @param[in]  lethal_events_mode  (algorithm specific) if true, allows to do calculations for cell survival
 * @param[out]  results            array of size 10 to be allocated by the user which will be used to return the results\n
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
void AT_run_SPIFF_method(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameters[],
    const long    er_model,
    const long    gamma_model,
    const double  gamma_parameters[],
    long          N2, // TODO investigate if this can be changed inside
    const double  fluence_factor,
    const bool    write_output,
    const bool    shrink_tails,
    const double  shrink_tails_under,
    const bool    adjust_N2,
    const bool    lethal_events_mode,
    double        results[]);


/**
 * Computes HCP response and RE/RBE using Katz' Ion-Gamma-Kill approach
 * according to Waligorski, 1988
 *
 * @param[in]  n                   number of particle types in the mixed particle field
 * @param[in]  E_MeV_u             energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no         type of the particles in the mixed particle field (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2_or_dose_Gy         fluences for the given particles, doses in Gy if negative (array of size n)
 * @param[in]  material_no         index number for detector material
 * @see          AT_DataMaterial.h for definition
 * @param[in]  RDD_model           index number for chosen radial dose distribution
 * @param[in]  RDD_parameters      parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @see          AT_RDD.h for definition
 * @param[in]  ER_model            index number for chosen electron-range model
 * @see          AT_ElectronRange.h for definition
 * @param[in]  gamma_model         index number for chosen gamma response
 * @param[in]  gamma_parameters    parameters for chosen gamma response (array of size depending on chosen model)
 * @see          AT_GammaResponse.h for definition
 * @param[in]  saturation_cross_section_factor  (algorithm specific)  scaling factor for the saturation cross section
 * @see          Waligorski, 1988
 * @param[in]  write_output        if true, a protocol is written to a file in the working directory
 * @param[out]  results            array of size 10 to be allocated by the user which will be used to return the results\n
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
void AT_run_IGK_method(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    RDD_model,
    const double  RDD_parameters[],
    const long    ER_model,
    const long    gamma_model,
    const double  gamma_parameters[],
    const double  saturation_cross_section_factor,
    const bool    write_output,
    double  results[]);


/**
 * Computes HCP response and RE/RBE using compound Poison process and
 * statistical sampling (CPP_SS, the 'SPISS' algorithm)
 *
 * @param[in]  n                   number of particle types in the mixed particle field
 * @param[in]  E_MeV_u             energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no         type of the particles in the mixed particle field (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2_or_dose_Gy         fluences for the given particles, doses in Gy if negative (array of size n)
 * @param[in]  material_no         index number for detector material
 * @see          AT_DataMaterial.h for definition
 * @param[in]  RDD_model           index number for chosen radial dose distribution
 * @param[in]  RDD_parameters      parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @see          AT_RDD.h for definition
 * @param[in]  ER_model            index number for chosen electron-range model
 * @see          AT_ElectronRange.h for definition
 * @param[in]  gamma_model         index number for chosen gamma response
 * @param[in]  gamma_parameters    parameters for chosen gamma response (array of size depending on chosen model)
 * @see          AT_GammaResponse.h for definition
 * @param[in] n_runs               (algorithm specific) number of points sampled for local dose distribution
 * @param[in]  N2                  (algorithm specific) number of bins per factor of two in local dose array
 * @param[in]  fluence_factor      factor to scale the fluences given as "fluence_cm2" with
 * @param[in]  write_output        if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory
 * @param[in]  importance_sampling if unequal zero importance sampling will be applied to the single impact local dose distribution
 * @param[in]  results             array of size 10 to be allocated by the user which will be used to return the results\n
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
void AT_run_SPISS_method(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    RDD_model,
    const double  RDD_parameters[],
    const long    ER_model,
    const long    gamma_model,
    const double  gamma_parameters[],
    const long    n_runs,
    const long    N2,
    const double  fluence_factor,
    const int     write_output,
    const long    importance_sampling,
    double        results[]);


/**
 * TODO
 * @param number_of_field_components
 * @param E_MeV_u
 * @param fluence_cm2
 * @param particle_no
 * @param material_no
 * @param rdd_model
 * @param rdd_parameter
 * @param er_model
 * @param nX
 * @param pixel_size_m
 * @param N_runs
 * @param N_repetitions
 * @param number_of_bins
 * @param dose_bin_centers_Gy
 * @param dose_bin_width_Gy
 * @param mean_d_check_Gy
 * @param sd_d_check_Gy
 * @param mean_zero_dose_fraction
 * @param sd_zero_dose_fraction
 * @param mean_dose_frequency_Gy
 * @param sd_dose_frequency_Gy
 */
void AT_GSM_calculate_multiple_dose_histograms( const long  number_of_field_components,
    const double   	E_MeV_u[],
    const double   	fluence_cm2[],
    const long     	particle_no[],
    const long     	material_no,
    const long     	rdd_model,
    const double   	rdd_parameter[],
    const long     	er_model,
    const long     	nX,
    const double   	pixel_size_m,
    const long		N_runs,
    const long		N_repetitions,
    const long     	number_of_bins,
    const double   	dose_bin_centers_Gy[],
    double    		dose_bin_width_Gy[],
    double *       	mean_d_check_Gy,
    double *       	sd_d_check_Gy,
    double *       	mean_zero_dose_fraction,
    double *       	sd_zero_dose_fraction,
    double        	mean_dose_frequency_Gy[],
    double        	sd_dose_frequency_Gy[]);

#endif /* AmTrack_H_ */
