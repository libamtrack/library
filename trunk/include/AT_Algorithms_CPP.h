#ifndef AT_ALGORITHMS_CPP_H_
#define AT_ALGORITHMS_CPP_H_

/**
 * @brief Algorithms for ATMs based on Compound Poisson Processes (CPP)
 */

/*
 *    AT_Algorithms.h
 *    ===============
 *
 *    Created on: 18.11.2010
 *    Creator: greilich
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
 * Computes HCP response and relative efficiency/RBE using compound Poison processes and
 * successive convolutions (CPP_SC, the 'SPIFF' algorithm)
 *
 * @param[in]      number_of_field_components     number of components in the mixed particle field
 * @param[in]      E_MeV_u                        particle energy for each component in the mixed particle field [MeV/u] (array of size number_of_field_components)
 * @param[in]      particle_no                    particle type for each component in the mixed particle field (array of size number_of_field_components)
 * @see AT_DataParticle.h for definition
 * @param[in]      fluence_cm2_or_dose_Gy         if positive, particle fluence for each component in the mixed particle field [1/cm2]; if negative, particle dose for each component in the mixed particle field [Gy] (array of size number_of_field_components)
 * @param[in]      material_no                    index number for detector material
 * @see AT_DataMaterial.h for definition
 * @param[in]      stopping_power_source_no       TODO
 * @param[in]      rdd_model                      index number for chosen radial dose distribution
 * @param[in]      rdd_parameters                 parameters for chosen radial dose distribution (array of size 4)
 * @see AT_RDD.h for definition
 * @param[in]      er_model                       index number for chosen electron-range model
 * @see AT_ElectronRange.h for definition
 * @param[in]      gamma_model                    index number for chosen gamma response
 * @param[in]      gamma_parameters               parameters for chosen gamma response (array of size 9)
 * @see AT_GammaResponse.h for definition
 * @param[in,out]  N2                             number of bins per factor of two for the dose scale of local dose histogram
 * @param[in]      fluence_factor                 factor to scale the fluences / doses given in "fluence_cm2_or_dose_Gy" with
 * @param[in]      write_output                   if true, a log-file is written to "SuccessiveConvolutions.txt" in the working directory
 * @param[in]      shrink_tails                   if true, tails of the local dose distribution, contributing less than "shrink_tails_under" are cut
 * @param[in]      shrink_tails_under             limit for tail cutting in local dose distribution
 * @param[in]      adjust_N2                      if true, "N2" will be increase if necessary at high fluence to ensure sufficient local dose histogram resolution
 * @param[in]      lethal_events_mode             if true, computations are done for dependent subtargets
 * @param[out]     relative_efficiency            particle response at dose D / gamma response at dose D
 * @param[out]     d_check                        sanity check:  total dose (in Gy) as returned by the algorithm
 * @param[out]     S_HCP                          absolute particle response
 * @param[out]     S_gamma                        absolute gamma response
 * @param[out]     mean_number_of_tracks_contrib  mean number of tracks contributing to representative point
 * @param[out]     start_number_of_tracks_contrib low fluence approximation for mean number of tracks contributing to representative point (start value for successive convolutions)
 * @param[out]     n_convolutions                 number of convolutions performed to reach requested dose/fluence
 * @param[out]     lower_Jensen_bound             lower bound for Jensen's inequity
 * @param[out]     upper_Jensen_bound             upper bound for Jensen's inequity
 */
void AT_run_CPPSC_method(  const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    stopping_power_source_no,
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
    double*       relative_efficiency,
    double*       d_check,
    double*       S_HCP,
    double*       S_gamma,
    double*       mean_number_of_tracks_contrib,
    double*       start_number_of_tracks_contrib,
    long*         n_convolutions,
    double*		  lower_Jensen_bound,
    double*       upper_Jensen_bound);

/**
 * Computes HCP response and RE/RBE using compound Poison process and
 * statistical sampling (CPP_SS, the 'SPISS' algorithm)
 *
 * @param[in]  number_of_field_components  number of particle types in the mixed particle field
 * @param[in]  E_MeV_u             energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no         type of the particles in the mixed particle field (array of size number_of_field_components)
 * @see          AT_DataParticle.h for definition
 * @param[in]  fluence_cm2_or_dose_Gy         fluences for the given particles, doses in Gy if negative (array of size number_of_field_components)
 * @param[in]  material_no         index number for detector material
 * @see          AT_DataMaterial.h for definition
 * @param[in]  stopping_power_source_no         TODO
 * @param[in]  RDD_model           index number for chosen radial dose distribution
 * @param[in]  RDD_parameters      parameters for chosen radial dose distribution (array of size 4)
 * @see          AT_RDD.h for definition
 * @param[in]  ER_model            index number for chosen electron-range model
 * @see          AT_ElectronRange.h for definition
 * @param[in]  gamma_model         index number for chosen gamma response
 * @param[in]  gamma_parameters    parameters for chosen gamma response (array of size 9)
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
void AT_run_CPPSS_method(  const long  number_of_field_components,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    stopping_power_source_no,
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

#endif /* AT_ALGORITHMS_CPP_H_ */
