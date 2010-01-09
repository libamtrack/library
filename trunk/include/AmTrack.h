#ifndef AmTrack_H_
#define AmTrack_H_

/**
 *    AmTrack.h
 *    =========
 *
 *    Created on: 28.07.2009
 *    Author: greilich
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
#include "AT_Data.h"
#include "AT_Utils.h"
#include "AT_Functions.h"
#include "AT_RDD.h"
#include "AT_SuccessiveConvolutions.h"
#include "AT_GammaResponse.h"
#include "AT_FileOperations.h"
#include "AT_NumericalRoutines.h"
#include "AT_Transport.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>

extern int indent_counter;
extern char isp[];
extern FILE * debf;

/**
 * Computes the efficiency using the SPIFF algorithm
 *
 * @param  n      number of particle types in the mixed particle field (pointer to single variable)
 * @param  E_MeV_u      energy of particles in the mixed particle field (pointer to array of size n)
 * @param  particle_no    type of the particles in the mixed particle field (pointer to array of size n)
 * @see          AT_Data.h for definition
 * @param  fluence_cm2    fluences for the given particles, doses in Gy if negative (pointer to array of size n)
 * @param  material_no    index number for detector material (pointer to single variable)
 * @see          AT_Constants.h for definition
 * @param  RDD_model    index number for chosen radial dose distribution (pointer to single variable)
 * @param  RDD_parameters    parameters for chosen radial dose distribution (pointer to array of size depending on chosen model)
 * @see          AT_Constants.h for definition
 * @param  ER_model    index number for chosen electron-range model (pointer to single variable)
 * @param  ER_parameters    parameters for chosen electron-range model (pointer to array of size depending on chosen model)
 * @see          AT_Constants.h for definition
 * @param  gamma_model    index number for chosen gamma response (pointer to single variable)
 * @param  gamma_parameters  parameters for chosen gamma response (pointer to array of size depending on chosen model)
 * @see          AT_Constants.h for definition
 * @param  N2      (algorithm specific) number of bins per factor of two in local dose array (pointer to single variable)
 * @param  fluence_factor    factor to scale the fluences given as "fluence_cm2" with (pointer to single variable)
 * @param  write_output    if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory (pointer to single variable)
 * @param  shrink_tails    (algorithm specific) if true, tails of the local dose distribution, contributing less than "shrink_tails_under" are cut (pointer to single variable)
 * @param  shrink_tails_under  (algorithm specific) limit for tail cutting in local dose distribution (pointer to single variable)
 * @param  adjust_N2    (algorithm specific) if true, "N2" will be increase if necessary at high fluence to ensure sufficient binning resolution
 * @param  results      pointer to array of size 10 to be allocated by the user which will be used to return the results
 *     results[0]    efficiency        (algorithm independent) main result:   particle response at dose D / gamma response at dose D
 *    results[1]    d_check          (algorithm independent) sanity check:  total dose (in Gy) as returned by the alogrithm
 *    results[2]    S_HCP          (algorithm independent)         absolute particle response
 *     results[3]    S_gamma          (algorithm independent)          absolute gamma response
 *     results[4]    not used        (algorithm independent)
 *    results[5]    u            (algorithm specific)          mean number of tracks contributing to representative point
 *    results[6]    u_start          (algorithm specific)          low starting value for mean number of tracks, where linearisation is applied
 *    results[7]    n_convolutions      (algorithm specific)          number of convolutions performed
 *     results[8]    not used        (algorithm specific)
 *     results[9]    not used        (algorithm specific)
 * @return  none
 */

void AT_SPIFF(  long*  n,
    float*  E_MeV_u,
    long*  particle_no,
    float*  fluence_cm2,
    long*  material_no,
    long*  RDD_model,
    float*  RDD_parameters,
    long*  ER_model,
    float*  ER_parameters,
    long*  gamma_model,
    float*  gamma_parameters,
    long*  N2,
    float*  fluence_factor,
    int*  write_output,
    int*  shrink_tails,
    float*  shrink_tails_under,
    int*  adjust_N2,
    bool*   lethal_events_mode,
    float*  results);

void AT_interparticleDistance_m(       long*   n,
                float*  LET_MeV_cm2_g,
                float*  fluence_cm2,
                float*  results_m
);

void AT_inv_interparticleDistance_Gy(  long*   n,
                float*  LET_MeV_cm2_g,
                float*  distance_m,
                float*  results_Gy
);

void AT_inv_interparticleDistance_cm2( long*   n,
                float*  distance_m,
                float*  results_cm2
);

void AT_GSM(  long*  n,
    float*  E_MeV_u,
    long*  particle_no,
    float*  fluence_cm2,
    long*  material_no,
    long*  RDD_model,
    float*  RDD_parameters,
    long*  ER_model,
    float*  ER_parameters,
    long*  gamma_model,
    float*  gamma_parameters,
    long*  N_runs,
    long*  N2,
    float*  fluence_factor,
    bool*  write_output,
    long*  nX,
    float*  grid_size_m,
    bool*  lethal_events_mode,
    float*  results);

void AT_IGK(  long*  n,
    float*  E_MeV_u,
    long*  particle_no,
    float*  fluence_cm2,
    long*  material_no,
    long*  RDD_model,
    float*  RDD_parameters,
    long*  ER_model,
    float*  ER_parameters,
    long*  gamma_model,
    float*  gamma_parameters,
    float*  results);

void AT_SPISS(	long*  	n,
		float*  E_MeV_u,
		long*  	particle_no,
		float*  fluence_cm2,
		long*  	material_no,
		long*  	RDD_model,
		float*  RDD_parameters,
		long*  	ER_model,
		float*  ER_parameters,
		long*  	gamma_model,
		float*  gamma_parameters,
		long*  	N2,
		long*	n_runs,
		float*  fluence_factor,
		int*  	write_output,
		long*  	importance_sampling,
		float*  results);

#endif /* AmTrack_H_ */
