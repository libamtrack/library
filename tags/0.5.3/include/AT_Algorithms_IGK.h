#ifndef AT_ALGORITHMS_IGK_H_
#define AT_ALGORITHMS_IGK_H_

/**
 * @brief Ion-gamma-kill model
 */

/*
 *    AT_Algorithms_IGK.h
 *    ===================
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
#include "AT_DataStoppingPower.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>


/**
 * Computes HCP response and relative efficiency/RBE using Katz' Ion-Gamma-Kill approach
 * according to Waligorski, 1988
 *
 * @param[in]      number_of_field_components       number of components in the mixed particle field
 * @param[in]      E_MeV_u                          particle energy for each component in the mixed particle field [MeV/u] (array of size number_of_field_components)
 * @param[in]      particle_no                      particle type for each component in the mixed particle field (array of size number_of_field_components)
 * @see AT_DataParticle.h for definition
 * @param[in]      fluence_cm2_or_dose_Gy           if positive, particle fluence for each component in the mixed particle field [1/cm2]; if negative, particle dose for each component in the mixed particle field [Gy] (array of size number_of_field_components)
 * @param[in]      material_no                      index number for detector material
 * @param[in]      stopping_power_source_no         stopping power source number (PSTAR,...)
 * @see AT_DataMaterial.h for definition
 * @param[in]      rdd_model                        index number for chosen radial dose distribution
 * @param[in]      rdd_parameters                   parameters for chosen radial dose distribution (array of size 4)
 * @see AT_RDD.h for definition
 * @param[in]      er_model                         index number for chosen electron-range model
 * @see AT_ElectronRange.h for definition
 * @param[in]      gamma_model                      index number for chosen gamma response
 * @param[in]      gamma_parameters                 parameters for chosen gamma response (array of size 9)
 * @see AT_GammaResponse.h for definition
 * @param[in]      saturation_cross_section_factor  scaling factor for the saturation cross section
 * @see Waligorski, 1988
 * @param[in]      write_output                     if true, a protocol is written to a file in the working directory
 * @param[out]     relative_efficiency              particle response at dose D / gamma response at dose D
 * @param[out]     S_HCP                            absolute particle response
 * @param[out]     S_gamma                          absolute gamma response
 * @param[out]     sI_cm2                           resulting ion saturation cross section in cm2
 * @param[out]     gamma_dose_Gy                    dose contribution from gamma kills
 * @param[out]     P_I                              ion kill probability
 * @param[out]     P_g                              gamma kill probability
 */
void AT_run_IGK_method(  const long  number_of_field_components,
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
    const double  saturation_cross_section_factor,
    const bool    write_output,
    double*       relative_efficiency,
    double*       S_HCP,
    double*       S_gamma,
    double*       sI_cm2,
    double*       gamma_dose_Gy,
    double*       P_I,
    double*       P_g);





#endif /* AT_ALGORITHMS_IGK_H_ */
