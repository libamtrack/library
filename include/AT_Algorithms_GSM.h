#ifndef AT_ALGORITHMS_GSM_H_
#define AT_ALGORITHMS_GSM_H_

/**
 * @brief GSM algorithm implementation
 */

/*
 *    AT_Algorithms_GSM.h
 *    =========
 *
 *    Created on: 30.09.2010
 *    Creator: kongruencja
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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>


/**
 * Computes how many particles will be in grid area (using Poissonian distribution
 * and fluence) and distribute particles randomly on grid.
 * Particles positions will be saved in x_position and y_position tables.
 * @param[in]   number_of_field_components               number of particle types in the mixed particle field
 * @param[in]   fluence_cm2                              fluences for the given particles (array of size number_of_field_components)
 * @param[in]   sample_grid_size_m                       linear size of grid (length of grid square side)
 * @param[in]   random_number_generator_seed             GSL random number generator seed
 * @param[out]  number_of_particles_in_field_component   table of number of particles of each type in the mixed particle field
 * @param[out]  x_position                               two-dimensional table (1st dimension - component of mixed field,
 *                                                         2nd dimension - particle index for given component)
 *                                                         of X-positions of particles of each type in the mixed particle field
 * @param[out]  y_position                               two-dimensional table (1st dimension - component of mixed field,
 *                                                         2nd dimension - particle index for given component)
 *                                                         of Y-positions of particles of each type in the mixed particle field
 */
void AT_GSM_shoot_particles_on_grid( const long  number_of_field_components,
    const double         fluence_cm2[],
    const double         sample_grid_size_m,
    unsigned long* 		 random_number_generator_seed,
    long                 number_of_particles_in_field_component[],
    double*              x_position[],
    double*              y_position[]);


/**
 * Computes dose distribution in grid
 *
 * @param[in]  number_of_field_components                number of particle types in the mixed particle field
 * @param[in]  E_MeV_u                                   energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  particle_no                               type of the particles in the mixed particle field (array of size number_of_field_components)
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no                               index number for detector material
 * @see          AT_DataMaterial.h for definition
 * @param[in]  rdd_model                                 index number for chosen radial dose distribution
 * @param[in]  rdd_parameter                             parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @see          AT_RDD.h for definition
 * @param[in]  er_model                                  index number for chosen electron-range model
 * @see          AT_ElectronRange.h for definition
 * @param[in]  number_of_particles_in_field_component    table of number of particles of each type in the mixed particle field
 * @param[in]  x_position                                two-dimensional table (1st dimension - component of mixed field,
 *                                                         2nd dimension - particle index for given component)
 *                                                         of X-positions of particles of each type in the mixed particle field
 * @param[in]  y_position                                two-dimensional table (1st dimension - component of mixed field,
 *                                                         2nd dimension - particle index for given component)
 *                                                         of Y-positions of particles of each type in the mixed particle field
 * @param[in]  nX                                        number of cells on grid side
 * @param[in]  pixel_size_m                              linear size of pixel in grid
 * @param[out] grid_D_Gy                                 grid dose pattern (2-D array of dimensions nX x nX, linearly allocated)
 */
void AT_GSM_calculate_dose_pattern( const long  number_of_field_components,
    const double   E_MeV_u[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     number_of_particles_in_field_component[],
    const double*  x_position[],
    const double*  y_position[],
    const long     nX,
    const double   pixel_size_m,
    double**       grid_D_Gy);


/**
 * Computes histogram from given grid
 *
 * @param[in]   nX                 number of cells on one grid side
 * @param[in]   number_of_bins     number of bins for resulting histogram
 * @param[in]   bin_centers_Gy     midpoints of bins for resulting histogram
 * @param[out]  zero_fraction      zero bin content - frequency of value zero
 * @param[out]  frequency          frequency vector for resulting histogram
 * @length number_of_bins
 */
void AT_GSM_calculate_histogram_from_grid( const long     nX,
    const double** grid,
    const long     number_of_bins,
    const double   bin_centers_Gy[],
    double *       zero_fraction,
    double         frequency[]);


/**
 * Computes local dose histogram for given mixed field
 *
 * @param[in]  number_of_field_components                number of particle types in the mixed particle field
 * @param[in]  E_MeV_u                                   energy of particles in the mixed particle field (array of size number_of_field_components)
 * @param[in]  fluence_cm2                               fluences for the given particles (array of size number_of_field_components)
 * @param[in]  particle_no                               type of the particles in the mixed particle field (array of size number_of_field_components)
 * @see          AT_DataParticle.h for definition
 * @param[in]  material_no                               index number for detector material
 * @see          AT_DataMaterial.h for definition
 * @param[in]  rdd_model                                 index number for chosen radial dose distribution
 * @param[in]  rdd_parameter                             parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @see          AT_RDD.h for definition
 * @param[in]  er_model                                  index number for chosen electron-range model
 * @see          AT_ElectronRange.h for definition
 */
void AT_GSM_calculate_dose_histogram( const long  number_of_field_components,
    const double   E_MeV_u[],
    const double   fluence_cm2[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     nX,
    const double   pixel_size_m,
    const long     number_of_bins,
    const double   dose_bin_centers_Gy[],
    unsigned long* random_number_generator_seed,
    double *       zero_dose_fraction,
    double         dose_frequency_Gy[]);


/**
 * Computes response distribution in grid
 *
 * @param[in]  nX                                        number of cells on grid side
 * @param[in]  gamma_model                               index number for chosen gamma response
 * @param[in]  gamma_parameters                          parameters for chosen gamma response (array of size depending on chosen model)
 * @see          AT_GammaResponse.h for definition
 * @param[in]  grid_D_Gy                                 grid dose pattern (2-D array of dimensions nX x nX)
 * @param[in]  lethal_events_mode                        if true, allows to do calculations for cell survival
 * @param[out] grid_response                             response pattern (2-D array of dimensions nX x nX)
 */
void AT_GSM_calculate_local_response_grid( const long      nX,
                const long      gamma_model,
                const double    gamma_parameters[],
                const double**  grid_D_Gy,
                const bool      lethal_events_mode,
                double**        grid_response);

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

/**
 * Computes HCP response and relative efficiency/RBE using summation of tracks
 * an a Cartesian grid (the 'GSM' algorithm)
 *
 * @param[in]      n                              number of components in the mixed particle field
 * @param[in]      E_MeV_u                        particle energy for each component in the mixed particle field [MeV/u] (array of size n)
 * @param[in]      particle_no                    particle type for each component in the mixed particle field (array of size n)
 * @see AT_DataParticle.h for definition
 * @param[in]      fluence_cm2_or_dose_Gy         if positive, particle fluence for each component in the mixed particle field [1/cm2]; if negative, particle dose for each component in the mixed particle field [Gy] (array of size n)
 * @param[in]      material_no                    index number for detector material
 * @see AT_DataMaterial.h for definition
 * @param[in]      rdd_model                      index number for chosen radial dose distribution
 * @param[in]      rdd_parameters                 parameters for chosen radial dose distribution (array of size 4)
 * @see AT_RDD.h for definition
 * @param[in]      er_model                       index number for chosen electron-range model
 * @see AT_ElectronRange.h for definition
 * @param[in]      gamma_model                    index number for chosen gamma response
 * @param[in]      gamma_parameters               parameters for chosen gamma response (array of size 4)
 * @see AT_GammaResponse.h for definition
 * @param[in]      N_runs                         number of runs within which track positions will be resampled
 * @param[in]      write_output                   if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory
 * @param[in]      nX                             number of voxels of the grid in x (and y as the grid is quadratic)
 * @param[in]      voxel_size_m                   side length of a voxel in m
 * @param[in]      lethal_events_mode             if true, allows to do calculations for cell survival
 * @param[out]     results                        results (array of size 10)
 *    results[0]    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D\n
 *    results[1]    d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm\n
 *    results[2]    S_HCP           (algorithm independent)  absolute particle response\n
 *    results[3]    S_gamma         (algorithm independent)  absolute gamma response\n
 *    results[4]    n_particles     (algorithm independent)  average number of particle tracks on the detector grid\n
 *    results[5]    sd_efficiency   (algorithm specific)     standard deviation for results[0]\n
 *    results[6]    sd_d_check      (algorithm specific)     standard deviation for results[1]\n
 *    results[7]    sd_S_HCP        (algorithm specific)     standard deviation for results[2]\n
 *    results[8]    sd_S_gamma      (algorithm specific)     standard deviation for results[3]\n
 *    results[9]    sd_n_particles  (algorithm specific)     standard deviation for results[4]
 */
void AT_run_GSM_method(  const long  n,
    const double   E_MeV_u[],
    const long     particle_no[],
    const double   fluence_cm2_or_dose_Gy[],
    const long     material_no,
    const long     RDD_model,
    const double   RDD_parameters[],
    const long     ER_model,
    const long     gamma_model,
    const double   gamma_parameters[],
    const long     N_runs,
    const bool     write_output,
    const long     nX,
    const double   voxel_size_m,
    const bool     lethal_events_mode,
    double  results[]);


#endif /* AT_ALGORITHMS_GSM_H_ */
