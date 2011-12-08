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
void AT_GSM_sample_particle_positions( const long  number_of_field_components,
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
 * @param[in]  stopping_power_source_no                  TODO
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
void AT_GSM_dose_grid_from_particles_positions( const long  number_of_field_components,
    const double   E_MeV_u[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     stopping_power_source_no,
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
 * @param[in]   grid               grid
 * @param[in]   number_of_bins     number of bins for resulting histogram
 * @param[in]   bin_centers_Gy     midpoints of bins for resulting histogram (array of size number_of_bins)
 * @param[out]  zero_fraction      zero bin content - frequency of value zero
 * @param[out]  frequency          frequency vector for resulting histogram (array of size number_of_bins)
 * @length number_of_bins
 */
void AT_GSM_local_dose_distrib_from_dose_grid( const long     nX,
    const double** grid,
    const long     number_of_bins,
    const double   bin_centers_Gy[],
    double*        zero_fraction,
    double         frequency[]);


/**
 * Computes response distribution in grid
 *
 * @param[in]  nX                                        number of cells on grid side
 * @param[in]  gamma_model                               index number for chosen gamma response
 * @param[in]  gamma_parameters                          parameters for chosen gamma response (array of size GR_MAX_NUMBER_OF_PARAMETERS)
 * @see          AT_GammaResponse.h for definition
 * @param[in]  grid_D_Gy                                 grid dose pattern (array of size nX*nX)
 * @param[in]  lethal_events_mode                        if true, allows to do calculations for cell survival
 * @param[out] grid_response                             response pattern (array of size nX*nX)
 */
void AT_GSM_response_grid_from_dose_grid( const long      nX,
                const long      gamma_model,
                const double    gamma_parameters[],
                const double**  grid_D_Gy,
                const bool      lethal_events_mode,
                double**        grid_response);


/**
 * Computes HCP response and relative efficiency/RBE using summation of tracks
 * an a Cartesian grid (the 'GSM' algorithm).
 *
 * Be aware that this routine can take considerable time to compute depending on
 * the arguments, esp. for higher energy (>10 MeV/u) particles. It is therefore
 * advantageous to test your settings with a low number of runs first.
 *
 * @param[in]      number_of_field_components     number of components in the mixed particle field
 * @param[in]      E_MeV_u                        particle energy for each component in the mixed particle field [MeV/u] (array of size number_of_field_components)
 * @param[in]      particle_no                    particle type for each component in the mixed particle field (array of size number_of_field_components)
 * @see AT_DataParticle.h for definition
 * @param[in]      fluence_cm2_or_dose_Gy         if positive, particle fluence for each component in the mixed particle field [1/cm2]; if negative, particle dose for each component in the mixed particle field [Gy] (array of size number_of_field_components)
 * @param[in]      material_no                    index number for detector material
 * @param[in]      stopping_power_source_no       TODO
 * @see AT_DataMaterial.h for definition
 * @param[in]      rdd_model                      index number for chosen radial dose distribution
 * @param[in]      rdd_parameters                 parameters for chosen radial dose distribution (array of size 4)
 * @see AT_RDD.h for definition
 * @param[in]      er_model                       index number for chosen electron-range model
 * @see AT_ElectronRange.h for definition
 * @param[in]      gamma_model                    index number for chosen gamma response
 * @param[in]      gamma_parameters               parameters for chosen gamma response (array of size 9)
 * @see AT_GammaResponse.h for definition
 * @param[in]      N_runs                         number of runs within which track positions will be resampled
 * @param[in]      write_output                   if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory
 * @param[in]      nX                             number of voxels of the grid in x (and y as the grid is quadratic)
 * @param[in]      voxel_size_m                   side length of a voxel in m
 * @param[in]      lethal_events_mode             if true, allows to do calculations for cell survival
 * @param[out]     relative_efficiency            particle response at dose D / gamma response at dose D
 * @param[out]     d_check                        sanity check:  total dose (in Gy) as returned by the algorithm
 * @param[out]     S_HCP                          absolute particle response
 * @param[out]     S_gamma                        absolute gamma response
 * @param[out]     n_particles                    average number of particle tracks on the detector grid
 * @param[out]     sd_relative_efficiency         standard deviation for relative_efficiency
 * @param[out]     sd_d_check                     standard deviation for d_check
 * @param[out]     sd_S_HCP                       standard deviation for S_HCP
 * @param[out]     sd_S_gamma                     standard deviation for S_gamma
 * @param[out]     sd_n_particles                 standard deviation for n_particles
 */
void AT_run_GSM_method(  const long  number_of_field_components,
    const double   E_MeV_u[],
    const long     particle_no[],
    const double   fluence_cm2_or_dose_Gy[],
    const long     material_no,
    const long     stopping_power_source_no,
    const long     rdd_model,
    const double   rdd_parameters[],
    const long     er_model,
    const long     gamma_model,
    const double   gamma_parameters[],
    const long     N_runs,
    const bool     write_output,
    const long     nX,
    const double   voxel_size_m,
    const bool     lethal_events_mode,
    double*        relative_efficiency,
    double*        d_check,
    double*        S_HCP,
    double*        S_gamma,
    double*        n_particles,
    double*        sd_relative_efficiency,
    double*        sd_d_check,
    double*        sd_S_HCP,
    double*        sd_S_gamma,
    double*        sd_n_particles);

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
 * @param[in] stopping_power_source_no                   TODO
 * @param[in] nX                                         TODO
 * @param[in] pixel_size_m                               TODO
 * @param[in] number_of_bins                             TODO
 * @param[in] dose_bin_centers_Gy                        TODO
 * @param[out] random_number_generator_seed              TODO
 * @param[out] zero_dose_fraction                        TODO
 * @param[out] dose_frequency_Gy                         TODO
 */
void AT_GSM_local_dose_distrib( const long  number_of_field_components,
    const double   E_MeV_u[],
    const double   fluence_cm2[],
    const long     particle_no[],
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const long     stopping_power_source_no,
    const long     nX,
    const double   pixel_size_m,
    const long     number_of_bins,
    const double   dose_bin_centers_Gy[],
    unsigned long* random_number_generator_seed,
    double*        zero_dose_fraction,
    double         dose_frequency_Gy[]);

/**
 * TODO
 * @param[in] number_of_field_components
 * @param[in] E_MeV_u (array of size number_of_field_components)
 * @param[in] fluence_cm2 (array of size number_of_field_components)
 * @param[in] particle_no (array of size number_of_field_components)
 * @param[in] material_no
 * @param[in] rdd_model
 * @param[in] rdd_parameter (array of size 4)
 * @param[in] er_model
 * @param[in] stopping_power_source_no
 * @param[in] nX
 * @param[in] pixel_size_m
 * @param[in] N_runs
 * @param[in] N_repetitions
 * @param[in] number_of_bins
 * @param[in] dose_bin_centers_Gy (array of size number_of_bins)
 * @param[out] dose_bin_width_Gy (array of size number_of_bins)
 * @param[out] mean_d_check_Gy
 * @param[out] sd_d_check_Gy
 * @param[out] mean_zero_dose_fraction
 * @param[out] sd_zero_dose_fraction
 * @param[out] mean_dose_frequency_Gy (array of size number_of_bins)
 * @param[out] sd_dose_frequency_Gy (array of size number_of_bins)
 */
void AT_GSM_multiple_local_dose_distrib( const long  number_of_field_components,
    const double   	E_MeV_u[],
    const double   	fluence_cm2[],
    const long     	particle_no[],
    const long     	material_no,
    const long     	rdd_model,
    const double   	rdd_parameter[],
    const long     	er_model,
    const long      stopping_power_source_no,
    const long     	nX,
    const double   	pixel_size_m,
    const long		N_runs,
    const long		N_repetitions,
    const long     	number_of_bins,
    const double   	dose_bin_centers_Gy[],
    double    		dose_bin_width_Gy[],
    double*       	mean_d_check_Gy,
    double*       	sd_d_check_Gy,
    double*       	mean_zero_dose_fraction,
    double*       	sd_zero_dose_fraction,
    double        	mean_dose_frequency_Gy[],
    double        	sd_dose_frequency_Gy[]);


#endif /* AT_ALGORITHMS_GSM_H_ */
