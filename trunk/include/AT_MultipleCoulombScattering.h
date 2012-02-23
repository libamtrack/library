#ifndef AT_MULTIPLECOULOMBSCATTERING_H_
#define AT_MULTIPLECOULOMBSCATTERING_H_

/**
 * @brief Moliere theory of multiple Coulomb scattering
 */

/*
 *    AT_MultipleCoulombScattering.h
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

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "AT_Constants.h"
#include "AT_DataMaterial.h"
#include "AT_DataParticle.h"
#include "AT_PhysicsRoutines.h"


/**
 *  Returns characteristic single scattering angle chi_c
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_charge_e    	charge number of particle
 * @param[in]  target_thickness_cm      thickness of the target material in cm
 * @param[in]  element_acronym		    elemental symbol of target material (array of size PARTICLE_NAME_NCHAR)
 * @return     chi_c in rad
 */
double AT_characteristic_single_scattering_angle_single( const double E_MeV_u,
														 const int    particle_charge_e,
														 const double target_thickness_cm,
														 const char   element_acronym[]);


/**
 *  Returns characteristic single scattering angles chi_c
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_charge_e    	vector of charge numbers of particle (array of size n)
 * @param[in]  target_thickness_cm      vector of thicknesses of target material in cm (array of size n)
 * @param[in]  element_acronym		    vector of elemental symbols of target material (array of size n)
 * @param[out] chi_c                    vector of characteristic single scattering angles in rad (array of size n)
 * @return     status code
 */
int AT_characteristic_single_scattering_angle( const long  		n,
											   const double 	E_MeV_u[],
											   const int		particle_charge_e[],
											   const double		target_thickness_cm[],
											   char*			element_acronym[],
											   double        	chi_c[]);


/**
 *  Returns screening angle chi_a
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_charge_e    	charge number of particle
 * @param[in]  element_acronym		    elemental symbol of target material (array of size PARTICLE_NAME_NCHAR)
 * @return     chi_a in rad
 */
double AT_screening_angle_single( const double E_MeV_u,
							      const int    particle_charge_e,
							      const char   element_acronym[]);


/**
 *  Returns screening angles chi_a
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_charge_e    	vector of charge numbers of particle (array of size n)
 * @param[in]  element_acronym		    vector of elemental symbols of target material (array of size n)
 * @param[out] chi_a                    vector of screening angles in rad (array of size n)
 * @return     status code
 */
int AT_screening_angle( const long 		 n,
						const double 	E_MeV_u[],
						const int		particle_charge_e[],
    					char*			element_acronym[],
    					double        	chi_a[]);


/**
 *  Returns effective collision number exp_b
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_charge_e    	charge number of particle
 * @param[in]  target_thickness_cm      thickness of the target material in cm
 * @param[in]  element_acronym		    elemental symbol of target material (array of size PARTICLE_NAME_NCHAR)
 * @return     exp_b
 */
double AT_effective_collision_number_single( const double E_MeV_u,
											 const int    particle_charge_e,
											 const double target_thickness_cm,
											 const char   element_acronym[]);


/**
 *  Returns effective collision numbers exp_b
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_charge_e    	vector of charge numbers of particle (array of size n)
 * @param[in]  target_thickness_cm      vector of thicknesses of target material in cm (array of size n)
 * @param[in]  element_acronym		    vector of elemental symbols of target material (array of size n)
 * @param[out] exp_b	                vector of effective collision numbers (array of size n)
 * @return     status code
 */
int AT_effective_collision_number( const long  		n,
								   const double 	E_MeV_u[],
								   const int		particle_charge_e[],
								   const double		target_thickness_cm[],
								   char*			element_acronym[],
								   double        	exp_b[]);


/**
 *  Returns reduced target thickness B
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_charge_e    	charge number of particle
 * @param[in]  target_thickness_cm      thickness of the target material in cm
 * @param[in]  element_acronym		    elemental symbol of target material (array of size PARTICLE_NAME_NCHAR)
 * @return     B
 */
double AT_reduced_target_thickness_single( const double E_MeV_u,
								           const int    particle_charge_e,
								           const double target_thickness_cm,
								           const char   element_acronym[]);


/**
 *  Returns reduced target thicknesses B
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_charge_e    	vector of charge numbers of particle (array of size n)
 * @param[in]  target_thickness_cm      vector of thicknesses of target material in cm (array of size n)
 * @param[in]  element_acronym		    vector of elemental symbols of target material (array of size n)
 * @param[out] B		                vector of reduced target thicknesses (array of size n)
 * @return     status code
 */
int AT_reduced_target_thickness( const long  	n,
								 const double 	E_MeV_u[],
								 const int		particle_charge_e[],
								 const double	target_thickness_cm[],
								 char*			element_acronym[],
								 double        	B[]);


/**
 *  Returns characteristic multiple scattering angle Theta_M
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_charge_e    	charge number of particle
 * @param[in]  target_thickness_cm      thickness of the target material in cm
 * @param[in]  element_acronym		    elemental symbol of target material (array of size PARTICLE_NAME_NCHAR)
 * @return     Theta_M in rad
 */
double AT_characteristic_multiple_scattering_angle_single( const double E_MeV_u,
														   const int    particle_charge_e,
														   const double target_thickness_cm,
														   const char   element_acronym[]);


/**
 *  Returns characteristic multiple scattering angles Theta_M
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_charge_e    	vector of charge numbers of particle (array of size n)
 * @param[in]  target_thickness_cm      vector of thicknesses of target material in cm (array of size n)
 * @param[in]  element_acronym		    vector of elemental symbols of target material (array of size n)
 * @param[out] Theta_M		            vector of characteristic multiple scattering angles in rad (array of size n)
 * @return     status code
 */
int AT_characteristic_multiple_scattering_angle( const long  	n,
												 const double 	E_MeV_u[],
												 const int		particle_charge_e[],
												 const double	target_thickness_cm[],
												 char*			element_acronym[],
												 double        	Theta_M[]);


/**
 *  Returns values f0(red_Theta) of the first Moliere function
 *
 * @param[in]  red_Theta                reduced angle variable
 * @return     first Moliere function
 */
double AT_Moliere_function_f0( double red_Theta );


/**
 *  Returns values f1(red_Theta) of the second Moliere function
 *
 * @param[in]  red_Theta                reduced angle variable
 * @return     second Moliere function
 */
double AT_Moliere_function_f1( double red_Theta );


/**
 *  Returns values f2(red_Theta) of the third Moliere function
 *
 * @param[in]  red_Theta                reduced angle variable
 * @return     third Moliere function
 */
double AT_Moliere_function_f2( double red_Theta );


/**
 *  Returns scattering angle distribution f(Theta)
 *  The distribution is not normalized because the energy loss in the target is not considered.
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_charge_e    	charge number of particle
 * @param[in]  target_thickness_cm      thickness of the target material in cm
 * @param[in]  element_acronym		    elemental symbol of target material (array of size PARTICLE_NAME_NCHAR)
 * @param[in]  Theta				   	polar angle in rad
 * @return     f(Theta)
 */
double AT_scattering_angle_distribution_single( const double E_MeV_u,
												const int    particle_charge_e,
												const double target_thickness_cm,
												const char   element_acronym[],
												double 		 Theta);


/**
 *  Returns scattering angle distribution f(Theta)
 *  The distribution is not normalized because the energy loss in the target is not considered.
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_charge_e    	charge number of particle
 * @param[in]  target_thickness_cm      thickness of the target material in cm
 * @param[in]  element_acronym		    elemental symbol of target material (array of size PARTICLE_NAME_NCHAR)
 * @param[in]  Theta				   	vector of polar angles in rad (array of size n)
 * @param[out] distribution		        vector of scattering angle distribution values (array of size n)
 * @return     status code
 */
int AT_scattering_angle_distribution( const long  	n,
									  const double 	E_MeV_u,
									  const int		particle_charge_e,
									  const double	target_thickness_cm,
									  const char	element_acronym[],
									  const double	Theta[],
									  double        distribution[]);


/**
 *  Returns Highland angle Theta0
 *
 * @param[in]  E_MeV_u                  energy of particle per nucleon [MeV]
 * @param[in]  particle_charge_e    	charge number of particle
 * @param[in]  l_over_lR     			thickness of the target material in radiation lengths
 * @return     Theta0 in rad
 */
double AT_Highland_angle_single( const double E_MeV_u,
							     const int    particle_charge_e,
							     const double l_over_lR);


/**
 *  Returns Highland angles Theta0
 *
 * @param[in]  n                        number of particles
 * @param[in]  E_MeV_u                  vector of energies of particle per nucleon [MeV] (array of size n)
 * @param[in]  particle_charge_e    	vector of charge numbers of particle (array of size n)
 * @param[in]  l_over_lR		        vector of thicknesses of target material in radiation lengths (array of size n)
 * @param[out] Theta0		            vector of Highland angles in rad (array of size n)
 * @return     status code
 */
int AT_Highland_angle( const long  		n,
					   const double 	E_MeV_u[],
					   const int		particle_charge_e[],
					   const double		l_over_lR[],
					   double        	Theta0[]);


#endif /* AT_MULTIPLECOULOMBSCATTERING_H_ */
