#ifndef AT_DATASTOPPINGPOWER_H_
#define AT_DATASTOPPINGPOWER_H_

/**
 * @brief Stopping power
 */

/*
 *    AT_DataStoppingPower.h
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

#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "AT_Error.h"
#include "AT_NumericalRoutines.h"
#include "AT_PhysicsRoutines.h"

/**
 * @enum stoppingPowerSource_no Stopping Power source code numbers
 */
enum stoppingPowerSource_no{
  PSTAR                = 0, /**< PSTAR */
  Bethe                = 1  /**< Bethe */
};


/**
 * TODO
 */
#define STOPPING_POWER_SOURCE_N    2


/**
 * TODO
 */
#define STOPPING_POWER_SOURCE_NAME_LENGTH    255

/**
 * @struct AT_table_of_stopping_power_sources
 * TODO
 */
typedef struct {
  const long    n;                                                     /**< number of sources */
  const long    stopping_power_source_no[STOPPING_POWER_SOURCE_N];     /**< TODO */
  const char*   stopping_power_source_name[STOPPING_POWER_SOURCE_N];   /**< TODO */
} AT_stopping_power_sources_struct;


/**
 * @struct AT_stopping_power_tabulated_source
 * TODO
 */
typedef struct {
  const long       number_of_data_points;                              /**< number of data points for given material and source */
  const long       stopping_power_source_no;
  const long       material_no;
  const double     E_MeV_u_and_stopping_power_total_MeV_cm2_g[][2];
} AT_stopping_power_tabulated_source_for_given_material_struct;


/**
 * @struct AT_stopping_power_tabulated_source_group_for_all_materials_struct
 * TODO
 */
typedef struct {
  const long                                                            number_of_materials;            /**< number of data points for given source */
  const long                                                            stopping_power_source_no;
  const long                                                            material_no[MATERIAL_DATA_N];
  const AT_stopping_power_tabulated_source_for_given_material_struct *  stopping_power_source_data[];
} AT_stopping_power_tabulated_source_group_for_all_materials_struct;


/**
 * @struct AT_stopping_power_tabulated_source_group_struct
 * TODO
 */
typedef struct {
  const long                                                                 number_of_sources;         /**< number of data sources */
  const long                                                                 stopping_power_source_no[STOPPING_POWER_SOURCE_N];
  const AT_stopping_power_tabulated_source_group_for_all_materials_struct *  stopping_power_source_data_group[];
} AT_stopping_power_tabulated_source_group_struct;


typedef double (*pointer_to_stopping_power_analytical_function)(const double, const long, const long);

/**
 * @struct AT_stopping_power_analytical_sources_struct
 * TODO
 */
typedef struct {
  const long                                     n;
  const long                                     stopping_power_source_no[STOPPING_POWER_SOURCE_N];
  pointer_to_stopping_power_analytical_function  access_function[STOPPING_POWER_SOURCE_N];
} AT_stopping_power_analytical_sources_struct;


/** ----------------------------------------------- FUNCTIONS --------------------------------------------- */

/**
 * TODO
 * @param[in]  source_no
 * @param[out] source_name
 */
void AT_stopping_power_source_model_name_from_number( const long source_no, char* source_name);


/**
 * TODO
 * @param[in] source_name
 * @return    source number
 */
long AT_stopping_power_source_model_number_from_name( const char* source_name );


/**
 * Wrapper for Bethe analytical function
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return
 */
double AT_Bethe_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no);


/**
 * Interpolation over tabulated Stopping Power data
 * @param[in] stopping_power_source_no
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return    gives stopping power for given energy
 */
double AT_Stopping_Power_data_interpolation(const long stopping_power_source_no,
		const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no);


/**
 * Main access method to stopping power data
 * @param[in] stopping_power_source_no
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return stopping power
 */
double AT_Stopping_Power_MeV_cm2_g_single( const long stopping_power_source_no,
		const double E_MeV_u,
		const long particle_no,
		const long material_no);


/**
 * Main access method to stopping power data - multiple field
 * @param[in] stopping_power_source_no
 * @param[in] number_of_particles
 * @param[in] E_MeV_u (array of size number_of_particles)
 * @param[in] particle_no (array of size number_of_particles)
 * @param[in] material_no
 * @param[out] Stopping_Power_MeV_cm2_g (array of size number_of_particles)
 */
void AT_Stopping_Power_MeV_cm2_g_multi( const long stopping_power_source_no,
		const long number_of_particles,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double Stopping_Power_MeV_cm2_g[]);


/**
 * Main access method to stopping power data
 * @param[in] stopping_power_source_no
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return stopping power
 */
double AT_Stopping_Power_keV_um_single( const long stopping_power_source_no,
		const double E_MeV_u,
		const long particle_no,
		const long material_no);


/**
 * Main access method to stopping power data - multiple field
 * @param[in] stopping_power_source_no
 * @param[in] number_of_particles
 * @param[in] E_MeV_u (array of size number_of_particles)
 * @param[in] particle_no (array of size number_of_particles)
 * @param[in] material_no
 * @param[out] Stopping_Power_keV_um (array of size number_of_particles)
 */
void AT_Stopping_Power_keV_um_multi( const long stopping_power_source_no,
		const long number_of_particles,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double Stopping_Power_keV_um[]);


/**
 * Computes the stopping number to be used with the Bethe formula
 * according to ICRU49, p.6, Eq. 2.3
 * BUT WITHOUT shell or density correction!
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      E_restricted_keV 	if positive and smaller than maximally transferable energy, the restricted stopping number will be computed
 * @return     result
 */
double AT_Stopping_Power_Bethe_Number(	const double 	E_MeV_u,
										const long      particle_no,
										const long 		material_no,
										const double	E_restricted_keV);


/**
 * Computes the mass stopping power using the Bethe formula
 * according to ICRU49, p.6, Eq. 2.1
 * BUT WITHOUT shell or density, Bloch or Barkas correction!
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index
 * @see          AT_DataMaterial.h for definition
 * @param[in]      E_restricted_keV 	if positive and smaller than maximally transferable energy, the restricted stopping power will be computed
 * @return     result
 */
double AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_single(	const double 	E_MeV_u,
														const long 	    particle_no,
														const long 		material_no,
														const double	E_restricted_keV);


/**
 * Computes the mass stopping power using the Bethe formula
 * for many particles according to ICRU49, p.6, Eq. 2.1
 * BUT WITHOUT shell or density, Bloch or Barkas correction!
 * @param[in]  	   n      		number of particles
 * @param[in]  	   E_MeV_u      energies of particle per nucleon (array of size n)
 * @param[in]  	   particle_no  particle indices (array of size n)
 * @see          AT_DataParticle.h for definition
 * @param[in]      material_no  material index (single value)
 * @see          AT_DataMaterial.h for definition
 * @param[in]      E_restricted_keV 	if positive and smaller than maximally transferable energy, the restricted stopping power will be computed (single value)
 * @param[out]     Mass_Stopping_Power_MeV_cm2_g (array of size n)
 */
void AT_Stopping_Power_Mass_Bethe_MeV_cm2_g_multi(	const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const double E_restricted_keV,
		double Mass_Stopping_Power_MeV_cm2_g[]);


/** ----------------------------------------------- PSTAR DATA --------------------------------------------- */

/** PSTAR data downloaded from NIST database: http://www.nist.gov/pml/data/star/index.cfm
 * Stopping-Power and Range Tables
 * for Electrons, Protons, and Helium Ions
 *
 * M.J. Berger,1) J.S. Coursey,2) M.A. Zucker2) and J. Chang2)
 *
 * 1) NIST, Physics Laboratory, Ionizing Radiation Division
 * 2) NIST, Physics Laboratory, ECSED
*/

static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_Water = {
  132,
  PSTAR,
  Water_Liquid,
  {
		  { 1.000E-03 , 1.769E+02 },
		  { 1.500E-03 , 1.984E+02 },
		  { 2.000E-03 , 2.184E+02 },
		  { 2.500E-03 , 2.370E+02 },
		  { 3.000E-03 , 2.544E+02 },
		  { 4.000E-03 , 2.864E+02 },
		  { 5.000E-03 , 3.153E+02 },
		  { 6.000E-03 , 3.420E+02 },
		  { 7.000E-03 , 3.667E+02 },
		  { 8.000E-03 , 3.900E+02 },
		  { 9.000E-03 , 4.120E+02 },
		  { 1.000E-02 , 4.329E+02 },
		  { 1.250E-02 , 4.745E+02 },
		  { 1.500E-02 , 5.110E+02 },
		  { 1.750E-02 , 5.437E+02 },
		  { 2.000E-02 , 5.733E+02 },
		  { 2.250E-02 , 6.001E+02 },
		  { 2.500E-02 , 6.245E+02 },
		  { 2.750E-02 , 6.467E+02 },
		  { 3.000E-02 , 6.671E+02 },
		  { 3.500E-02 , 7.028E+02 },
		  { 4.000E-02 , 7.324E+02 },
		  { 4.500E-02 , 7.569E+02 },
		  { 5.000E-02 , 7.768E+02 },
		  { 5.500E-02 , 7.927E+02 },
		  { 6.000E-02 , 8.050E+02 },
		  { 6.500E-02 , 8.142E+02 },
		  { 7.000E-02 , 8.205E+02 },
		  { 7.500E-02 , 8.243E+02 },
		  { 8.000E-02 , 8.260E+02 },
		  { 8.500E-02 , 8.258E+02 },
		  { 9.000E-02 , 8.239E+02 },
		  { 9.500E-02 , 8.206E+02 },
		  { 1.000E-01 , 8.161E+02 },
		  { 1.250E-01 , 7.814E+02 },
		  { 1.500E-01 , 7.371E+02 },
		  { 1.750E-01 , 6.969E+02 },
		  { 2.000E-01 , 6.613E+02 },
		  { 2.250E-01 , 6.294E+02 },
		  { 2.500E-01 , 6.006E+02 },
		  { 2.750E-01 , 5.744E+02 },
		  { 3.000E-01 , 5.504E+02 },
		  { 3.500E-01 , 5.080E+02 },
		  { 4.000E-01 , 4.719E+02 },
		  { 4.500E-01 , 4.406E+02 },
		  { 5.000E-01 , 4.132E+02 },
		  { 5.500E-01 , 3.891E+02 },
		  { 6.000E-01 , 3.680E+02 },
		  { 6.500E-01 , 3.492E+02 },
		  { 7.000E-01 , 3.325E+02 },
		  { 7.500E-01 , 3.175E+02 },
		  { 8.000E-01 , 3.039E+02 },
		  { 8.500E-01 , 2.917E+02 },
		  { 9.000E-01 , 2.805E+02 },
		  { 9.500E-01 , 2.702E+02 },
		  { 1.000E+00 , 2.608E+02 },
		  { 1.250E+00 , 2.229E+02 },
		  { 1.500E+00 , 1.957E+02 },
		  { 1.750E+00 , 1.749E+02 },
		  { 2.000E+00 , 1.586E+02 },
		  { 2.250E+00 , 1.454E+02 },
		  { 2.500E+00 , 1.344E+02 },
		  { 2.750E+00 , 1.251E+02 },
		  { 3.000E+00 , 1.172E+02 },
		  { 3.500E+00 , 1.042E+02 },
		  { 4.000E+00 , 9.404E+01 },
		  { 4.500E+00 , 8.586E+01 },
		  { 5.000E+00 , 7.911E+01 },
		  { 5.500E+00 , 7.343E+01 },
		  { 6.000E+00 , 6.858E+01 },
		  { 6.500E+00 , 6.438E+01 },
		  { 7.000E+00 , 6.071E+01 },
		  { 7.500E+00 , 5.747E+01 },
		  { 8.000E+00 , 5.460E+01 },
		  { 8.500E+00 , 5.202E+01 },
		  { 9.000E+00 , 4.969E+01 },
		  { 9.500E+00 , 4.759E+01 },
		  { 1.000E+01 , 4.567E+01 },
		  { 1.250E+01 , 3.815E+01 },
		  { 1.500E+01 , 3.292E+01 },
		  { 1.750E+01 , 2.905E+01 },
		  { 2.000E+01 , 2.607E+01 },
		  { 2.500E+01 , 2.175E+01 },
		  { 2.750E+01 , 2.013E+01 },
		  { 3.000E+01 , 1.876E+01 },
		  { 3.500E+01 , 1.656E+01 },
		  { 4.000E+01 , 1.488E+01 },
		  { 4.500E+01 , 1.354E+01 },
		  { 5.000E+01 , 1.245E+01 },
		  { 5.500E+01 , 1.154E+01 },
		  { 6.000E+01 , 1.078E+01 },
		  { 6.500E+01 , 1.013E+01 },
		  { 7.000E+01 , 9.559E+00 },
		  { 7.500E+01 , 9.063E+00 },
		  { 8.000E+01 , 8.625E+00 },
		  { 8.500E+01 , 8.236E+00 },
		  { 9.000E+01 , 7.888E+00 },
		  { 9.500E+01 , 7.573E+00 },
		  { 1.000E+02 , 7.289E+00 },
		  { 1.250E+02 , 6.192E+00 },
		  { 1.500E+02 , 5.445E+00 },
		  { 1.750E+02 , 4.903E+00 },
		  { 2.000E+02 , 4.492E+00 },
		  { 2.250E+02 , 4.170E+00 },
		  { 2.500E+02 , 3.911E+00 },
		  { 2.750E+02 , 3.698E+00 },
		  { 3.000E+02 , 3.520E+00 },
		  { 3.500E+02 , 3.241E+00 },
		  { 4.000E+02 , 3.032E+00 },
		  { 4.500E+02 , 2.871E+00 },
		  { 5.000E+02 , 2.743E+00 },
		  { 5.500E+02 , 2.640E+00 },
		  { 6.000E+02 , 2.556E+00 },
		  { 6.500E+02 , 2.485E+00 },
		  { 7.000E+02 , 2.426E+00 },
		  { 7.500E+02 , 2.376E+00 },
		  { 8.000E+02 , 2.333E+00 },
		  { 8.500E+02 , 2.296E+00 },
		  { 9.000E+02 , 2.264E+00 },
		  { 9.500E+02 , 2.236E+00 },
		  { 1.000E+03 , 2.211E+00 },
		  { 1.500E+03 , 2.070E+00 },
		  { 2.000E+03 , 2.021E+00 },
		  { 2.500E+03 , 2.004E+00 },
		  { 3.000E+03 , 2.001E+00 },
		  { 4.000E+03 , 2.012E+00 },
		  { 5.000E+03 , 2.031E+00 },
		  { 6.000E+03 , 2.052E+00 },
		  { 7.000E+03 , 2.072E+00 },
		  { 8.000E+03 , 2.091E+00 },
		  { 9.000E+03 , 2.109E+00 },
		  { 1.000E+04 , 2.126E+00 }
  }
};


static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_Aluminum_Oxide = {
  132,
  PSTAR,
  Aluminum_Oxide,
  {
		  { 1.000E-03 , 8.929E+01 },
		  { 1.500E-03 , 1.038E+02 },
		  { 2.000E-03 , 1.163E+02 },
		  { 2.500E-03 , 1.274E+02 },
		  { 3.000E-03 , 1.376E+02 },
		  { 4.000E-03 , 1.559E+02 },
		  { 5.000E-03 , 1.723E+02 },
		  { 6.000E-03 , 1.872E+02 },
		  { 7.000E-03 , 2.010E+02 },
		  { 8.000E-03 , 2.139E+02 },
		  { 9.000E-03 , 2.261E+02 },
		  { 1.000E-02 , 2.377E+02 },
		  { 1.250E-02 , 2.607E+02 },
		  { 1.500E-02 , 2.809E+02 },
		  { 1.750E-02 , 2.989E+02 },
		  { 2.000E-02 , 3.152E+02 },
		  { 2.250E-02 , 3.300E+02 },
		  { 2.500E-02 , 3.436E+02 },
		  { 2.750E-02 , 3.560E+02 },
		  { 3.000E-02 , 3.674E+02 },
		  { 3.500E-02 , 3.875E+02 },
		  { 4.000E-02 , 4.045E+02 },
		  { 4.500E-02 , 4.188E+02 },
		  { 5.000E-02 , 4.308E+02 },
		  { 5.500E-02 , 4.408E+02 },
		  { 6.000E-02 , 4.490E+02 },
		  { 6.500E-02 , 4.557E+02 },
		  { 7.000E-02 , 4.609E+02 },
		  { 7.500E-02 , 4.649E+02 },
		  { 8.000E-02 , 4.678E+02 },
		  { 8.500E-02 , 4.697E+02 },
		  { 9.000E-02 , 4.708E+02 },
		  { 9.500E-02 , 4.712E+02 },
		  { 1.000E-01 , 4.709E+02 },
		  { 1.250E-01 , 4.627E+02 },
		  { 1.500E-01 , 4.479E+02 },
		  { 1.750E-01 , 4.307E+02 },
		  { 2.000E-01 , 4.133E+02 },
		  { 2.250E-01 , 3.970E+02 },
		  { 2.500E-01 , 3.819E+02 },
		  { 2.750E-01 , 3.679E+02 },
		  { 3.000E-01 , 3.549E+02 },
		  { 3.500E-01 , 3.314E+02 },
		  { 4.000E-01 , 3.109E+02 },
		  { 4.500E-01 , 2.928E+02 },
		  { 5.000E-01 , 2.768E+02 },
		  { 5.500E-01 , 2.625E+02 },
		  { 6.000E-01 , 2.497E+02 },
		  { 6.500E-01 , 2.383E+02 },
		  { 7.000E-01 , 2.281E+02 },
		  { 7.500E-01 , 2.188E+02 },
		  { 8.000E-01 , 2.104E+02 },
		  { 8.500E-01 , 2.028E+02 },
		  { 9.000E-01 , 1.958E+02 },
		  { 9.500E-01 , 1.894E+02 },
		  { 1.000E+00 , 1.834E+02 },
		  { 1.250E+00 , 1.590E+02 },
		  { 1.500E+00 , 1.410E+02 },
		  { 1.750E+00 , 1.271E+02 },
		  { 2.000E+00 , 1.160E+02 },
		  { 2.250E+00 , 1.069E+02 },
		  { 2.500E+00 , 9.933E+01 },
		  { 2.750E+00 , 9.285E+01 },
		  { 3.000E+00 , 8.726E+01 },
		  { 3.500E+00 , 7.807E+01 },
		  { 4.000E+00 , 7.082E+01 },
		  { 4.500E+00 , 6.494E+01 },
		  { 5.000E+00 , 6.005E+01 },
		  { 5.500E+00 , 5.593E+01 },
		  { 6.000E+00 , 5.238E+01 },
		  { 6.500E+00 , 4.931E+01 },
		  { 7.000E+00 , 4.661E+01 },
		  { 7.500E+00 , 4.422E+01 },
		  { 8.000E+00 , 4.208E+01 },
		  { 8.500E+00 , 4.017E+01 },
		  { 9.000E+00 , 3.844E+01 },
		  { 9.500E+00 , 3.686E+01 },
		  { 1.000E+01 , 3.543E+01 },
		  { 1.250E+01 , 2.977E+01 },
		  { 1.500E+01 , 2.580E+01 },
		  { 1.750E+01 , 2.285E+01 },
		  { 2.000E+01 , 2.056E+01 },
		  { 2.500E+01 , 1.723E+01 },
		  { 2.750E+01 , 1.598E+01 },
		  { 3.000E+01 , 1.491E+01 },
		  { 3.500E+01 , 1.320E+01 },
		  { 4.000E+01 , 1.188E+01 },
		  { 4.500E+01 , 1.083E+01 },
		  { 5.000E+01 , 9.977E+00 },
		  { 5.500E+01 , 9.265E+00 },
		  { 6.000E+01 , 8.662E+00 },
		  { 6.500E+01 , 8.145E+00 },
		  { 7.000E+01 , 7.696E+00 },
		  { 7.500E+01 , 7.303E+00 },
		  { 8.000E+01 , 6.956E+00 },
		  { 8.500E+01 , 6.647E+00 },
		  { 9.000E+01 , 6.370E+00 },
		  { 9.500E+01 , 6.120E+00 },
		  { 1.000E+02 , 5.893E+00 },
		  { 1.250E+02 , 5.018E+00 },
		  { 1.500E+02 , 4.421E+00 },
		  { 1.750E+02 , 3.987E+00 },
		  { 2.000E+02 , 3.657E+00 },
		  { 2.250E+02 , 3.398E+00 },
		  { 2.500E+02 , 3.190E+00 },
		  { 2.750E+02 , 3.019E+00 },
		  { 3.000E+02 , 2.876E+00 },
		  { 3.500E+02 , 2.651E+00 },
		  { 4.000E+02 , 2.483E+00 },
		  { 4.500E+02 , 2.353E+00 },
		  { 5.000E+02 , 2.248E+00 },
		  { 5.500E+02 , 2.164E+00 },
		  { 6.000E+02 , 2.094E+00 },
		  { 6.500E+02 , 2.036E+00 },
		  { 7.000E+02 , 1.988E+00 },
		  { 7.500E+02 , 1.946E+00 },
		  { 8.000E+02 , 1.910E+00 },
		  { 8.500E+02 , 1.879E+00 },
		  { 9.000E+02 , 1.853E+00 },
		  { 9.500E+02 , 1.829E+00 },
		  { 1.000E+03 , 1.809E+00 },
		  { 1.500E+03 , 1.697E+00 },
		  { 2.000E+03 , 1.660E+00 },
		  { 2.500E+03 , 1.651E+00 },
		  { 3.000E+03 , 1.653E+00 },
		  { 4.000E+03 , 1.669E+00 },
		  { 5.000E+03 , 1.690E+00 },
		  { 6.000E+03 , 1.711E+00 },
		  { 7.000E+03 , 1.731E+00 },
		  { 8.000E+03 , 1.749E+00 },
		  { 9.000E+03 , 1.766E+00 },
		  { 1.000E+04 , 1.782E+00 }
  }
};


static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_Aluminum = {
  132,
  PSTAR,
  Aluminum,
  {
		  { 1.000E-03 , 1.043E+02 },
		  { 1.500E-03 , 1.239E+02 },
		  { 2.000E-03 , 1.404E+02 },
		  { 2.500E-03 , 1.550E+02 },
		  { 3.000E-03 , 1.683E+02 },
		  { 4.000E-03 , 1.921E+02 },
		  { 5.000E-03 , 2.131E+02 },
		  { 6.000E-03 , 2.323E+02 },
		  { 7.000E-03 , 2.499E+02 },
		  { 8.000E-03 , 2.664E+02 },
		  { 9.000E-03 , 2.819E+02 },
		  { 1.000E-02 , 2.966E+02 },
		  { 1.250E-02 , 3.245E+02 },
		  { 1.500E-02 , 3.483E+02 },
		  { 1.750E-02 , 3.689E+02 },
		  { 2.000E-02 , 3.867E+02 },
		  { 2.250E-02 , 4.022E+02 },
		  { 2.500E-02 , 4.157E+02 },
		  { 2.750E-02 , 4.273E+02 },
		  { 3.000E-02 , 4.373E+02 },
		  { 3.500E-02 , 4.529E+02 },
		  { 4.000E-02 , 4.638E+02 },
		  { 4.500E-02 , 4.709E+02 },
		  { 5.000E-02 , 4.749E+02 },
		  { 5.500E-02 , 4.766E+02 },
		  { 6.000E-02 , 4.764E+02 },
		  { 6.500E-02 , 4.749E+02 },
		  { 7.000E-02 , 4.724E+02 },
		  { 7.500E-02 , 4.691E+02 },
		  { 8.000E-02 , 4.653E+02 },
		  { 8.500E-02 , 4.611E+02 },
		  { 9.000E-02 , 4.567E+02 },
		  { 9.500E-02 , 4.522E+02 },
		  { 1.000E-01 , 4.477E+02 },
		  { 1.250E-01 , 4.253E+02 },
		  { 1.500E-01 , 4.051E+02 },
		  { 1.750E-01 , 3.873E+02 },
		  { 2.000E-01 , 3.715E+02 },
		  { 2.250E-01 , 3.573E+02 },
		  { 2.500E-01 , 3.444E+02 },
		  { 2.750E-01 , 3.327E+02 },
		  { 3.000E-01 , 3.218E+02 },
		  { 3.500E-01 , 3.020E+02 },
		  { 4.000E-01 , 2.844E+02 },
		  { 4.500E-01 , 2.689E+02 },
		  { 5.000E-01 , 2.550E+02 },
		  { 5.500E-01 , 2.427E+02 },
		  { 6.000E-01 , 2.316E+02 },
		  { 6.500E-01 , 2.216E+02 },
		  { 7.000E-01 , 2.126E+02 },
		  { 7.500E-01 , 2.043E+02 },
		  { 8.000E-01 , 1.968E+02 },
		  { 8.500E-01 , 1.899E+02 },
		  { 9.000E-01 , 1.835E+02 },
		  { 9.500E-01 , 1.775E+02 },
		  { 1.000E+00 , 1.720E+02 },
		  { 1.250E+00 , 1.495E+02 },
		  { 1.500E+00 , 1.328E+02 },
		  { 1.750E+00 , 1.199E+02 },
		  { 2.000E+00 , 1.095E+02 },
		  { 2.250E+00 , 1.010E+02 },
		  { 2.500E+00 , 9.383E+01 },
		  { 2.750E+00 , 8.775E+01 },
		  { 3.000E+00 , 8.250E+01 },
		  { 3.500E+00 , 7.388E+01 },
		  { 4.000E+00 , 6.707E+01 },
		  { 4.500E+00 , 6.154E+01 },
		  { 5.000E+00 , 5.695E+01 },
		  { 5.500E+00 , 5.306E+01 },
		  { 6.000E+00 , 4.973E+01 },
		  { 6.500E+00 , 4.684E+01 },
		  { 7.000E+00 , 4.430E+01 },
		  { 7.500E+00 , 4.205E+01 },
		  { 8.000E+00 , 4.004E+01 },
		  { 8.500E+00 , 3.824E+01 },
		  { 9.000E+00 , 3.660E+01 },
		  { 9.500E+00 , 3.512E+01 },
		  { 1.000E+01 , 3.376E+01 },
		  { 1.250E+01 , 2.842E+01 },
		  { 1.500E+01 , 2.466E+01 },
		  { 1.750E+01 , 2.186E+01 },
		  { 2.000E+01 , 1.969E+01 },
		  { 2.500E+01 , 1.652E+01 },
		  { 2.750E+01 , 1.532E+01 },
		  { 3.000E+01 , 1.431E+01 },
		  { 3.500E+01 , 1.268E+01 },
		  { 4.000E+01 , 1.142E+01 },
		  { 4.500E+01 , 1.041E+01 },
		  { 5.000E+01 , 9.594E+00 },
		  { 5.500E+01 , 8.911E+00 },
		  { 6.000E+01 , 8.334E+00 },
		  { 6.500E+01 , 7.838E+00 },
		  { 7.000E+01 , 7.408E+00 },
		  { 7.500E+01 , 7.031E+00 },
		  { 8.000E+01 , 6.698E+00 },
		  { 8.500E+01 , 6.401E+00 },
		  { 9.000E+01 , 6.135E+00 },
		  { 9.500E+01 , 5.895E+00 },
		  { 1.000E+02 , 5.678E+00 },
		  { 1.250E+02 , 4.837E+00 },
		  { 1.500E+02 , 4.262E+00 },
		  { 1.750E+02 , 3.844E+00 },
		  { 2.000E+02 , 3.526E+00 },
		  { 2.250E+02 , 3.277E+00 },
		  { 2.500E+02 , 3.076E+00 },
		  { 2.750E+02 , 2.911E+00 },
		  { 3.000E+02 , 2.773E+00 },
		  { 3.500E+02 , 2.555E+00 },
		  { 4.000E+02 , 2.393E+00 },
		  { 4.500E+02 , 2.267E+00 },
		  { 5.000E+02 , 2.167E+00 },
		  { 5.500E+02 , 2.086E+00 },
		  { 6.000E+02 , 2.020E+00 },
		  { 6.500E+02 , 1.965E+00 },
		  { 7.000E+02 , 1.918E+00 },
		  { 7.500E+02 , 1.879E+00 },
		  { 8.000E+02 , 1.845E+00 },
		  { 8.500E+02 , 1.816E+00 },
		  { 9.000E+02 , 1.791E+00 },
		  { 9.500E+02 , 1.769E+00 },
		  { 1.000E+03 , 1.750E+00 },
		  { 1.500E+03 , 1.647E+00 },
		  { 2.000E+03 , 1.618E+00 },
		  { 2.500E+03 , 1.613E+00 },
		  { 3.000E+03 , 1.619E+00 },
		  { 4.000E+03 , 1.642E+00 },
		  { 5.000E+03 , 1.668E+00 },
		  { 6.000E+03 , 1.692E+00 },
		  { 7.000E+03 , 1.714E+00 },
		  { 8.000E+03 , 1.734E+00 },
		  { 9.000E+03 , 1.752E+00 },
		  { 1.000E+04 , 1.768E+00 }
  }
};


static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_PMMA = {
  0,
  PSTAR,
  PMMA,
  {
  }
};


static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_Alanine = {
  0,
  PSTAR,
  Alanine,
  {
  }
};


static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_LiF = {
  132,
  PSTAR,
  LiF,
  {
		  { 1.000E-03 , 1.008E+02 },
		  { 1.500E-03 , 1.157E+02 },
		  { 2.000E-03 , 1.289E+02 },
		  { 2.500E-03 , 1.408E+02 },
		  { 3.000E-03 , 1.518E+02 },
		  { 4.000E-03 , 1.717E+02 },
		  { 5.000E-03 , 1.895E+02 },
		  { 6.000E-03 , 2.058E+02 },
		  { 7.000E-03 , 2.210E+02 },
		  { 8.000E-03 , 2.352E+02 },
		  { 9.000E-03 , 2.485E+02 },
		  { 1.000E-02 , 2.613E+02 },
		  { 1.250E-02 , 2.866E+02 },
		  { 1.500E-02 , 3.090E+02 },
		  { 1.750E-02 , 3.291E+02 },
		  { 2.000E-02 , 3.474E+02 },
		  { 2.250E-02 , 3.640E+02 },
		  { 2.500E-02 , 3.793E+02 },
		  { 2.750E-02 , 3.934E+02 },
		  { 3.000E-02 , 4.065E+02 },
		  { 3.500E-02 , 4.298E+02 },
		  { 4.000E-02 , 4.498E+02 },
		  { 4.500E-02 , 4.671E+02 },
		  { 5.000E-02 , 4.820E+02 },
		  { 5.500E-02 , 4.947E+02 },
		  { 6.000E-02 , 5.055E+02 },
		  { 6.500E-02 , 5.147E+02 },
		  { 7.000E-02 , 5.223E+02 },
		  { 7.500E-02 , 5.286E+02 },
		  { 8.000E-02 , 5.337E+02 },
		  { 8.500E-02 , 5.376E+02 },
		  { 9.000E-02 , 5.406E+02 },
		  { 9.500E-02 , 5.427E+02 },
		  { 1.000E-01 , 5.440E+02 },
		  { 1.250E-01 , 5.411E+02 },
		  { 1.500E-01 , 5.281E+02 },
		  { 1.750E-01 , 5.098E+02 },
		  { 2.000E-01 , 4.893E+02 },
		  { 2.250E-01 , 4.679E+02 },
		  { 2.500E-01 , 4.469E+02 },
		  { 2.750E-01 , 4.273E+02 },
		  { 3.000E-01 , 4.091E+02 },
		  { 3.500E-01 , 3.773E+02 },
		  { 4.000E-01 , 3.507E+02 },
		  { 4.500E-01 , 3.283E+02 },
		  { 5.000E-01 , 3.093E+02 },
		  { 5.500E-01 , 2.930E+02 },
		  { 6.000E-01 , 2.783E+02 },
		  { 6.500E-01 , 2.652E+02 },
		  { 7.000E-01 , 2.534E+02 },
		  { 7.500E-01 , 2.427E+02 },
		  { 8.000E-01 , 2.329E+02 },
		  { 8.500E-01 , 2.240E+02 },
		  { 9.000E-01 , 2.158E+02 },
		  { 9.500E-01 , 2.083E+02 },
		  { 1.000E+00 , 2.014E+02 },
		  { 1.250E+00 , 1.732E+02 },
		  { 1.500E+00 , 1.526E+02 },
		  { 1.750E+00 , 1.369E+02 },
		  { 2.000E+00 , 1.244E+02 },
		  { 2.250E+00 , 1.142E+02 },
		  { 2.500E+00 , 1.058E+02 },
		  { 2.750E+00 , 9.859E+01 },
		  { 3.000E+00 , 9.243E+01 },
		  { 3.500E+00 , 8.237E+01 },
		  { 4.000E+00 , 7.447E+01 },
		  { 4.500E+00 , 6.807E+01 },
		  { 5.000E+00 , 6.279E+01 },
		  { 5.500E+00 , 5.834E+01 },
		  { 6.000E+00 , 5.453E+01 },
		  { 6.500E+00 , 5.124E+01 },
		  { 7.000E+00 , 4.836E+01 },
		  { 7.500E+00 , 4.581E+01 },
		  { 8.000E+00 , 4.354E+01 },
		  { 8.500E+00 , 4.151E+01 },
		  { 9.000E+00 , 3.967E+01 },
		  { 9.500E+00 , 3.801E+01 },
		  { 1.000E+01 , 3.649E+01 },
		  { 1.250E+01 , 3.054E+01 },
		  { 1.500E+01 , 2.639E+01 },
		  { 1.750E+01 , 2.331E+01 },
		  { 2.000E+01 , 2.093E+01 },
		  { 2.500E+01 , 1.749E+01 },
		  { 2.750E+01 , 1.620E+01 },
		  { 3.000E+01 , 1.510E+01 },
		  { 3.500E+01 , 1.334E+01 },
		  { 4.000E+01 , 1.199E+01 },
		  { 4.500E+01 , 1.092E+01 },
		  { 5.000E+01 , 1.005E+01 },
		  { 5.500E+01 , 9.320E+00 },
		  { 6.000E+01 , 8.707E+00 },
		  { 6.500E+01 , 8.181E+00 },
		  { 7.000E+01 , 7.725E+00 },
		  { 7.500E+01 , 7.327E+00 },
		  { 8.000E+01 , 6.975E+00 },
		  { 8.500E+01 , 6.661E+00 },
		  { 9.000E+01 , 6.381E+00 },
		  { 9.500E+01 , 6.128E+00 },
		  { 1.000E+02 , 5.899E+00 },
		  { 1.250E+02 , 5.015E+00 },
		  { 1.500E+02 , 4.412E+00 },
		  { 1.750E+02 , 3.975E+00 },
		  { 2.000E+02 , 3.643E+00 },
		  { 2.250E+02 , 3.383E+00 },
		  { 2.500E+02 , 3.174E+00 },
		  { 2.750E+02 , 3.002E+00 },
		  { 3.000E+02 , 2.858E+00 },
		  { 3.500E+02 , 2.633E+00 },
		  { 4.000E+02 , 2.464E+00 },
		  { 4.500E+02 , 2.333E+00 },
		  { 5.000E+02 , 2.229E+00 },
		  { 5.500E+02 , 2.144E+00 },
		  { 6.000E+02 , 2.075E+00 },
		  { 6.500E+02 , 2.017E+00 },
		  { 7.000E+02 , 1.968E+00 },
		  { 7.500E+02 , 1.927E+00 },
		  { 8.000E+02 , 1.891E+00 },
		  { 8.500E+02 , 1.860E+00 },
		  { 9.000E+02 , 1.833E+00 },
		  { 9.500E+02 , 1.810E+00 },
		  { 1.000E+03 , 1.789E+00 },
		  { 1.500E+03 , 1.675E+00 },
		  { 2.000E+03 , 1.635E+00 },
		  { 2.500E+03 , 1.622E+00 },
		  { 3.000E+03 , 1.621E+00 },
		  { 4.000E+03 , 1.631E+00 },
		  { 5.000E+03 , 1.648E+00 },
		  { 6.000E+03 , 1.665E+00 },
		  { 7.000E+03 , 1.682E+00 },
		  { 8.000E+03 , 1.697E+00 },
		  { 9.000E+03 , 1.712E+00 },
		  { 1.000E+04 , 1.726E+00 }
  }
};


static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_Air = {
  132,
  PSTAR,
  Air,
  {
		  { 1.000E-03 , 1.414E+02 },
		  { 1.500E-03 , 1.651E+02 },
		  { 2.000E-03 , 1.855E+02 },
		  { 2.500E-03 , 2.038E+02 },
		  { 3.000E-03 , 2.206E+02 },
		  { 4.000E-03 , 2.507E+02 },
		  { 5.000E-03 , 2.776E+02 },
		  { 6.000E-03 , 3.021E+02 },
		  { 7.000E-03 , 3.248E+02 },
		  { 8.000E-03 , 3.460E+02 },
		  { 9.000E-03 , 3.660E+02 },
		  { 1.000E-02 , 3.850E+02 },
		  { 1.250E-02 , 4.224E+02 },
		  { 1.500E-02 , 4.552E+02 },
		  { 1.750E-02 , 4.843E+02 },
		  { 2.000E-02 , 5.106E+02 },
		  { 2.250E-02 , 5.343E+02 },
		  { 2.500E-02 , 5.558E+02 },
		  { 2.750E-02 , 5.755E+02 },
		  { 3.000E-02 , 5.934E+02 },
		  { 3.500E-02 , 6.246E+02 },
		  { 4.000E-02 , 6.506E+02 },
		  { 4.500E-02 , 6.721E+02 },
		  { 5.000E-02 , 6.897E+02 },
		  { 5.500E-02 , 7.038E+02 },
		  { 6.000E-02 , 7.149E+02 },
		  { 6.500E-02 , 7.233E+02 },
		  { 7.000E-02 , 7.293E+02 },
		  { 7.500E-02 , 7.333E+02 },
		  { 8.000E-02 , 7.355E+02 },
		  { 8.500E-02 , 7.360E+02 },
		  { 9.000E-02 , 7.352E+02 },
		  { 9.500E-02 , 7.332E+02 },
		  { 1.000E-01 , 7.301E+02 },
		  { 1.250E-01 , 7.038E+02 },
		  { 1.500E-01 , 6.680E+02 },
		  { 1.750E-01 , 6.298E+02 },
		  { 2.000E-01 , 5.928E+02 },
		  { 2.250E-01 , 5.589E+02 },
		  { 2.500E-01 , 5.284E+02 },
		  { 2.750E-01 , 5.011E+02 },
		  { 3.000E-01 , 4.767E+02 },
		  { 3.500E-01 , 4.353E+02 },
		  { 4.000E-01 , 4.015E+02 },
		  { 4.500E-01 , 3.736E+02 },
		  { 5.000E-01 , 3.501E+02 },
		  { 5.500E-01 , 3.300E+02 },
		  { 6.000E-01 , 3.123E+02 },
		  { 6.500E-01 , 2.967E+02 },
		  { 7.000E-01 , 2.826E+02 },
		  { 7.500E-01 , 2.701E+02 },
		  { 8.000E-01 , 2.589E+02 },
		  { 8.500E-01 , 2.486E+02 },
		  { 9.000E-01 , 2.393E+02 },
		  { 9.500E-01 , 2.308E+02 },
		  { 1.000E+00 , 2.229E+02 },
		  { 1.250E+00 , 1.912E+02 },
		  { 1.500E+00 , 1.683E+02 },
		  { 1.750E+00 , 1.509E+02 },
		  { 2.000E+00 , 1.371E+02 },
		  { 2.250E+00 , 1.258E+02 },
		  { 2.500E+00 , 1.165E+02 },
		  { 2.750E+00 , 1.086E+02 },
		  { 3.000E+00 , 1.018E+02 },
		  { 3.500E+00 , 9.068E+01 },
		  { 4.000E+00 , 8.197E+01 },
		  { 4.500E+00 , 7.492E+01 },
		  { 5.000E+00 , 6.909E+01 },
		  { 5.500E+00 , 6.417E+01 },
		  { 6.000E+00 , 5.997E+01 },
		  { 6.500E+00 , 5.633E+01 },
		  { 7.000E+00 , 5.315E+01 },
		  { 7.500E+00 , 5.033E+01 },
		  { 8.000E+00 , 4.783E+01 },
		  { 8.500E+00 , 4.559E+01 },
		  { 9.000E+00 , 4.357E+01 },
		  { 9.500E+00 , 4.173E+01 },
		  { 1.000E+01 , 4.006E+01 },
		  { 1.250E+01 , 3.351E+01 },
		  { 1.500E+01 , 2.894E+01 },
		  { 1.750E+01 , 2.555E+01 },
		  { 2.000E+01 , 2.294E+01 },
		  { 2.500E+01 , 1.915E+01 },
		  { 2.750E+01 , 1.773E+01 },
		  { 3.000E+01 , 1.653E+01 },
		  { 3.500E+01 , 1.460E+01 },
		  { 4.000E+01 , 1.312E+01 },
		  { 4.500E+01 , 1.194E+01 },
		  { 5.000E+01 , 1.099E+01 },
		  { 5.500E+01 , 1.019E+01 },
		  { 6.000E+01 , 9.517E+00 },
		  { 6.500E+01 , 8.942E+00 },
		  { 7.000E+01 , 8.443E+00 },
		  { 7.500E+01 , 8.006E+00 },
		  { 8.000E+01 , 7.620E+00 },
		  { 8.500E+01 , 7.277E+00 },
		  { 9.000E+01 , 6.970E+00 },
		  { 9.500E+01 , 6.693E+00 },
		  { 1.000E+02 , 6.443E+00 },
		  { 1.250E+02 , 5.475E+00 },
		  { 1.500E+02 , 4.816E+00 },
		  { 1.750E+02 , 4.338E+00 },
		  { 2.000E+02 , 3.976E+00 },
		  { 2.250E+02 , 3.691E+00 },
		  { 2.500E+02 , 3.462E+00 },
		  { 2.750E+02 , 3.275E+00 },
		  { 3.000E+02 , 3.118E+00 },
		  { 3.500E+02 , 2.871E+00 },
		  { 4.000E+02 , 2.687E+00 },
		  { 4.500E+02 , 2.544E+00 },
		  { 5.000E+02 , 2.431E+00 },
		  { 5.500E+02 , 2.340E+00 },
		  { 6.000E+02 , 2.266E+00 },
		  { 6.500E+02 , 2.203E+00 },
		  { 7.000E+02 , 2.151E+00 },
		  { 7.500E+02 , 2.107E+00 },
		  { 8.000E+02 , 2.069E+00 },
		  { 8.500E+02 , 2.037E+00 },
		  { 9.000E+02 , 2.008E+00 },
		  { 9.500E+02 , 1.984E+00 },
		  { 1.000E+03 , 1.963E+00 },
		  { 1.500E+03 , 1.850E+00 },
		  { 2.000E+03 , 1.820E+00 },
		  { 2.500E+03 , 1.818E+00 },
		  { 3.000E+03 , 1.828E+00 },
		  { 4.000E+03 , 1.861E+00 },
		  { 5.000E+03 , 1.898E+00 },
		  { 6.000E+03 , 1.934E+00 },
		  { 7.000E+03 , 1.967E+00 },
		  { 8.000E+03 , 1.998E+00 },
		  { 9.000E+03 , 2.026E+00 },
		  { 1.000E+04 , 2.052E+00 }
  }
};


/** ----------------------------------------------- Static objects  --------------------------------------- */


static const AT_stopping_power_sources_struct AT_stopping_power_sources = {
	STOPPING_POWER_SOURCE_N,
    {  PSTAR,         Bethe           }, // source_no
    {  "PSTAR data",  "Bethe formula" }  // source_name
};


static const AT_stopping_power_tabulated_source_group_for_all_materials_struct AT_data_PSTAR_source = {
  MATERIAL_DATA_N,
  PSTAR,
  { User_Defined_Material,  Water_Liquid,                         Aluminum_Oxide,                                 Aluminum,                                 PMMA,                                Alanine,                                LiF,                                Air },
  { NULL,                   &AT_stopping_power_data_PSTAR_Water,  &AT_stopping_power_data_PSTAR_Aluminum_Oxide,   &AT_stopping_power_data_PSTAR_Aluminum,   &AT_stopping_power_data_PSTAR_PMMA,  &AT_stopping_power_data_PSTAR_Alanine,  &AT_stopping_power_data_PSTAR_LiF,  &AT_stopping_power_data_PSTAR_Air}
};


static const AT_stopping_power_tabulated_source_group_struct AT_stopping_power_tabulated_source = {
  STOPPING_POWER_SOURCE_N,
  { PSTAR,                  Bethe},
  { &AT_data_PSTAR_source,  NULL}
};


static const AT_stopping_power_analytical_sources_struct AT_stopping_power_analytical_source = {
  STOPPING_POWER_SOURCE_N,
  { PSTAR,   Bethe},
  { NULL,    &AT_Bethe_wrapper}
};

#endif /* AT_DATASTOPPINGPOWER_H_ */
