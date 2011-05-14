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
  Bethe                = 1, /**< Bethe */
  ShieldHit			   = 2  /**< ShieldHit code: extended Bethe formula */
};


/**
 * TODO
 */
#define STOPPING_POWER_SOURCE_N    3


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
 * Wrapper for ShieldHit stopping powers
 * momentarily read from hardwired struct for water
 * TODO: augment by routine to read data file into
 * TODO: struct and look-up value
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return stopping power (MeV cm2 per g)
 */
double AT_ShieldHit_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no);


/**
 * TODO
 * @param[in] stopping_power_source_no
 * @param[in] particle_no
 * @param[in] material_no
 * @param[out] source_for_given_material
 * @return
 */
double _AT_Stopping_Power_get_data(const long stopping_power_source_no,
		const long 	    particle_no,
		const long 		material_no,
		AT_stopping_power_tabulated_source_for_given_material_struct ** source_for_given_material);

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
 * TODO
 * @param[in] stopping_power_source_no
 * @param[in] Stopping_Power_MeV_cm2_g
 * @param[in] particle_no
 * @param[in] material_no
 * @return range [m]
 */
double AT_Energy_MeV_u_from_Stopping_Power_single( const long stopping_power_source_no,
		const double Stopping_Power_MeV_cm2_g,
		const long particle_no,
		const long material_no);


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

/**
 * Computes the Rutherford single differential cross section
 * for the energy spectrum of secondary electrons produced by
 * an HCP
 * @param[in]  	   E_MeV_u      energy of particle per nucleon
 * @param[in]  	   particle_no  particle index
 * @param[in]      material_no  material index
 * @param[in]  	   n      		number of secondary electron energies
 * @param[in]      T_MeV 	    electron energies (array of size n)
 * @param[out]     dsdT_m2_MeV  Rutherford SDCS for given electron energies (array of size n)
 * @return         status code
 */
int AT_Rutherford_SDCS( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const long n,
		const double T_MeV[],
		double dsdT_m2_MeV[]);

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
  132,
  PSTAR,
  PMMA,
  {
		  { 0.001 , 214.7 },
		  { 0.0015 , 246.3 },
		  { 0.002 , 274.6 },
		  { 0.0025 , 300.4 },
		  { 0.003 , 324.2 },
		  { 0.004 , 363.5 },
		  { 0.005 , 399.4 },
		  { 0.006 , 432.2 },
		  { 0.007 , 462.4 },
		  { 0.008 , 489.9 },
		  { 0.009 , 515.3 },
		  { 0.01 , 539.1 },
		  { 0.0125 , 587 },
		  { 0.015 , 627.7 },
		  { 0.0175 , 663.2 },
		  { 0.02 , 694.8 },
		  { 0.0225 , 722.9 },
		  { 0.025 , 747.9 },
		  { 0.0275 , 770 },
		  { 0.03 , 789.8 },
		  { 0.035 , 824.4 },
		  { 0.04 , 853.7 },
		  { 0.045 , 877.9 },
		  { 0.05 , 897.2 },
		  { 0.055 , 911.8 },
		  { 0.06 , 922.4 },
		  { 0.065 , 929.6 },
		  { 0.07 , 934 },
		  { 0.075 , 935.8 },
		  { 0.08 , 935.5 },
		  { 0.085 , 933.4 },
		  { 0.09 , 929.6 },
		  { 0.095 , 924.5 },
		  { 0.1 , 918.3 },
		  { 0.125 , 875.7 },
		  { 0.15 , 824.3 },
		  { 0.175 , 772.3 },
		  { 0.2 , 723.2 },
		  { 0.225 , 676.9 },
		  { 0.25 , 634.6 },
		  { 0.275 , 596.9 },
		  { 0.3 , 563.4 },
		  { 0.35 , 507.7 },
		  { 0.4 , 463.9 },
		  { 0.45 , 428.9 },
		  { 0.5 , 400.6 },
		  { 0.55 , 377 },
		  { 0.6 , 356.4 },
		  { 0.65 , 338.3 },
		  { 0.7 , 322.2 },
		  { 0.75 , 307.8 },
		  { 0.8 , 294.8 },
		  { 0.85 , 283 },
		  { 0.9 , 272.2 },
		  { 0.95 , 262.3 },
		  { 1 , 253.2 },
		  { 1.25 , 216.8 },
		  { 1.5 , 190.5 },
		  { 1.75 , 170.5 },
		  { 2 , 154.6 },
		  { 2.25 , 141.8 },
		  { 2.5 , 131.1 },
		  { 2.75 , 122.1 },
		  { 3 , 114.3 },
		  { 3.5 , 101.7 },
		  { 4 , 91.79 },
		  { 4.5 , 83.79 },
		  { 5 , 77.19 },
		  { 5.5 , 71.64 },
		  { 6 , 66.9 },
		  { 6.5 , 62.8 },
		  { 7 , 59.21 },
		  { 7.5 , 56.05 },
		  { 8 , 53.24 },
		  { 8.5 , 50.73 },
		  { 9 , 48.45 },
		  { 9.5 , 46.4 },
		  { 10 , 44.52 },
		  { 12.5 , 37.19 },
		  { 15 , 32.08 },
		  { 17.5 , 28.31 },
		  { 20 , 25.39 },
		  { 25 , 21.18 },
		  { 27.5 , 19.61 },
		  { 30 , 18.27 },
		  { 35 , 16.13 },
		  { 40 , 14.49 },
		  { 45 , 13.18 },
		  { 50 , 12.12 },
		  { 55 , 11.24 },
		  { 60 , 10.5 },
		  { 65 , 9.858 },
		  { 70 , 9.306 },
		  { 75 , 8.823 },
		  { 80 , 8.397 },
		  { 85 , 8.018 },
		  { 90 , 7.678 },
		  { 95 , 7.372 },
		  { 100 , 7.095 },
		  { 125 , 6.027 },
		  { 150 , 5.3 },
		  { 175 , 4.772 },
		  { 200 , 4.372 },
		  { 225 , 4.058 },
		  { 250 , 3.806 },
		  { 275 , 3.599 },
		  { 300 , 3.426 },
		  { 350 , 3.154 },
		  { 400 , 2.951 },
		  { 450 , 2.794 },
		  { 500 , 2.67 },
		  { 550 , 2.569 },
		  { 600 , 2.487 },
		  { 650 , 2.418 },
		  { 700 , 2.361 },
		  { 750 , 2.312 },
		  { 800 , 2.269 },
		  { 850 , 2.232 },
		  { 900 , 2.199 },
		  { 950 , 2.17 },
		  { 1000 , 2.145 },
		  { 1500 , 2.004 },
		  { 2000 , 1.954 },
		  { 2500 , 1.937 },
		  { 3000 , 1.934 },
		  { 4000 , 1.945 },
		  { 5000 , 1.964 },
		  { 6000 , 1.985 },
		  { 7000 , 2.005 },
		  { 8000 , 2.024 },
		  { 9000 , 2.042 },
		  { 10000 , 2.059 }
  }
};


static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_Alanine = {
  132,
  PSTAR,
  Alanine,
  {
		  { 0.001 , 209.5 },
		  { 0.0015 , 240.4 },
		  { 0.002 , 268 },
		  { 0.0025 , 293.2 },
		  { 0.003 , 316.5 },
		  { 0.004 , 356.1 },
		  { 0.005 , 392 },
		  { 0.006 , 424.7 },
		  { 0.007 , 455.1 },
		  { 0.008 , 483 },
		  { 0.009 , 508.8 },
		  { 0.01 , 533.3 },
		  { 0.0125 , 581.7 },
		  { 0.015 , 623.1 },
		  { 0.0175 , 659.5 },
		  { 0.02 , 691.7 },
		  { 0.0225 , 720.6 },
		  { 0.025 , 746.2 },
		  { 0.0275 , 769.1 },
		  { 0.03 , 789.6 },
		  { 0.035 , 825.2 },
		  { 0.04 , 854.8 },
		  { 0.045 , 879 },
		  { 0.05 , 898.1 },
		  { 0.055 , 912.6 },
		  { 0.06 , 923.1 },
		  { 0.065 , 930.3 },
		  { 0.07 , 934.5 },
		  { 0.075 , 936.3 },
		  { 0.08 , 935.9 },
		  { 0.085 , 933.6 },
		  { 0.09 , 929.7 },
		  { 0.095 , 924.6 },
		  { 0.1 , 918.3 },
		  { 0.125 , 875.3 },
		  { 0.15 , 823.7 },
		  { 0.175 , 771.5 },
		  { 0.2 , 722.4 },
		  { 0.225 , 677.3 },
		  { 0.25 , 636.9 },
		  { 0.275 , 601.1 },
		  { 0.3 , 569.4 },
		  { 0.35 , 517.2 },
		  { 0.4 , 475.4 },
		  { 0.45 , 441.2 },
		  { 0.5 , 412.3 },
		  { 0.55 , 387.5 },
		  { 0.6 , 365.9 },
		  { 0.65 , 346.9 },
		  { 0.7 , 330.1 },
		  { 0.75 , 315 },
		  { 0.8 , 301.5 },
		  { 0.85 , 289.2 },
		  { 0.9 , 278 },
		  { 0.95 , 267.7 },
		  { 1 , 258.4 },
		  { 1.25 , 220.7 },
		  { 1.5 , 193.7 },
		  { 1.75 , 173.1 },
		  { 2 , 157 },
		  { 2.25 , 143.8 },
		  { 2.5 , 133 },
		  { 2.75 , 123.7 },
		  { 3 , 115.8 },
		  { 3.5 , 103 },
		  { 4 , 92.89 },
		  { 4.5 , 84.78 },
		  { 5 , 78.08 },
		  { 5.5 , 72.45 },
		  { 6 , 67.64 },
		  { 6.5 , 63.48 },
		  { 7 , 59.85 },
		  { 7.5 , 56.64 },
		  { 8 , 53.79 },
		  { 8.5 , 51.25 },
		  { 9 , 48.95 },
		  { 9.5 , 46.87 },
		  { 10 , 44.96 },
		  { 12.5 , 37.55 },
		  { 15 , 32.38 },
		  { 17.5 , 28.56 },
		  { 20 , 25.62 },
		  { 25 , 21.36 },
		  { 27.5 , 19.77 },
		  { 30 , 18.42 },
		  { 35 , 16.26 },
		  { 40 , 14.6 },
		  { 45 , 13.28 },
		  { 50 , 12.21 },
		  { 55 , 11.32 },
		  { 60 , 10.57 },
		  { 65 , 9.929 },
		  { 70 , 9.372 },
		  { 75 , 8.886 },
		  { 80 , 8.455 },
		  { 85 , 8.073 },
		  { 90 , 7.73 },
		  { 95 , 7.422 },
		  { 100 , 7.142 },
		  { 125 , 6.065 },
		  { 150 , 5.331 },
		  { 175 , 4.8 },
		  { 200 , 4.396 },
		  { 225 , 4.08 },
		  { 250 , 3.825 },
		  { 275 , 3.616 },
		  { 300 , 3.441 },
		  { 350 , 3.167 },
		  { 400 , 2.961 },
		  { 450 , 2.802 },
		  { 500 , 2.677 },
		  { 550 , 2.575 },
		  { 600 , 2.491 },
		  { 650 , 2.421 },
		  { 700 , 2.362 },
		  { 750 , 2.313 },
		  { 800 , 2.27 },
		  { 850 , 2.233 },
		  { 900 , 2.2 },
		  { 950 , 2.173 },
		  { 1000 , 2.148 },
		  { 1500 , 2.016 },
		  { 2000 , 1.976 },
		  { 2500 , 1.968 },
		  { 3000 , 1.973 },
		  { 4000 , 2 },
		  { 5000 , 2.032 },
		  { 6000 , 2.064 },
		  { 7000 , 2.094 },
		  { 8000 , 2.122 },
		  { 9000 , 2.149 },
		  { 10000 , 2.173 }
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

/** ----------------------------------------------- SHIELD HIT DATA --------------------------------------------- */

/** For the time being this array holds the data
 * as given by the extended Bethe formular of
 * SHIELD-HIT for water
 * TODO: Later use this array to read in data
 * TODO: from table
*/

/**
 * @struct AT_stopping_power_ShieldHit_table
 * TODO
 */
typedef struct {
  const long       material_no;
  const long       number_of_data_points;
  const double     E_MeV_u_and_stopping_power_total_MeV_cm2_g[53][19];
} AT_stopping_power_ShieldHit_table_struct;

static const AT_stopping_power_ShieldHit_table_struct AT_stopping_power_ShieldHit_table = {
  Water_Liquid,
  53,
  {
	  {0.025, 6.2450E+02, 1.151E+03, 2.3193E+03, 2.8720E+03, 3.2922E+03, 3.6037E+03, 3.8821E+03, 4.1215E+03, 4.2951E+03, 4.4513E+03, 4.5914E+03, 4.7537E+03, 4.9110E+03, 5.0697E+03, 5.2616E+03, 5.4129E+03, 5.6171E+03, 5.7158E+03},
	  {0.03,  6.6710E+02, 1.257E+03, 2.5198E+03, 3.1439E+03, 3.6315E+03, 4.0035E+03, 4.3278E+03, 4.6063E+03, 4.8107E+03, 4.9910E+03, 5.1553E+03, 5.3178E+03, 5.4856E+03, 5.6529E+03, 5.8538E+03, 6.0193E+03, 6.2307E+03, 6.3394E+03},
	  {0.04,  7.3240E+02, 1.440E+03, 2.8539E+03, 3.6102E+03, 4.2226E+03, 4.7125E+03, 5.1314E+03, 5.4945E+03, 5.7678E+03, 6.0040E+03, 6.2244E+03, 6.3962E+03, 6.5852E+03, 6.7700E+03, 6.9867E+03, 7.1761E+03, 7.3984E+03, 7.5204E+03},
	  {0.05,  7.7680E+02, 1.593E+03, 3.1164E+03, 3.9967E+03, 4.7242E+03, 5.3248E+03, 5.8371E+03, 6.2868E+03, 6.6342E+03, 6.9338E+03, 7.2203E+03, 7.4137E+03, 7.6302E+03, 7.8376E+03, 8.0740E+03, 8.2871E+03, 8.5209E+03, 8.6525E+03},
	  {0.06,  8.0500E+02, 1.723E+03, 3.3203E+03, 4.3169E+03, 5.1543E+03, 5.8601E+03, 6.4632E+03, 6.9978E+03, 7.4196E+03, 7.7852E+03, 8.1435E+03, 8.3663E+03, 8.6193E+03, 8.8564E+03, 9.1192E+03, 9.3600E+03, 9.6097E+03, 9.7528E+03},
	  {0.07,  8.2050E+02, 1.833E+03, 3.4756E+03, 4.5788E+03, 5.5219E+03, 6.3287E+03, 7.0210E+03, 7.6391E+03, 8.1343E+03, 8.5662E+03, 8.9982E+03, 9.2554E+03, 9.5510E+03, 9.8229E+03, 1.0117E+04, 1.0390E+04, 1.0661E+04, 1.0820E+04},
	  {0.08,  8.2600E+02, 1.925E+03, 3.5914E+03, 4.7897E+03, 5.8323E+03, 6.7363E+03, 7.5168E+03, 8.2175E+03, 8.7858E+03, 9.2840E+03, 9.7906E+03, 1.0085E+04, 1.0426E+04, 1.0736E+04, 1.1065E+04, 1.1373E+04, 1.1669E+04, 1.1849E+04},
	  {0.09,  8.2390E+02, 2.004E+03, 3.6755E+03, 4.9568E+03, 6.0911E+03, 7.0875E+03, 7.9545E+03, 8.7372E+03, 9.3786E+03, 9.9435E+03, 1.0526E+04, 1.0859E+04, 1.1248E+04, 1.1597E+04, 1.1964E+04, 1.2308E+04, 1.2632E+04, 1.2836E+04},
	  {0.1,   8.1610E+02, 2.069E+03, 3.7347E+03, 5.0872E+03, 6.3043E+03, 7.3872E+03, 8.3378E+03, 9.2010E+03, 9.9152E+03, 1.0547E+04, 1.1206E+04, 1.1581E+04, 1.2018E+04, 1.2409E+04, 1.2815E+04, 1.3196E+04, 1.3551E+04, 1.3780E+04},
	  {0.15,  7.3710E+02, 2.245E+03, 3.8125E+03, 5.3870E+03, 6.8888E+03, 8.2986E+03, 9.5943E+03, 1.0808E+04, 1.1857E+04, 1.2813E+04, 1.3848E+04, 1.4455E+04, 1.5152E+04, 1.5773E+04, 1.6396E+04, 1.6986E+04, 1.7518E+04, 1.7904E+04},
	  {0.2,   6.6130E+02, 2.260E+03, 3.7349E+03, 5.4028E+03, 7.0451E+03, 8.6350E+03, 1.0145E+04, 1.1596E+04, 1.2890E+04, 1.4099E+04, 1.5453E+04, 1.6279E+04, 1.7230E+04, 1.8090E+04, 1.8945E+04, 1.9762E+04, 2.0497E+04, 2.1070E+04},
	  {0.25,  6.0060E+02, 2.193E+03, 3.6134E+03, 5.3185E+03, 7.0309E+03, 8.7189E+03, 1.0356E+04, 1.1955E+04, 1.3408E+04, 1.4791E+04, 1.6384E+04, 1.7373E+04, 1.8531E+04, 1.9594E+04, 2.0657E+04, 2.1683E+04, 2.2615E+04, 2.3375E+04},
	  {0.3,   5.5040E+02, 2.083E+03, 3.4818E+03, 5.1971E+03, 6.9445E+03, 8.6879E+03, 1.0402E+04, 1.2096E+04, 1.3652E+04, 1.5151E+04, 1.6916E+04, 1.8018E+04, 1.9330E+04, 2.0549E+04, 2.1779E+04, 2.2976E+04, 2.4076E+04, 2.4999E+04},
	  {0.4,   4.7190E+02, 1.839E+03, 3.2258E+03, 4.9243E+03, 6.6925E+03, 8.4850E+03, 1.0278E+04, 1.2077E+04, 1.3749E+04, 1.5382E+04, 1.7369E+04, 1.8598E+04, 2.0110E+04, 2.1535E+04, 2.2993E+04, 2.4431E+04, 2.5769E+04, 2.6928E+04},
	  {0.5,   4.1320E+02, 1.625E+03, 2.9949E+03, 4.6549E+03, 6.4129E+03, 8.2162E+03, 1.0041E+04, 1.1890E+04, 1.3620E+04, 1.5321E+04, 1.7442E+04, 1.8727E+04, 2.0354E+04, 2.1900E+04, 2.3502E+04, 2.5091E+04, 2.6580E+04, 2.7889E+04},
	  {0.6,   3.6800E+02, 1.452E+03, 2.7909E+03, 4.4036E+03, 6.1372E+03, 7.9331E+03, 9.7677E+03, 1.1639E+04, 1.3398E+04, 1.5136E+04, 1.7344E+04, 1.8654E+04, 2.0352E+04, 2.1975E+04, 2.3672E+04, 2.5363E+04, 2.6953E+04, 2.8361E+04},
	  {0.7,   3.3250E+02, 1.315E+03, 2.6110E+03, 4.1732E+03, 5.8747E+03, 7.6534E+03, 9.4845E+03, 1.1364E+04, 1.3136E+04, 1.4895E+04, 1.7160E+04, 1.8479E+04, 2.0224E+04, 2.1896E+04, 2.3659E+04, 2.5421E+04, 2.7082E+04, 2.8559E+04},
	  {0.8,   3.0390E+02, 1.204E+03, 2.4517E+03, 3.9629E+03, 5.6283E+03, 7.3839E+03, 9.2035E+03, 1.1081E+04, 1.2857E+04, 1.4626E+04, 1.6930E+04, 1.8250E+04, 2.0024E+04, 2.1730E+04, 2.3540E+04, 2.5354E+04, 2.7066E+04, 2.8592E+04},
	  {0.9,   2.8050E+02, 1.112E+03, 2.3103E+03, 3.7710E+03, 5.3983E+03, 7.1272E+03, 8.9301E+03, 1.0799E+04, 1.2573E+04, 1.4345E+04, 1.6675E+04, 1.7990E+04, 1.9785E+04, 2.1513E+04, 2.3356E+04, 2.5208E+04, 2.6958E+04, 2.8521E+04},
	  {1,     2.6080E+02, 1.035E+03, 2.1841E+03, 3.5957E+03, 5.1841E+03, 6.8840E+03, 8.6668E+03, 1.0523E+04, 1.2291E+04, 1.4061E+04, 1.6407E+04, 1.7716E+04, 1.9521E+04, 2.1265E+04, 2.3132E+04, 2.5013E+04, 2.6791E+04, 2.8381E+04},
	  {1.5,   1.9570E+02, 7.777E+02, 1.7151E+03, 2.9117E+03, 4.3121E+03, 5.8573E+03, 7.5173E+03, 9.2787E+03, 1.0982E+04, 1.2705E+04, 1.5045E+04, 1.6313E+04, 1.8120E+04, 1.9878E+04, 2.1792E+04, 2.3734E+04, 2.5579E+04, 2.7228E+04},
	  {2,     1.5860E+02, 6.306E+02, 1.4139E+03, 2.4439E+03, 3.6826E+03, 5.0814E+03, 6.6126E+03, 8.2615E+03, 9.8800E+03, 1.1530E+04, 1.3799E+04, 1.5026E+04, 1.6795E+04, 1.8525E+04, 2.0425E+04, 2.2366E+04, 2.4214E+04, 2.5869E+04},
	  {2.5,   1.3440E+02, 5.344E+02, 1.2053E+03, 2.1062E+03, 3.2109E+03, 4.4808E+03, 5.8919E+03, 7.4307E+03, 8.9601E+03, 1.0532E+04, 1.2706E+04, 1.3895E+04, 1.5611E+04, 1.7298E+04, 1.9162E+04, 2.1074E+04, 2.2901E+04, 2.4541E+04},
	  {3,     1.1720E+02, 4.659E+02, 1.0525E+03, 1.8518E+03, 2.8459E+03, 4.0044E+03, 5.3080E+03, 6.7431E+03, 8.1871E+03, 9.6823E+03, 1.1753E+04, 1.2907E+04, 1.4565E+04, 1.6202E+04, 1.8020E+04, 1.9891E+04, 2.1685E+04, 2.3299E+04},
	  {4,     9.4040E+01, 3.739E+02, 8.4417E+02, 1.4956E+03, 2.3203E+03, 3.3005E+03, 4.4226E+03, 5.6787E+03, 6.9687E+03, 8.3208E+03, 1.0187E+04, 1.1277E+04, 1.2819E+04, 1.4354E+04, 1.6066E+04, 1.7840E+04, 1.9553E+04, 2.1103E+04},
	  {5,     7.9110E+01, 3.146E+02, 7.0862E+02, 1.2587E+03, 1.9619E+03, 2.8080E+03, 3.7883E+03, 4.8981E+03, 6.0574E+03, 7.2846E+03, 8.9672E+03, 9.9981E+03, 1.1431E+04, 1.2867E+04, 1.4472E+04, 1.6146E+04, 1.7771E+04, 1.9255E+04},
	  {6,     6.8580E+01, 2.727E+02, 6.1317E+02, 1.0901E+03, 1.7028E+03, 2.4458E+03, 3.3138E+03, 4.3045E+03, 5.3535E+03, 6.4735E+03, 7.9956E+03, 8.9727E+03, 1.0305E+04, 1.1651E+04, 1.3155E+04, 1.4731E+04, 1.6271E+04, 1.7686E+04},
	  {7,     6.0710E+01, 2.414E+02, 5.4214E+02, 9.6393E+02, 1.5072E+03, 2.1691E+03, 2.9467E+03, 3.8393E+03, 4.7954E+03, 5.8234E+03, 7.2072E+03, 8.1344E+03, 9.3772E+03, 1.0640E+04, 1.2051E+04, 1.3536E+04, 1.4996E+04, 1.6346E+04},
	  {8,     5.4600E+01, 2.171E+02, 4.8708E+02, 8.6589E+02, 1.3543E+03, 1.9511E+03, 2.6550E+03, 3.4663E+03, 4.3434E+03, 5.2919E+03, 6.5562E+03, 7.4376E+03, 8.6003E+03, 9.7876E+03, 1.1114E+04, 1.2516E+04, 1.3900E+04, 1.5190E+04},
	  {9,     4.9690E+01, 1.976E+02, 4.4305E+02, 7.8742E+02, 1.2315E+03, 1.7751E+03, 2.4179E+03, 3.1610E+03, 3.9704E+03, 4.8502E+03, 6.0112E+03, 6.8507E+03, 7.9413E+03, 9.0607E+03, 1.0311E+04, 1.1636E+04, 1.2949E+04, 1.4183E+04},
	  {10,    4.5670E+01, 1.816E+02, 4.0697E+02, 7.2313E+02, 1.1307E+03, 1.6302E+03, 2.2217E+03, 2.9069E+03, 3.6580E+03, 4.4779E+03, 5.5493E+03, 6.3500E+03, 7.3761E+03, 8.4341E+03, 9.6149E+03, 1.0869E+04, 1.2118E+04, 1.3300E+04},
	  {15,    3.2920E+01, 1.309E+02, 2.9312E+02, 5.2053E+02, 8.1305E+02, 1.1714E+03, 1.5965E+03, 2.0903E+03, 2.6420E+03, 3.2520E+03, 4.0154E+03, 4.6625E+03, 5.4500E+03, 6.2756E+03, 7.1917E+03, 8.1727E+03, 9.1688E+03, 1.0137E+04},
	  {20,    2.6070E+01, 1.037E+02, 2.3208E+02, 4.1214E+02, 6.4344E+02, 9.2630E+02, 1.2614E+03, 1.6501E+03, 2.0878E+03, 2.5745E+03, 3.1635E+03, 3.7049E+03, 4.3420E+03, 5.0169E+03, 5.7614E+03, 6.5617E+03, 7.3850E+03, 8.2021E+03},
	  {25,    2.1750E+01, 8.649E+01, 1.9364E+02, 3.4394E+02, 5.3693E+02, 7.7269E+02, 1.0516E+03, 1.3745E+03, 1.7393E+03, 2.1460E+03, 2.6262E+03, 3.0916E+03, 3.6272E+03, 4.1981E+03, 4.8249E+03, 5.4997E+03, 6.2003E+03, 6.9062E+03},
	  {30,    1.8760E+01, 7.461E+01, 1.6706E+02, 2.9680E+02, 4.6337E+02, 6.6676E+02, 9.0715E+02, 1.1851E+03, 1.4994E+03, 1.8503E+03, 2.2573E+03, 2.6660E+03, 3.1290E+03, 3.6247E+03, 4.1667E+03, 4.7505E+03, 5.3604E+03, 5.9819E+03},
	  {40,    1.4880E+01, 5.917E+01, 1.3252E+02, 2.3550E+02, 3.6777E+02, 5.2925E+02, 7.1995E+02, 9.4004E+02, 1.1892E+03, 1.4675E+03, 1.7832E+03, 2.1135E+03, 2.4808E+03, 2.8758E+03, 3.3050E+03, 3.7671E+03, 4.2535E+03, 4.7570E+03},
	  {50,    1.2450E+01, 4.952E+01, 1.1092E+02, 1.9717E+02, 3.0797E+02, 4.4328E+02, 6.0305E+02, 7.8733E+02, 9.9601E+02, 1.2291E+03, 1.4904E+03, 1.7695E+03, 2.0767E+03, 2.4078E+03, 2.7661E+03, 3.1514E+03, 3.5588E+03, 3.9841E+03},
	  {60,    1.0780E+01, 4.289E+01, 9.6080E+01, 1.7081E+02, 2.6684E+02, 3.8415E+02, 5.2268E+02, 6.8244E+02, 8.6336E+02, 1.0654E+03, 1.2907E+03, 1.5336E+03, 1.7997E+03, 2.0868E+03, 2.3966E+03, 2.7293E+03, 3.0821E+03, 3.4523E+03},
	  {70,    9.5590E+00, 3.803E+01, 8.5220E+01, 1.5152E+02, 2.3674E+02, 3.4086E+02, 4.6384E+02, 6.0568E+02, 7.6632E+02, 9.4575E+02, 1.1451E+03, 1.3613E+03, 1.5975E+03, 1.8523E+03, 2.1269E+03, 2.4215E+03, 2.7343E+03, 3.0636E+03},
	  {80,    8.6250E+00, 3.432E+01, 7.6915E+01, 1.3676E+02, 2.1371E+02, 3.0773E+02, 4.1882E+02, 5.4694E+02, 6.9207E+02, 8.5419E+02, 1.0339E+03, 1.2296E+03, 1.4429E+03, 1.6731E+03, 1.9215E+03, 2.1866E+03, 2.4689E+03, 2.7668E+03},
	  {90,    7.8880E+00, 3.139E+01, 7.0350E+01, 1.2510E+02, 1.9550E+02, 2.8154E+02, 3.8320E+02, 5.0048E+02, 6.3334E+02, 7.8178E+02, 9.4615E+02, 1.1255E+03, 1.3207E+03, 1.5315E+03, 1.7603E+03, 2.0012E+03, 2.2595E+03, 2.5324E+03},
	  {100,   7.2890E+00, 2.901E+01, 6.5026E+01, 1.1564E+02, 1.8072E+02, 2.6028E+02, 3.5431E+02, 4.6278E+02, 5.8568E+02, 7.2301E+02, 8.7498E+02, 1.0410E+03, 1.2216E+03, 1.4167E+03, 1.6263E+03, 1.8509E+03, 2.0898E+03, 2.3425E+03},
	  {150,   5.4450E+00, 2.168E+01, 4.8615E+01, 8.6471E+01, 1.3517E+02, 1.9473E+02, 2.6515E+02, 3.4643E+02, 4.3857E+02, 5.4158E+02, 6.5548E+02, 7.8020E+02, 9.1581E+02, 1.0623E+03, 1.2197E+03, 1.3881E+03, 1.5673E+03, 1.7573E+03},
	  {200,   4.4920E+00, 1.789E+01, 4.0137E+01, 7.1399E+01, 1.1163E+02, 1.6083E+02, 2.1903E+02, 2.8622E+02, 3.6242E+02, 4.4763E+02, 5.4186E+02, 6.4511E+02, 7.5739E+02, 8.7870E+02, 1.0091E+03, 1.1485E+03, 1.2970E+03, 1.4546E+03},
	  {250,   3.9110E+00, 1.558E+01, 3.4964E+01, 6.2202E+01, 9.7256E+01, 1.4014E+02, 1.9087E+02, 2.4945E+02, 3.1590E+02, 3.9022E+02, 4.7243E+02, 5.6252E+02, 6.6053E+02, 7.6644E+02, 8.8027E+02, 1.0020E+03, 1.1318E+03, 1.2694E+03},
	  {300,   3.5200E+00, 1.410E+01, 3.1490E+01, 5.6023E+01, 8.7601E+01, 1.2624E+02, 1.7194E+02, 2.2474E+02, 2.8462E+02, 3.5162E+02, 4.2574E+02, 5.0698E+02, 5.9537E+02, 6.9091E+02, 7.9361E+02, 9.0348E+02, 1.0205E+03, 1.1448E+03},
	  {400,   3.0320E+00, 1.210E+01, 2.7148E+01, 4.8303E+01, 7.5535E+01, 1.0886E+02, 1.4829E+02, 1.9384E+02, 2.4552E+02, 3.0335E+02, 3.6734E+02, 4.3749E+02, 5.1383E+02, 5.9637E+02, 6.8512E+02, 7.8008E+02, 8.8128E+02, 9.8872E+02},
	  {500,   2.7430E+00, 1.100E+01, 2.4579E+01, 4.3734E+01, 6.8395E+01, 9.8575E+01, 1.3429E+02, 1.7555E+02, 2.2237E+02, 2.7477E+02, 3.3275E+02, 3.9634E+02, 4.6554E+02, 5.4036E+02, 6.2083E+02, 7.0695E+02, 7.9874E+02, 8.9620E+02},
	  {600,   2.5560E+00, 1.020E+01, 2.2911E+01, 4.0767E+01, 6.3757E+01, 9.1894E+01, 1.2519E+02, 1.6367E+02, 2.0733E+02, 2.5619E+02, 3.1027E+02, 3.6958E+02, 4.3413E+02, 5.0394E+02, 5.7902E+02, 6.5938E+02, 7.4504E+02, 8.3601E+02},
	  {700,   2.4260E+00, 9.700E+00, 2.1761E+01, 3.8723E+01, 6.0561E+01, 8.7290E+01, 1.1892E+02, 1.5547E+02, 1.9696E+02, 2.4339E+02, 2.9478E+02, 3.5114E+02, 4.1248E+02, 4.7883E+02, 5.5019E+02, 6.2657E+02, 7.0800E+02, 7.9448E+02},
	  {800,   2.3330E+00, 9.330E+00, 2.0937E+01, 3.7256E+01, 5.8268E+01, 8.3986E+01, 1.1442E+02, 1.4959E+02, 1.8951E+02, 2.3419E+02, 2.8364E+02, 3.3789E+02, 3.9693E+02, 4.6078E+02, 5.2947E+02, 6.0299E+02, 6.8138E+02, 7.6464E+02},
	  {900,   2.2640E+00, 9.040E+00, 2.0327E+01, 3.6172E+01, 5.6574E+01, 8.1545E+01, 1.1110E+02, 1.4525E+02, 1.8401E+02, 2.2740E+02, 2.7542E+02, 3.2810E+02, 3.8544E+02, 4.4746E+02, 5.1417E+02, 5.8558E+02, 6.6172E+02, 7.4259E+02},
	  {1000,  2.2110E+00, 8.810E+00, 1.9862E+01, 3.5345E+01, 5.5280E+01, 7.9682E+01, 1.0856E+02, 1.4194E+02, 1.7981E+02, 2.2222E+02, 2.6915E+02, 3.2063E+02, 3.7667E+02, 4.3728E+02, 5.0248E+02, 5.7229E+02, 6.4671E+02, 7.2575E+02}
  }
};

/** ----------------------------------------------- Static objects  --------------------------------------- */


static const AT_stopping_power_sources_struct AT_stopping_power_sources = {
	STOPPING_POWER_SOURCE_N,
    {  PSTAR,         Bethe, 			ShieldHit           }, // source_no
    {  "PSTAR data",  "Bethe formula",  "ShieldHit (Bethe)" }  // source_name
};


static const AT_stopping_power_tabulated_source_group_for_all_materials_struct AT_data_PSTAR_source = {
  MATERIAL_DATA_N,
  PSTAR,
  { User_Defined_Material,  Water_Liquid,                         Aluminum_Oxide,                                 Aluminum,                                 PMMA,                                Alanine,                                LiF,                                Air },
  { NULL,                   &AT_stopping_power_data_PSTAR_Water,  &AT_stopping_power_data_PSTAR_Aluminum_Oxide,   &AT_stopping_power_data_PSTAR_Aluminum,   &AT_stopping_power_data_PSTAR_PMMA,  &AT_stopping_power_data_PSTAR_Alanine,  &AT_stopping_power_data_PSTAR_LiF,  &AT_stopping_power_data_PSTAR_Air}
};


static const AT_stopping_power_tabulated_source_group_struct AT_stopping_power_tabulated_source = {
  STOPPING_POWER_SOURCE_N,
  { PSTAR,                  Bethe, 	ShieldHit},
  { &AT_data_PSTAR_source,  NULL,	NULL}
};


static const AT_stopping_power_analytical_sources_struct AT_stopping_power_analytical_source = {
  STOPPING_POWER_SOURCE_N,
  { PSTAR,   Bethe,				ShieldHit},
  { NULL,    &AT_Bethe_wrapper, &AT_ShieldHit_wrapper}
};

#endif /* AT_DATASTOPPINGPOWER_H_ */
