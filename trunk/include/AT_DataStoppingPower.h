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
  ShieldHit			   = 2, /**< ShieldHit code: extended Bethe formula */
  ICRU                 = 3  /**< ICRU 49 and 73: for liquid water */
};


#define BETHE_LOWER_LIMIT_E_MEV_U 1.0

/**
 * TODO
 */
#define STOPPING_POWER_SOURCE_N    4


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
 * @return statuc code
 */
int AT_stopping_power_source_model_name_from_number( const long source_no,
		char* source_name);


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
 * Wrapper for ICRU stopping powers
 * momentarily read from hardwired struct for water
 * TODO: augment by routine to read data file into
 * TODO: struct and look-up value
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return stopping power (MeV cm2 per g)
 */
double AT_ICRU_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no);

/**
 * TODO
 */
double _AT_Stopping_Power_get_data(const long stopping_power_source_no,
		const long 	    particle_no,
		const long 		material_no,
		AT_stopping_power_tabulated_source_for_given_material_struct** source_for_given_material);

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

static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_Silicon = {
  132,
  PSTAR,
  Silicon,
  {
		  { 1.000e-03 , 1.169e+02 },
		  { 1.500e-03 , 1.393e+02 },
		  { 2.000e-03 , 1.583e+02 },
		  { 2.500e-03 , 1.750e+02 },
		  { 3.000e-03 , 1.903e+02 },
		  { 4.000e-03 , 2.174e+02 },
		  { 5.000e-03 , 2.414e+02 },
		  { 6.000e-03 , 2.633e+02 },
		  { 7.000e-03 , 2.834e+02 },
		  { 8.000e-03 , 3.022e+02 },
		  { 9.000e-03 , 3.199e+02 },
		  { 1.000e-02 , 3.366e+02 },
		  { 1.250e-02 , 3.686e+02 },
		  { 1.500e-02 , 3.961e+02 },
		  { 1.750e-02 , 4.199e+02 },
		  { 2.000e-02 , 4.408e+02 },
		  { 2.250e-02 , 4.591e+02 },
		  { 2.500e-02 , 4.751e+02 },
		  { 2.750e-02 , 4.890e+02 },
		  { 3.000e-02 , 5.011e+02 },
		  { 3.500e-02 , 5.204e+02 },
		  { 4.000e-02 , 5.342e+02 },
		  { 4.500e-02 , 5.433e+02 },
		  { 5.000e-02 , 5.485e+02 },
		  { 5.500e-02 , 5.507e+02 },
		  { 6.000e-02 , 5.504e+02 },
		  { 6.500e-02 , 5.480e+02 },
		  { 7.000e-02 , 5.442e+02 },
		  { 7.500e-02 , 5.392e+02 },
		  { 8.000e-02 , 5.333e+02 },
		  { 8.500e-02 , 5.267e+02 },
		  { 9.000e-02 , 5.198e+02 },
		  { 9.500e-02 , 5.126e+02 },
		  { 1.000e-01 , 5.053e+02 },
		  { 1.250e-01 , 4.695e+02 },
		  { 1.500e-01 , 4.379e+02 },
		  { 1.750e-01 , 4.113e+02 },
		  { 2.000e-01 , 3.889e+02 },
		  { 2.250e-01 , 3.699e+02 },
		  { 2.500e-01 , 3.535e+02 },
		  { 2.750e-01 , 3.391e+02 },
		  { 3.000e-01 , 3.263e+02 },
		  { 3.500e-01 , 3.045e+02 },
		  { 4.000e-01 , 2.862e+02 },
		  { 4.500e-01 , 2.706e+02 },
		  { 5.000e-01 , 2.569e+02 },
		  { 5.500e-01 , 2.454e+02 },
		  { 6.000e-01 , 2.353e+02 },
		  { 6.500e-01 , 2.261e+02 },
		  { 7.000e-01 , 2.174e+02 },
		  { 7.500e-01 , 2.091e+02 },
		  { 8.000e-01 , 2.012e+02 },
		  { 8.500e-01 , 1.938e+02 },
		  { 9.000e-01 , 1.872e+02 },
		  { 9.500e-01 , 1.811e+02 },
		  { 1.000e+00 , 1.754e+02 },
		  { 1.250e+00 , 1.524e+02 },
		  { 1.500e+00 , 1.355e+02 },
		  { 1.750e+00 , 1.223e+02 },
		  { 2.000e+00 , 1.118e+02 },
		  { 2.250e+00 , 1.031e+02 },
		  { 2.500e+00 , 9.586e+01 },
		  { 2.750e+00 , 8.968e+01 },
		  { 3.000e+00 , 8.433e+01 },
		  { 3.500e+00 , 7.554e+01 },
		  { 4.000e+00 , 6.860e+01 },
		  { 4.500e+00 , 6.296e+01 },
		  { 5.000e+00 , 5.828e+01 },
		  { 5.500e+00 , 5.431e+01 },
		  { 6.000e+00 , 5.091e+01 },
		  { 6.500e+00 , 4.795e+01 },
		  { 7.000e+00 , 4.536e+01 },
		  { 7.500e+00 , 4.306e+01 },
		  { 8.000e+00 , 4.101e+01 },
		  { 8.500e+00 , 3.916e+01 },
		  { 9.000e+00 , 3.750e+01 },
		  { 9.500e+00 , 3.598e+01 },
		  { 1.000e+01 , 3.459e+01 },
		  { 1.250e+01 , 2.913e+01 },
		  { 1.500e+01 , 2.529e+01 },
		  { 1.750e+01 , 2.243e+01 },
		  { 2.000e+01 , 2.020e+01 },
		  { 2.500e+01 , 1.695e+01 },
		  { 2.750e+01 , 1.573e+01 },
		  { 3.000e+01 , 1.469e+01 },
		  { 3.500e+01 , 1.302e+01 },
		  { 4.000e+01 , 1.173e+01 },
		  { 4.500e+01 , 1.070e+01 },
		  { 5.000e+01 , 9.856e+00 },
		  { 5.500e+01 , 9.157e+00 },
		  { 6.000e+01 , 8.564e+00 },
		  { 6.500e+01 , 8.056e+00 },
		  { 7.000e+01 , 7.614e+00 },
		  { 7.500e+01 , 7.227e+00 },
		  { 8.000e+01 , 6.885e+00 },
		  { 8.500e+01 , 6.581e+00 },
		  { 9.000e+01 , 6.308e+00 },
		  { 9.500e+01 , 6.061e+00 },
		  { 1.000e+02 , 5.838e+00 },
		  { 1.250e+02 , 4.974e+00 },
		  { 1.500e+02 , 4.383e+00 },
		  { 1.750e+02 , 3.954e+00 },
		  { 2.000e+02 , 3.628e+00 },
		  { 2.250e+02 , 3.371e+00 },
		  { 2.500e+02 , 3.165e+00 },
		  { 2.750e+02 , 2.995e+00 },
		  { 3.000e+02 , 2.853e+00 },
		  { 3.500e+02 , 2.630e+00 },
		  { 4.000e+02 , 2.462e+00 },
		  { 4.500e+02 , 2.333e+00 },
		  { 5.000e+02 , 2.230e+00 },
		  { 5.500e+02 , 2.147e+00 },
		  { 6.000e+02 , 2.079e+00 },
		  { 6.500e+02 , 2.022e+00 },
		  { 7.000e+02 , 1.975e+00 },
		  { 7.500e+02 , 1.934e+00 },
		  { 8.000e+02 , 1.899e+00 },
		  { 8.500e+02 , 1.870e+00 },
		  { 9.000e+02 , 1.844e+00 },
		  { 9.500e+02 , 1.821e+00 },
		  { 1.000e+03 , 1.801e+00 },
		  { 1.500e+03 , 1.696e+00 },
		  { 2.000e+03 , 1.666e+00 },
		  { 2.500e+03 , 1.661e+00 },
		  { 3.000e+03 , 1.668e+00 },
		  { 4.000e+03 , 1.692e+00 },
		  { 5.000e+03 , 1.720e+00 },
		  { 6.000e+03 , 1.746e+00 },
		  { 7.000e+03 , 1.770e+00 },
		  { 8.000e+03 , 1.792e+00 },
		  { 9.000e+03 , 1.812e+00 },
		  { 1.000e+04 , 1.830e+00 }
  }
};

static const AT_stopping_power_tabulated_source_for_given_material_struct AT_stopping_power_data_PSTAR_Copper = {
  132,
  PSTAR,
  Copper,
  {
		  { 1.000e-03 , 3.931e+01 },
		  { 1.500e-03 , 4.697e+01 },
		  { 2.000e-03 , 5.340e+01 },
		  { 2.500e-03 , 5.906e+01 },
		  { 3.000e-03 , 6.416e+01 },
		  { 4.000e-03 , 7.324e+01 },
		  { 5.000e-03 , 8.124e+01 },
		  { 6.000e-03 , 8.849e+01 },
		  { 7.000e-03 , 9.517e+01 },
		  { 8.000e-03 , 1.014e+02 },
		  { 9.000e-03 , 1.073e+02 },
		  { 1.000e-02 , 1.128e+02 },
		  { 1.250e-02 , 1.236e+02 },
		  { 1.500e-02 , 1.329e+02 },
		  { 1.750e-02 , 1.411e+02 },
		  { 2.000e-02 , 1.484e+02 },
		  { 2.250e-02 , 1.548e+02 },
		  { 2.500e-02 , 1.606e+02 },
		  { 2.750e-02 , 1.657e+02 },
		  { 3.000e-02 , 1.703e+02 },
		  { 3.500e-02 , 1.781e+02 },
		  { 4.000e-02 , 1.844e+02 },
		  { 4.500e-02 , 1.894e+02 },
		  { 5.000e-02 , 1.935e+02 },
		  { 5.500e-02 , 1.968e+02 },
		  { 6.000e-02 , 1.996e+02 },
		  { 6.500e-02 , 2.018e+02 },
		  { 7.000e-02 , 2.036e+02 },
		  { 7.500e-02 , 2.051e+02 },
		  { 8.000e-02 , 2.063e+02 },
		  { 8.500e-02 , 2.073e+02 },
		  { 9.000e-02 , 2.082e+02 },
		  { 9.500e-02 , 2.088e+02 },
		  { 1.000e-01 , 2.093e+02 },
		  { 1.250e-01 , 2.104e+02 },
		  { 1.500e-01 , 2.097e+02 },
		  { 1.750e-01 , 2.079e+02 },
		  { 2.000e-01 , 2.055e+02 },
		  { 2.250e-01 , 2.026e+02 },
		  { 2.500e-01 , 1.994e+02 },
		  { 2.750e-01 , 1.961e+02 },
		  { 3.000e-01 , 1.926e+02 },
		  { 3.500e-01 , 1.856e+02 },
		  { 4.000e-01 , 1.788e+02 },
		  { 4.500e-01 , 1.722e+02 },
		  { 5.000e-01 , 1.660e+02 },
		  { 5.500e-01 , 1.598e+02 },
		  { 6.000e-01 , 1.538e+02 },
		  { 6.500e-01 , 1.481e+02 },
		  { 7.000e-01 , 1.428e+02 },
		  { 7.500e-01 , 1.379e+02 },
		  { 8.000e-01 , 1.334e+02 },
		  { 8.500e-01 , 1.292e+02 },
		  { 9.000e-01 , 1.253e+02 },
		  { 9.500e-01 , 1.218e+02 },
		  { 1.000e+00 , 1.184e+02 },
		  { 1.250e+00 , 1.048e+02 },
		  { 1.500e+00 , 9.450e+01 },
		  { 1.750e+00 , 8.652e+01 },
		  { 2.000e+00 , 7.992e+01 },
		  { 2.250e+00 , 7.430e+01 },
		  { 2.500e+00 , 6.960e+01 },
		  { 2.750e+00 , 6.556e+01 },
		  { 3.000e+00 , 6.199e+01 },
		  { 3.500e+00 , 5.610e+01 },
		  { 4.000e+00 , 5.136e+01 },
		  { 4.500e+00 , 4.746e+01 },
		  { 5.000e+00 , 4.418e+01 },
		  { 5.500e+00 , 4.139e+01 },
		  { 6.000e+00 , 3.897e+01 },
		  { 6.500e+00 , 3.685e+01 },
		  { 7.000e+00 , 3.498e+01 },
		  { 7.500e+00 , 3.331e+01 },
		  { 8.000e+00 , 3.181e+01 },
		  { 8.500e+00 , 3.046e+01 },
		  { 9.000e+00 , 2.924e+01 },
		  { 9.500e+00 , 2.812e+01 },
		  { 1.000e+01 , 2.709e+01 },
		  { 1.250e+01 , 2.301e+01 },
		  { 1.500e+01 , 2.010e+01 },
		  { 1.750e+01 , 1.791e+01 },
		  { 2.000e+01 , 1.620e+01 },
		  { 2.500e+01 , 1.368e+01 },
		  { 2.750e+01 , 1.273e+01 },
		  { 3.000e+01 , 1.191e+01 },
		  { 3.500e+01 , 1.060e+01 },
		  { 4.000e+01 , 9.577e+00 },
		  { 4.500e+01 , 8.760e+00 },
		  { 5.000e+01 , 8.090e+00 },
		  { 5.500e+01 , 7.530e+00 },
		  { 6.000e+01 , 7.055e+00 },
		  { 6.500e+01 , 6.646e+00 },
		  { 7.000e+01 , 6.291e+00 },
		  { 7.500e+01 , 5.979e+00 },
		  { 8.000e+01 , 5.702e+00 },
		  { 8.500e+01 , 5.455e+00 },
		  { 9.000e+01 , 5.234e+00 },
		  { 9.500e+01 , 5.034e+00 },
		  { 1.000e+02 , 4.852e+00 },
		  { 1.250e+02 , 4.147e+00 },
		  { 1.500e+02 , 3.664e+00 },
		  { 1.750e+02 , 3.311e+00 },
		  { 2.000e+02 , 3.042e+00 },
		  { 2.250e+02 , 2.830e+00 },
		  { 2.500e+02 , 2.659e+00 },
		  { 2.750e+02 , 2.518e+00 },
		  { 3.000e+02 , 2.400e+00 },
		  { 3.500e+02 , 2.214e+00 },
		  { 4.000e+02 , 2.074e+00 },
		  { 4.500e+02 , 1.966e+00 },
		  { 5.000e+02 , 1.880e+00 },
		  { 5.500e+02 , 1.811e+00 },
		  { 6.000e+02 , 1.754e+00 },
		  { 6.500e+02 , 1.706e+00 },
		  { 7.000e+02 , 1.666e+00 },
		  { 7.500e+02 , 1.632e+00 },
		  { 8.000e+02 , 1.603e+00 },
		  { 8.500e+02 , 1.578e+00 },
		  { 9.000e+02 , 1.556e+00 },
		  { 9.500e+02 , 1.538e+00 },
		  { 1.000e+03 , 1.521e+00 },
		  { 1.500e+03 , 1.433e+00 },
		  { 2.000e+03 , 1.407e+00 },
		  { 2.500e+03 , 1.404e+00 },
		  { 3.000e+03 , 1.409e+00 },
		  { 4.000e+03 , 1.430e+00 },
		  { 5.000e+03 , 1.454e+00 },
		  { 6.000e+03 , 1.476e+00 },
		  { 7.000e+03 , 1.498e+00 },
		  { 8.000e+03 , 1.517e+00 },
		  { 9.000e+03 , 1.535e+00 },
		  { 1.000e+04 , 1.552e+00 }
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
  const double     E_MeV_u_and_stopping_power_total_MeV_cm2_g[19][53];
} AT_stopping_power_ShieldHit_table_struct;

static const AT_stopping_power_ShieldHit_table_struct AT_stopping_power_ShieldHit_table = {
  Water_Liquid,
  53,
  {
			{0.025000,  0.030000,  0.040000,  0.050000,  0.060000,  0.070000,  0.080000,  0.090000,  0.100000,  0.150000,  0.200000,  0.250000,  0.300000,  0.400000,  0.500000,  0.600000,  0.700000,  0.800000,  0.900000,  1.000000,  1.500000,  2.000000,  2.500000,  3.000000,  4.000000,  5.000000,  6.000000,  7.000000,  8.000000,  9.000000,  10.000000, 15.000000, 20.000000, 25.000000, 30.000000, 40.000000, 50.000000, 60.000000, 70.000000, 80.000000, 90.000000, 100.000000, 150.000000, 200.000000, 250.000000, 300.000000, 400.000000, 500.000000, 600.000000, 700.000000, 800.000000, 900.000000, 1000.000000},
			{6.245e+02, 6.671e+02, 7.324e+02, 7.768e+02, 8.050e+02, 8.205e+02, 8.260e+02, 8.239e+02, 8.161e+02, 7.371e+02, 6.613e+02, 6.006e+02, 5.504e+02, 4.719e+02, 4.132e+02, 3.680e+02, 3.325e+02, 3.039e+02, 2.805e+02, 2.608e+02, 1.957e+02, 1.586e+02, 1.344e+02, 1.172e+02, 9.404e+01, 7.911e+01, 6.858e+01, 6.071e+01, 5.460e+01, 4.969e+01, 4.567e+01, 3.292e+01, 2.607e+01, 2.175e+01, 1.876e+01, 1.488e+01, 1.245e+01, 1.078e+01, 9.559e+00, 8.625e+00, 7.888e+00, 7.289e+00, 5.445e+00, 4.492e+00, 3.911e+00, 3.520e+00, 3.032e+00, 2.743e+00, 2.556e+00, 2.426e+00, 2.333e+00, 2.264e+00, 2.211e+00},
			{1.151e+03, 1.257e+03, 1.440e+03, 1.593e+03, 1.723e+03, 1.833e+03, 1.925e+03, 2.004e+03, 2.069e+03, 2.245e+03, 2.260e+03, 2.193e+03, 2.083e+03, 1.839e+03, 1.625e+03, 1.452e+03, 1.315e+03, 1.204e+03, 1.112e+03, 1.035e+03, 7.777e+02, 6.306e+02, 5.344e+02, 4.659e+02, 3.739e+02, 3.146e+02, 2.727e+02, 2.414e+02, 2.171e+02, 1.976e+02, 1.816e+02, 1.309e+02, 1.037e+02, 8.649e+01, 7.461e+01, 5.917e+01, 4.952e+01, 4.289e+01, 3.803e+01, 3.432e+01, 3.139e+01, 2.901e+01, 2.168e+01, 1.789e+01, 1.558e+01, 1.410e+01, 1.210e+01, 1.100e+01, 1.020e+01, 9.700e+00, 9.330e+00, 9.040e+00, 8.810e+00},
			{2.319e+03, 2.520e+03, 2.854e+03, 3.116e+03, 3.320e+03, 3.476e+03, 3.591e+03, 3.675e+03, 3.735e+03, 3.812e+03, 3.735e+03, 3.613e+03, 3.482e+03, 3.226e+03, 2.995e+03, 2.791e+03, 2.611e+03, 2.452e+03, 2.310e+03, 2.184e+03, 1.715e+03, 1.414e+03, 1.205e+03, 1.052e+03, 8.442e+02, 7.086e+02, 6.132e+02, 5.421e+02, 4.871e+02, 4.431e+02, 4.070e+02, 2.931e+02, 2.321e+02, 1.936e+02, 1.671e+02, 1.325e+02, 1.109e+02, 9.608e+01, 8.522e+01, 7.692e+01, 7.035e+01, 6.503e+01, 4.862e+01, 4.014e+01, 3.496e+01, 3.149e+01, 2.715e+01, 2.458e+01, 2.291e+01, 2.176e+01, 2.094e+01, 2.033e+01, 1.986e+01},
			{2.872e+03, 3.144e+03, 3.610e+03, 3.997e+03, 4.317e+03, 4.579e+03, 4.790e+03, 4.957e+03, 5.087e+03, 5.387e+03, 5.403e+03, 5.318e+03, 5.197e+03, 4.924e+03, 4.655e+03, 4.404e+03, 4.173e+03, 3.963e+03, 3.771e+03, 3.596e+03, 2.912e+03, 2.444e+03, 2.106e+03, 1.852e+03, 1.496e+03, 1.259e+03, 1.090e+03, 9.639e+02, 8.659e+02, 7.874e+02, 7.231e+02, 5.205e+02, 4.121e+02, 3.439e+02, 2.968e+02, 2.355e+02, 1.972e+02, 1.708e+02, 1.515e+02, 1.368e+02, 1.251e+02, 1.156e+02, 8.647e+01, 7.140e+01, 6.220e+01, 5.602e+01, 4.830e+01, 4.373e+01, 4.077e+01, 3.872e+01, 3.726e+01, 3.617e+01, 3.534e+01},
			{3.292e+03, 3.631e+03, 4.223e+03, 4.724e+03, 5.154e+03, 5.522e+03, 5.832e+03, 6.091e+03, 6.304e+03, 6.889e+03, 7.045e+03, 7.031e+03, 6.944e+03, 6.692e+03, 6.413e+03, 6.137e+03, 5.875e+03, 5.628e+03, 5.398e+03, 5.184e+03, 4.312e+03, 3.683e+03, 3.211e+03, 2.846e+03, 2.320e+03, 1.962e+03, 1.703e+03, 1.507e+03, 1.354e+03, 1.231e+03, 1.131e+03, 8.130e+02, 6.434e+02, 5.369e+02, 4.634e+02, 3.678e+02, 3.080e+02, 2.668e+02, 2.367e+02, 2.137e+02, 1.955e+02, 1.807e+02, 1.352e+02, 1.116e+02, 9.726e+01, 8.760e+01, 7.553e+01, 6.839e+01, 6.376e+01, 6.056e+01, 5.827e+01, 5.657e+01, 5.528e+01},
			{3.604e+03, 4.003e+03, 4.712e+03, 5.325e+03, 5.860e+03, 6.329e+03, 6.736e+03, 7.087e+03, 7.387e+03, 8.299e+03, 8.635e+03, 8.719e+03, 8.688e+03, 8.485e+03, 8.216e+03, 7.933e+03, 7.653e+03, 7.384e+03, 7.127e+03, 6.884e+03, 5.857e+03, 5.081e+03, 4.481e+03, 4.004e+03, 3.300e+03, 2.808e+03, 2.446e+03, 2.169e+03, 1.951e+03, 1.775e+03, 1.630e+03, 1.171e+03, 9.263e+02, 7.727e+02, 6.668e+02, 5.292e+02, 4.433e+02, 3.841e+02, 3.409e+02, 3.077e+02, 2.815e+02, 2.603e+02, 1.947e+02, 1.608e+02, 1.401e+02, 1.262e+02, 1.089e+02, 9.858e+01, 9.189e+01, 8.729e+01, 8.399e+01, 8.155e+01, 7.968e+01},
			{3.882e+03, 4.328e+03, 5.131e+03, 5.837e+03, 6.463e+03, 7.021e+03, 7.517e+03, 7.954e+03, 8.338e+03, 9.594e+03, 1.014e+04, 1.036e+04, 1.040e+04, 1.028e+04, 1.004e+04, 9.768e+03, 9.484e+03, 9.203e+03, 8.930e+03, 8.667e+03, 7.517e+03, 6.613e+03, 5.892e+03, 5.308e+03, 4.423e+03, 3.788e+03, 3.314e+03, 2.947e+03, 2.655e+03, 2.418e+03, 2.222e+03, 1.596e+03, 1.261e+03, 1.052e+03, 9.071e+02, 7.200e+02, 6.030e+02, 5.227e+02, 4.638e+02, 4.188e+02, 3.832e+02, 3.543e+02, 2.651e+02, 2.190e+02, 1.909e+02, 1.719e+02, 1.483e+02, 1.343e+02, 1.252e+02, 1.189e+02, 1.144e+02, 1.111e+02, 1.086e+02},
			{4.121e+03, 4.606e+03, 5.494e+03, 6.287e+03, 6.998e+03, 7.639e+03, 8.217e+03, 8.737e+03, 9.201e+03, 1.081e+04, 1.160e+04, 1.196e+04, 1.210e+04, 1.208e+04, 1.189e+04, 1.164e+04, 1.136e+04, 1.108e+04, 1.080e+04, 1.052e+04, 9.279e+03, 8.261e+03, 7.431e+03, 6.743e+03, 5.679e+03, 4.898e+03, 4.304e+03, 3.839e+03, 3.466e+03, 3.161e+03, 2.907e+03, 2.090e+03, 1.650e+03, 1.374e+03, 1.185e+03, 9.400e+02, 7.873e+02, 6.824e+02, 6.057e+02, 5.469e+02, 5.005e+02, 4.628e+02, 3.464e+02, 2.862e+02, 2.494e+02, 2.247e+02, 1.938e+02, 1.756e+02, 1.637e+02, 1.555e+02, 1.496e+02, 1.452e+02, 1.419e+02},
			{4.295e+03, 4.811e+03, 5.768e+03, 6.634e+03, 7.420e+03, 8.134e+03, 8.786e+03, 9.379e+03, 9.915e+03, 1.186e+04, 1.289e+04, 1.341e+04, 1.365e+04, 1.375e+04, 1.362e+04, 1.340e+04, 1.314e+04, 1.286e+04, 1.257e+04, 1.229e+04, 1.098e+04, 9.880e+03, 8.960e+03, 8.187e+03, 6.969e+03, 6.057e+03, 5.353e+03, 4.795e+03, 4.343e+03, 3.970e+03, 3.658e+03, 2.642e+03, 2.088e+03, 1.739e+03, 1.499e+03, 1.189e+03, 9.960e+02, 8.634e+02, 7.663e+02, 6.921e+02, 6.333e+02, 5.857e+02, 4.386e+02, 3.624e+02, 3.159e+02, 2.846e+02, 2.455e+02, 2.224e+02, 2.073e+02, 1.970e+02, 1.895e+02, 1.840e+02, 1.798e+02},
			{4.451e+03, 4.991e+03, 6.004e+03, 6.934e+03, 7.785e+03, 8.566e+03, 9.284e+03, 9.943e+03, 1.055e+04, 1.281e+04, 1.410e+04, 1.479e+04, 1.515e+04, 1.538e+04, 1.532e+04, 1.514e+04, 1.490e+04, 1.463e+04, 1.434e+04, 1.406e+04, 1.270e+04, 1.153e+04, 1.053e+04, 9.682e+03, 8.321e+03, 7.285e+03, 6.473e+03, 5.823e+03, 5.292e+03, 4.850e+03, 4.478e+03, 3.252e+03, 2.574e+03, 2.146e+03, 1.850e+03, 1.467e+03, 1.229e+03, 1.065e+03, 9.457e+02, 8.542e+02, 7.818e+02, 7.230e+02, 5.416e+02, 4.476e+02, 3.902e+02, 3.516e+02, 3.034e+02, 2.748e+02, 2.562e+02, 2.434e+02, 2.342e+02, 2.274e+02, 2.222e+02},
			{4.591e+03, 5.155e+03, 6.224e+03, 7.220e+03, 8.143e+03, 8.998e+03, 9.791e+03, 1.053e+04, 1.121e+04, 1.385e+04, 1.545e+04, 1.638e+04, 1.692e+04, 1.737e+04, 1.744e+04, 1.734e+04, 1.716e+04, 1.693e+04, 1.668e+04, 1.641e+04, 1.504e+04, 1.380e+04, 1.271e+04, 1.175e+04, 1.019e+04, 8.967e+03, 7.996e+03, 7.207e+03, 6.556e+03, 6.011e+03, 5.549e+03, 4.015e+03, 3.163e+03, 2.626e+03, 2.257e+03, 1.783e+03, 1.490e+03, 1.291e+03, 1.145e+03, 1.034e+03, 9.461e+02, 8.750e+02, 6.555e+02, 5.419e+02, 4.724e+02, 4.257e+02, 3.673e+02, 3.327e+02, 3.103e+02, 2.948e+02, 2.836e+02, 2.754e+02, 2.691e+02},
			{4.754e+03, 5.318e+03, 6.396e+03, 7.414e+03, 8.366e+03, 9.255e+03, 1.008e+04, 1.086e+04, 1.158e+04, 1.446e+04, 1.628e+04, 1.737e+04, 1.802e+04, 1.860e+04, 1.873e+04, 1.865e+04, 1.848e+04, 1.825e+04, 1.799e+04, 1.772e+04, 1.631e+04, 1.503e+04, 1.389e+04, 1.291e+04, 1.128e+04, 9.998e+03, 8.973e+03, 8.134e+03, 7.438e+03, 6.851e+03, 6.350e+03, 4.662e+03, 3.705e+03, 3.092e+03, 2.666e+03, 2.113e+03, 1.769e+03, 1.534e+03, 1.361e+03, 1.230e+03, 1.125e+03, 1.041e+03, 7.802e+02, 6.451e+02, 5.625e+02, 5.070e+02, 4.375e+02, 3.963e+02, 3.696e+02, 3.511e+02, 3.379e+02, 3.281e+02, 3.206e+02},
			{4.911e+03, 5.486e+03, 6.585e+03, 7.630e+03, 8.619e+03, 9.551e+03, 1.043e+04, 1.125e+04, 1.202e+04, 1.515e+04, 1.723e+04, 1.853e+04, 1.933e+04, 2.011e+04, 2.035e+04, 2.035e+04, 2.022e+04, 2.002e+04, 1.978e+04, 1.952e+04, 1.812e+04, 1.680e+04, 1.561e+04, 1.456e+04, 1.282e+04, 1.143e+04, 1.030e+04, 9.377e+03, 8.600e+03, 7.941e+03, 7.376e+03, 5.450e+03, 4.342e+03, 3.627e+03, 3.129e+03, 2.481e+03, 2.077e+03, 1.800e+03, 1.597e+03, 1.443e+03, 1.321e+03, 1.222e+03, 9.158e+02, 7.574e+02, 6.605e+02, 5.954e+02, 5.138e+02, 4.655e+02, 4.341e+02, 4.125e+02, 3.969e+02, 3.854e+02, 3.767e+02},
			{5.070e+03, 5.653e+03, 6.770e+03, 7.838e+03, 8.856e+03, 9.823e+03, 1.074e+04, 1.160e+04, 1.241e+04, 1.577e+04, 1.809e+04, 1.959e+04, 2.055e+04, 2.153e+04, 2.190e+04, 2.197e+04, 2.190e+04, 2.173e+04, 2.151e+04, 2.126e+04, 1.988e+04, 1.852e+04, 1.730e+04, 1.620e+04, 1.435e+04, 1.287e+04, 1.165e+04, 1.064e+04, 9.788e+03, 9.061e+03, 8.434e+03, 6.276e+03, 5.017e+03, 4.198e+03, 3.625e+03, 2.876e+03, 2.408e+03, 2.087e+03, 1.852e+03, 1.673e+03, 1.531e+03, 1.417e+03, 1.062e+03, 8.787e+02, 7.664e+02, 6.909e+02, 5.964e+02, 5.404e+02, 5.039e+02, 4.788e+02, 4.608e+02, 4.475e+02, 4.373e+02},
			{5.262e+03, 5.854e+03, 6.987e+03, 8.074e+03, 9.119e+03, 1.012e+04, 1.106e+04, 1.196e+04, 1.282e+04, 1.640e+04, 1.894e+04, 2.066e+04, 2.178e+04, 2.299e+04, 2.350e+04, 2.367e+04, 2.366e+04, 2.354e+04, 2.336e+04, 2.313e+04, 2.179e+04, 2.042e+04, 1.916e+04, 1.802e+04, 1.607e+04, 1.447e+04, 1.316e+04, 1.205e+04, 1.111e+04, 1.031e+04, 9.615e+03, 7.192e+03, 5.761e+03, 4.825e+03, 4.167e+03, 3.305e+03, 2.766e+03, 2.397e+03, 2.127e+03, 1.921e+03, 1.760e+03, 1.626e+03, 1.220e+03, 1.009e+03, 8.803e+02, 7.936e+02, 6.851e+02, 6.208e+02, 5.790e+02, 5.502e+02, 5.295e+02, 5.142e+02, 5.025e+02},
			{5.413e+03, 6.019e+03, 7.176e+03, 8.287e+03, 9.360e+03, 1.039e+04, 1.137e+04, 1.231e+04, 1.320e+04, 1.699e+04, 1.976e+04, 2.168e+04, 2.298e+04, 2.443e+04, 2.509e+04, 2.536e+04, 2.542e+04, 2.535e+04, 2.521e+04, 2.501e+04, 2.373e+04, 2.237e+04, 2.107e+04, 1.989e+04, 1.784e+04, 1.615e+04, 1.473e+04, 1.354e+04, 1.252e+04, 1.164e+04, 1.087e+04, 8.173e+03, 6.562e+03, 5.500e+03, 4.750e+03, 3.767e+03, 3.151e+03, 2.729e+03, 2.421e+03, 2.187e+03, 2.001e+03, 1.851e+03, 1.388e+03, 1.148e+03, 1.002e+03, 9.035e+02, 7.801e+02, 7.070e+02, 6.594e+02, 6.266e+02, 6.030e+02, 5.856e+02, 5.723e+02},
			{5.617e+03, 6.231e+03, 7.398e+03, 8.521e+03, 9.610e+03, 1.066e+04, 1.167e+04, 1.263e+04, 1.355e+04, 1.752e+04, 2.050e+04, 2.262e+04, 2.408e+04, 2.577e+04, 2.658e+04, 2.695e+04, 2.708e+04, 2.707e+04, 2.696e+04, 2.679e+04, 2.558e+04, 2.421e+04, 2.290e+04, 2.168e+04, 1.955e+04, 1.777e+04, 1.627e+04, 1.500e+04, 1.390e+04, 1.295e+04, 1.212e+04, 9.169e+03, 7.385e+03, 6.200e+03, 5.360e+03, 4.253e+03, 3.559e+03, 3.082e+03, 2.734e+03, 2.469e+03, 2.259e+03, 2.090e+03, 1.567e+03, 1.297e+03, 1.132e+03, 1.020e+03, 8.813e+02, 7.987e+02, 7.450e+02, 7.080e+02, 6.814e+02, 6.617e+02, 6.467e+02},
			{5.716e+03, 6.339e+03, 7.520e+03, 8.652e+03, 9.753e+03, 1.082e+04, 1.185e+04, 1.284e+04, 1.378e+04, 1.790e+04, 2.107e+04, 2.338e+04, 2.500e+04, 2.693e+04, 2.789e+04, 2.836e+04, 2.856e+04, 2.859e+04, 2.852e+04, 2.838e+04, 2.723e+04, 2.587e+04, 2.454e+04, 2.330e+04, 2.110e+04, 1.926e+04, 1.769e+04, 1.635e+04, 1.519e+04, 1.418e+04, 1.330e+04, 1.014e+04, 8.202e+03, 6.906e+03, 5.982e+03, 4.757e+03, 3.984e+03, 3.452e+03, 3.064e+03, 2.767e+03, 2.532e+03, 2.342e+03, 1.757e+03, 1.455e+03, 1.269e+03, 1.145e+03, 9.887e+02, 8.962e+02, 8.360e+02, 7.945e+02, 7.646e+02, 7.426e+02, 7.257e+02},
  }
};

/** ----------------------------------------------- ICRU49/73 DATA --------------------------------------------- */

/** For the time being this array holds the data
 * as given by ICRU49 (H+He) and ICRU 73 (> He)
 * for liquid water
 * TODO: Later use this array to read in data
 * TODO: from table
*/

/**
 * @struct AT_stopping_power_ICRU_table
 * TODO
 */
typedef struct {
  const long       material_no;
  const long       number_of_data_points;
  const double     E_MeV_u_and_stopping_power_total_MeV_cm2_g[19][53];
} AT_stopping_power_ICRU_table_struct;

static const AT_stopping_power_ICRU_table_struct AT_stopping_power_ICRU_table = {
  Water_Liquid,
  53,
  {
			{0.025000, 0.030000, 0.040000, 0.050000, 0.060000, 0.070000, 0.080000, 0.090000, 0.100000, 0.150000, 0.200000, 0.250000, 0.300000, 0.400000, 0.500000, 0.600000, 0.700000, 0.800000, 0.900000, 1.000000, 1.500000, 2.000000, 2.500000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000, 10.000000, 15.000000, 20.000000, 25.000000, 30.000000, 40.000000, 50.000000, 60.000000, 70.000000, 80.000000, 90.000000, 100.000000, 150.000000, 200.000000, 250.000000, 300.000000, 400.000000, 500.000000, 600.000000, 700.000000, 800.000000, 900.000000, 1000.000000},
			{6.245e-01, 6.671e-01, 7.324e-01, 7.768e-01, 8.050e-01, 8.205e-01, 8.260e-01, 8.239e-01, 8.161e-01, 7.371e-01, 6.613e-01, 6.006e-01, 5.504e-01, 4.719e-01, 4.132e-01, 3.680e-01, 3.325e-01, 3.039e-01, 2.805e-01, 2.608e-01, 1.957e-01, 1.586e-01, 1.344e-01, 1.172e-01, 9.404e-02, 7.911e-02, 6.858e-02, 6.071e-02, 5.460e-02, 4.969e-02, 4.567e-02, 3.292e-02, 2.607e-02, 2.175e-02, 1.876e-02, 1.488e-02, 1.245e-02, 1.078e-02, 9.559e-03, 8.625e-03, 7.888e-03, 7.289e-03, 5.445e-03, 4.492e-03, 3.911e-03, 3.520e-03, 3.032e-03, 2.743e-03, 2.556e-03, 2.426e-03, 2.333e-03, 2.264e-03, 2.211e-03},
			{1.151e+00, 1.255e+00, 1.438e+00, 1.593e+00, 1.722e+00, 1.832e+00, 1.923e+00, 2.002e+00, 2.069e+00, 2.245e+00, 2.260e+00, 2.193e+00, 2.080e+00, 1.840e+00, 1.625e+00, 1.454e+00, 1.315e+00, 1.207e+00, 1.113e+00, 1.035e+00, 7.777e-01, 6.306e-01, 5.344e-01, 4.682e-01, 3.754e-01, 3.146e-01, 2.732e-01, 2.416e-01, 2.180e-01, 1.980e-01, 1.816e-01, 1.309e-01, 1.037e-01, 8.649e-02, 7.505e-02, 5.942e-02, 4.952e-02, 4.297e-02, 3.807e-02, 3.446e-02, 3.145e-02, 2.901e-02, 2.168e-02, 1.789e-02, 1.558e-02,   0,   0,   0,   0,   0,   0,   0,   0},
			{2.626e+00, 2.840e+00, 3.191e+00, 3.461e+00, 3.665e+00, 3.817e+00, 3.927e+00, 4.004e+00, 4.056e+00, 4.102e+00, 3.998e+00, 3.853e+00, 3.702e+00, 3.413e+00, 3.158e+00, 2.934e+00, 2.739e+00, 2.567e+00, 2.415e+00, 2.280e+00, 1.782e+00, 1.465e+00, 1.247e+00, 1.087e+00, 8.706e-01, 7.299e-01, 6.310e-01, 5.575e-01, 5.005e-01, 4.550e-01, 4.178e-01, 3.004e-01, 2.376e-01, 1.981e-01, 1.708e-01, 1.354e-01, 1.132e-01, 9.803e-02, 8.692e-02, 7.842e-02, 7.171e-02, 6.627e-02, 4.950e-02, 4.085e-02, 3.557e-02, 3.203e-02, 2.760e-02, 2.498e-02, 2.327e-02, 2.210e-02, 2.126e-02, 2.064e-02, 2.016e-02},
			{3.272e+00, 3.565e+00, 4.061e+00, 4.463e+00, 4.790e+00, 5.052e+00, 5.258e+00, 5.419e+00, 5.542e+00, 5.803e+00, 5.787e+00, 5.675e+00, 5.529e+00, 5.215e+00, 4.912e+00, 4.634e+00, 4.381e+00, 4.152e+00, 3.944e+00, 3.756e+00, 3.026e+00, 2.533e+00, 2.179e+00, 1.913e+00, 1.542e+00, 1.296e+00, 1.122e+00, 9.911e-01, 8.898e-01, 8.087e-01, 7.423e-01, 5.335e-01, 4.219e-01, 3.518e-01, 3.034e-01, 2.406e-01, 2.013e-01, 1.743e-01, 1.545e-01, 1.394e-01, 1.275e-01, 1.178e-01, 8.805e-02, 7.266e-02, 6.328e-02, 5.698e-02, 4.910e-02, 4.444e-02, 4.141e-02, 3.933e-02, 3.783e-02, 3.672e-02, 3.588e-02},
			{3.773e+00, 4.142e+00, 4.776e+00, 5.304e+00, 5.749e+00, 6.122e+00, 6.431e+00, 6.684e+00, 6.890e+00, 7.432e+00, 7.551e+00, 7.505e+00, 7.391e+00, 7.091e+00, 6.772e+00, 6.463e+00, 6.172e+00, 5.901e+00, 5.650e+00, 5.418e+00, 4.484e+00, 3.817e+00, 3.322e+00, 2.940e+00, 2.392e+00, 2.020e+00, 1.752e+00, 1.549e+00, 1.391e+00, 1.265e+00, 1.161e+00, 8.332e-01, 6.587e-01, 5.492e-01, 4.737e-01, 3.757e-01, 3.144e-01, 2.723e-01, 2.415e-01, 2.179e-01, 1.993e-01, 1.842e-01, 1.376e-01, 1.136e-01, 9.894e-02, 8.909e-02, 7.678e-02, 6.950e-02, 6.477e-02, 6.151e-02, 5.916e-02, 5.743e-02, 5.611e-02},
			{4.154e+00, 4.593e+00, 5.358e+00, 6.009e+00, 6.568e+00, 7.049e+00, 7.460e+00, 7.809e+00, 8.103e+00, 8.968e+00, 9.262e+00, 9.311e+00, 9.250e+00, 8.994e+00, 8.680e+00, 8.358e+00, 8.045e+00, 7.747e+00, 7.465e+00, 7.199e+00, 6.093e+00, 5.269e+00, 4.636e+00, 4.137e+00, 3.403e+00, 2.891e+00, 2.516e+00, 2.230e+00, 2.004e+00, 1.823e+00, 1.673e+00, 1.200e+00, 9.483e-01, 7.904e-01, 6.817e-01, 5.406e-01, 4.525e-01, 3.920e-01, 3.477e-01, 3.138e-01, 2.870e-01, 2.653e-01, 1.983e-01, 1.637e-01, 1.426e-01, 1.284e-01, 1.107e-01, 1.002e-01, 9.335e-02, 8.865e-02, 8.528e-02, 8.278e-02, 8.088e-02},
			{4.490e+00, 4.984e+00, 5.860e+00, 6.616e+00, 7.276e+00, 7.854e+00, 8.360e+00, 8.799e+00, 9.179e+00, 1.039e+01, 1.089e+01, 1.107e+01, 1.108e+01, 1.090e+01, 1.061e+01, 1.030e+01, 9.974e+00, 9.660e+00, 9.357e+00, 9.068e+00, 7.823e+00, 6.859e+00, 6.097e+00, 5.484e+00, 4.560e+00, 3.900e+00, 3.408e+00, 3.029e+00, 2.727e+00, 2.482e+00, 2.280e+00, 1.636e+00, 1.291e+00, 1.076e+00, 9.274e-01, 7.354e-01, 6.156e-01, 5.333e-01, 4.731e-01, 4.270e-01, 3.906e-01, 3.611e-01, 2.700e-01, 2.229e-01, 1.942e-01, 1.749e-01, 1.507e-01, 1.365e-01, 1.272e-01, 1.208e-01, 1.162e-01, 1.128e-01, 1.102e-01},
			{4.778e+00, 5.321e+00, 6.298e+00, 7.152e+00, 7.907e+00, 8.578e+00, 9.173e+00, 9.700e+00, 1.016e+01, 1.173e+01, 1.246e+01, 1.278e+01, 1.289e+01, 1.281e+01, 1.257e+01, 1.227e+01, 1.195e+01, 1.163e+01, 1.132e+01, 1.101e+01, 9.659e+00, 8.571e+00, 7.691e+00, 6.967e+00, 5.854e+00, 5.042e+00, 4.427e+00, 3.945e+00, 3.560e+00, 3.245e+00, 2.983e+00, 2.142e+00, 1.689e+00, 1.406e+00, 1.212e+00, 9.602e-01, 8.037e-01, 6.963e-01, 6.178e-01, 5.577e-01, 5.102e-01, 4.716e-01, 3.527e-01, 2.913e-01, 2.538e-01, 2.285e-01, 1.970e-01, 1.784e-01, 1.663e-01, 1.579e-01, 1.519e-01, 1.475e-01, 1.441e-01},
			{4.992e+00, 5.575e+00, 6.637e+00, 7.578e+00, 8.418e+00, 9.171e+00, 9.847e+00, 1.045e+01, 1.100e+01, 1.290e+01, 1.388e+01, 1.435e+01, 1.456e+01, 1.459e+01, 1.440e+01, 1.413e+01, 1.383e+01, 1.351e+01, 1.319e+01, 1.287e+01, 1.144e+01, 1.026e+01, 9.279e+00, 8.463e+00, 7.187e+00, 6.237e+00, 5.506e+00, 4.928e+00, 4.461e+00, 4.076e+00, 3.753e+00, 2.707e+00, 2.137e+00, 1.779e+00, 1.533e+00, 1.215e+00, 1.017e+00, 8.809e-01, 7.816e-01, 7.056e-01, 6.456e-01, 5.969e-01, 4.466e-01, 3.688e-01, 3.213e-01, 2.894e-01, 2.496e-01, 2.259e-01, 2.106e-01, 2.000e-01, 1.924e-01, 1.868e-01, 1.825e-01},
			{5.182e+00, 5.797e+00, 6.931e+00, 7.948e+00, 8.865e+00, 9.693e+00, 1.044e+01, 1.112e+01, 1.174e+01, 1.398e+01, 1.521e+01, 1.585e+01, 1.617e+01, 1.633e+01, 1.621e+01, 1.598e+01, 1.569e+01, 1.538e+01, 1.506e+01, 1.474e+01, 1.324e+01, 1.198e+01, 1.091e+01, 1.001e+01, 8.584e+00, 7.503e+00, 6.660e+00, 5.986e+00, 5.436e+00, 4.979e+00, 4.595e+00, 3.332e+00, 2.635e+00, 2.195e+00, 1.892e+00, 1.499e+00, 1.255e+00, 1.087e+00, 9.646e-01, 8.709e-01, 7.969e-01, 7.368e-01, 5.514e-01, 4.555e-01, 3.969e-01, 3.576e-01, 3.083e-01, 2.792e-01, 2.602e-01, 2.472e-01, 2.378e-01, 2.308e-01, 2.255e-01},
			{5.352e+00, 5.998e+00, 7.203e+00, 8.300e+00, 9.298e+00, 1.021e+01, 1.104e+01, 1.181e+01, 1.250e+01, 1.513e+01, 1.668e+01, 1.756e+01, 1.804e+01, 1.843e+01, 1.844e+01, 1.829e+01, 1.805e+01, 1.778e+01, 1.748e+01, 1.718e+01, 1.567e+01, 1.432e+01, 1.315e+01, 1.214e+01, 1.050e+01, 9.226e+00, 8.218e+00, 7.401e+00, 6.728e+00, 6.166e+00, 5.690e+00, 4.112e+00, 3.237e+00, 2.686e+00, 2.307e+00, 1.821e+00, 1.521e+00, 1.317e+00, 1.168e+00, 1.054e+00, 9.644e-01, 8.917e-01, 6.674e-01, 5.514e-01, 4.806e-01, 4.329e-01, 3.734e-01, 3.381e-01, 3.152e-01, 2.993e-01, 2.880e-01, 2.796e-01, 2.732e-01},
			{5.542e+00, 6.193e+00, 7.420e+00, 8.551e+00, 9.590e+00, 1.054e+01, 1.142e+01, 1.223e+01, 1.298e+01, 1.585e+01, 1.762e+01, 1.866e+01, 1.926e+01, 1.976e+01, 1.983e+01, 1.970e+01, 1.947e+01, 1.920e+01, 1.889e+01, 1.858e+01, 1.702e+01, 1.562e+01, 1.441e+01, 1.336e+01, 1.164e+01, 1.030e+01, 9.233e+00, 8.362e+00, 7.640e+00, 7.033e+00, 6.516e+00, 4.777e+00, 3.792e+00, 3.162e+00, 2.725e+00, 2.159e+00, 1.806e+00, 1.565e+00, 1.388e+00, 1.254e+00, 1.147e+00, 1.061e+00, 7.944e-01, 6.565e-01, 5.722e-01, 5.156e-01, 4.447e-01, 4.027e-01, 3.754e-01, 3.566e-01, 3.430e-01, 3.330e-01, 3.254e-01},
			{5.724e+00, 6.390e+00, 7.649e+00, 8.820e+00, 9.905e+00, 1.091e+01, 1.184e+01, 1.271e+01, 1.351e+01, 1.666e+01, 1.869e+01, 1.993e+01, 2.068e+01, 2.138e+01, 2.156e+01, 2.150e+01, 2.132e+01, 2.107e+01, 2.079e+01, 2.048e+01, 1.891e+01, 1.747e+01, 1.619e+01, 1.508e+01, 1.324e+01, 1.178e+01, 1.061e+01, 9.641e+00, 8.835e+00, 8.153e+00, 7.569e+00, 5.583e+00, 4.444e+00, 3.710e+00, 3.199e+00, 2.534e+00, 2.120e+00, 1.836e+00, 1.629e+00, 1.471e+00, 1.346e+00, 1.245e+00, 9.325e-01, 7.707e-01, 6.719e-01, 6.054e-01, 5.223e-01, 4.730e-01, 4.410e-01, 4.189e-01, 4.030e-01, 3.912e-01, 3.823e-01},
			{5.905e+00, 6.583e+00, 7.868e+00, 9.073e+00, 1.020e+01, 1.125e+01, 1.223e+01, 1.314e+01, 1.399e+01, 1.740e+01, 1.966e+01, 2.110e+01, 2.201e+01, 2.291e+01, 2.321e+01, 2.322e+01, 2.309e+01, 2.287e+01, 2.261e+01, 2.232e+01, 2.076e+01, 1.928e+01, 1.795e+01, 1.678e+01, 1.483e+01, 1.326e+01, 1.199e+01, 1.094e+01, 1.006e+01, 9.304e+00, 8.656e+00, 6.430e+00, 5.135e+00, 4.294e+00, 3.705e+00, 2.938e+00, 2.458e+00, 2.129e+00, 1.889e+00, 1.706e+00, 1.561e+00, 1.444e+00, 1.082e+00, 8.942e-01, 7.796e-01, 7.026e-01, 6.061e-01, 5.490e-01, 5.119e-01, 4.862e-01, 4.678e-01, 4.542e-01, 4.438e-01},
			{6.120e+00, 6.810e+00, 8.118e+00, 9.352e+00, 1.051e+01, 1.161e+01, 1.263e+01, 1.358e+01, 1.448e+01, 1.813e+01, 2.063e+01, 2.228e+01, 2.334e+01, 2.447e+01, 2.491e+01, 2.502e+01, 2.495e+01, 2.478e+01, 2.455e+01, 2.428e+01, 2.276e+01, 2.126e+01, 1.989e+01, 1.867e+01, 1.659e+01, 1.492e+01, 1.354e+01, 1.239e+01, 1.142e+01, 1.059e+01, 9.867e+00, 7.367e+00, 5.896e+00, 4.935e+00, 4.259e+00, 3.376e+00, 2.824e+00, 2.445e+00, 2.169e+00, 1.959e+00, 1.792e+00, 1.657e+00, 1.242e+00, 1.027e+00, 8.954e-01, 8.070e-01, 6.963e-01, 6.308e-01, 5.881e-01, 5.587e-01, 5.375e-01, 5.219e-01, 5.100e-01},
			{6.294e+00, 7.000e+00, 8.338e+00, 9.604e+00, 1.080e+01, 1.194e+01, 1.300e+01, 1.400e+01, 1.494e+01, 1.882e+01, 2.155e+01, 2.341e+01, 2.465e+01, 2.601e+01, 2.660e+01, 2.681e+01, 2.681e+01, 2.669e+01, 2.650e+01, 2.626e+01, 2.479e+01, 2.328e+01, 2.188e+01, 2.061e+01, 1.843e+01, 1.664e+01, 1.516e+01, 1.392e+01, 1.286e+01, 1.195e+01, 1.115e+01, 8.371e+00, 6.715e+00, 5.624e+00, 4.856e+00, 3.847e+00, 3.217e+00, 2.785e+00, 2.470e+00, 2.229e+00, 2.040e+00, 1.886e+00, 1.413e+00, 1.169e+00, 1.019e+00, 9.187e-01, 7.929e-01, 7.183e-01, 6.697e-01, 6.362e-01, 6.122e-01, 5.944e-01, 5.808e-01},
			{6.522e+00, 7.237e+00, 8.590e+00, 9.875e+00, 1.110e+01, 1.226e+01, 1.336e+01, 1.439e+01, 1.537e+01, 1.945e+01, 2.240e+01, 2.445e+01, 2.586e+01, 2.746e+01, 2.819e+01, 2.850e+01, 2.857e+01, 2.850e+01, 2.834e+01, 2.813e+01, 2.672e+01, 2.521e+01, 2.378e+01, 2.247e+01, 2.020e+01, 1.832e+01, 1.675e+01, 1.542e+01, 1.428e+01, 1.330e+01, 1.244e+01, 9.392e+00, 7.557e+00, 6.340e+00, 5.479e+00, 4.344e+00, 3.633e+00, 3.145e+00, 2.789e+00, 2.517e+00, 2.303e+00, 2.130e+00, 1.596e+00, 1.320e+00, 1.151e+00, 1.038e+00, 8.957e-01, 8.115e-01, 7.567e-01, 7.189e-01, 6.917e-01, 6.717e-01, 6.563e-01},
			{6.642e+00, 7.369e+00, 8.739e+00, 1.004e+01, 1.128e+01, 1.247e+01, 1.359e+01, 1.466e+01, 1.566e+01, 1.993e+01, 2.308e+01, 2.532e+01, 2.689e+01, 2.872e+01, 2.960e+01, 3.001e+01, 3.015e+01, 3.013e+01, 3.000e+01, 2.982e+01, 2.846e+01, 2.694e+01, 2.549e+01, 2.415e+01, 2.181e+01, 1.986e+01, 1.822e+01, 1.682e+01, 1.561e+01, 1.457e+01, 1.365e+01, 1.039e+01, 8.394e+00, 7.063e+00, 6.114e+00, 4.859e+00, 4.067e+00, 3.522e+00, 3.125e+00, 2.821e+00, 2.581e+00, 2.387e+00, 1.789e+00, 1.480e+00, 1.291e+00, 1.164e+00, 1.005e+00, 9.105e-01, 8.491e-01, 8.067e-01, 7.763e-01, 7.537e-01, 7.365e-01},
  }
};

/** ----------------------------------------------- Static objects  --------------------------------------- */


static const AT_stopping_power_sources_struct AT_stopping_power_sources = {
	STOPPING_POWER_SOURCE_N,
    {  PSTAR,         Bethe, 			ShieldHit, 				ICRU           }, // source_no
    {  "PSTAR data",  "Bethe formula",  "ShieldHit (Bethe)",    "ICRU 49&73" }    // source_name
};


static const AT_stopping_power_tabulated_source_group_for_all_materials_struct AT_data_PSTAR_source = {
  MATERIAL_DATA_N,
  PSTAR,
  { User_Defined_Material,  Water_Liquid,                         Aluminum_Oxide,                                 Aluminum,                                 PMMA,                                Alanine,                                LiF,                                Air,								Silicon,								Copper },
  { NULL,                   &AT_stopping_power_data_PSTAR_Water,  &AT_stopping_power_data_PSTAR_Aluminum_Oxide,   &AT_stopping_power_data_PSTAR_Aluminum,   &AT_stopping_power_data_PSTAR_PMMA,  &AT_stopping_power_data_PSTAR_Alanine,  &AT_stopping_power_data_PSTAR_LiF,  &AT_stopping_power_data_PSTAR_Air,	&AT_stopping_power_data_PSTAR_Silicon,	&AT_stopping_power_data_PSTAR_Copper}
};


static const AT_stopping_power_tabulated_source_group_struct AT_stopping_power_tabulated_source = {
  STOPPING_POWER_SOURCE_N,
  { PSTAR,                  Bethe, 	ShieldHit,	ICRU},
  { &AT_data_PSTAR_source,  NULL,	NULL,       NULL}
};


static const AT_stopping_power_analytical_sources_struct AT_stopping_power_analytical_source = {
  STOPPING_POWER_SOURCE_N,
  { PSTAR,   Bethe,				ShieldHit,             ICRU},
  { NULL,    &AT_Bethe_wrapper, &AT_ShieldHit_wrapper, &AT_ICRU_wrapper}
};

#endif /* AT_DATASTOPPINGPOWER_H_ */
