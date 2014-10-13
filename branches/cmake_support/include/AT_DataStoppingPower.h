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
#include "AT_EnergyLoss.h"

/**
 * @enum stoppingPowerSource_no Stopping Power source code numbers
 */
enum stoppingPowerSource_no{
  PSTAR                = 0, /**< PSTAR */
  Bethe                = 1, /**< Bethe */
  ShieldHit			   = 2, /**< ShieldHit code: extended Bethe formula */
  ICRU                 = 3, /**< ICRU 49 and 73: for liquid water */
  FLUKA                = 4, /**< DEDX FLUKA WATER I=77.3eV */
  ATIMA                = 5  /**< Data from ATIMA for TLD sim Daria Bascolo (GSI) */
};


/**
 * TODO
 */
#define STOPPING_POWER_SOURCE_N    6


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
 * Wrapper for TRiP98 stopping powers
 * momentarily read from hardwired struct for water
 * TODO: augment by routine to read data file into
 * TODO: struct and look-up value
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return stopping power (MeV cm2 per g)
 */
double AT_FLUKA_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no);

/**
 * Wrapper for ATIMA stopping powers
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return stopping power (MeV cm2 per g)
 */
double AT_ATIMA_wrapper(const double 	E_MeV_u,
		const long 	    particle_no,
		const long 		material_no);

/**
 * Wrapper for stopping powers from Andrea Mairani's Radical Diffusion Model
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @return stopping power (MeV cm2 per g)
 */
double AT_AMRadDiff_wrapper(const double 	E_MeV_u,
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

/** ----------------------------------------------- FLUKA DATA --------------------------------------------- */

/**
 * @struct AT_stopping_power_FLUKA_table
 * TODO
 */
#define NDATA_FLUKA	40
typedef struct {
  const long       material_no;
  const long       number_of_data_points;
  const double     E_MeV_u_and_stopping_power_total_MeV_cm2_g[1+8][NDATA_FLUKA];
} AT_stopping_power_FLUKA_table_struct;

/** PROPRIETARY DATA --- NOT TO BE COMMITTED INTo PUBLIC SVN
	BEGIN */
static const AT_stopping_power_FLUKA_table_struct AT_stopping_power_FLUKA_table = {
  Water_Liquid,
  NDATA_FLUKA,
  {
		  // E / (MeV/u)
		  {1.260000e-01, 1.580000e-01, 2.000000e-01, 2.510000e-01, 3.160000e-01, 3.980000e-01, 5.010000e-01, 6.310000e-01, 7.940000e-01, 1.000000e+00,
		  1.260000e+00, 1.580000e+00, 2.000000e+00, 2.510000e+00, 3.160000e+00, 3.980000e+00, 5.010000e+00, 6.310000e+00, 7.940000e+00, 1.000000e+01,
		  1.260000e+01, 1.580000e+01, 2.000000e+01, 2.510000e+01, 3.160000e+01, 3.980000e+01, 5.010000e+01, 6.310000e+01, 7.940000e+01, 1.000000e+02,
		  1.260000e+02, 1.580000e+02, 2.000000e+02, 2.510000e+02, 3.160000e+02, 3.980000e+02, 5.010000e+02, 6.310000e+02, 7.940000e+02, 1.000000e+03},
		  // Z = 1
		  {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		  },
		  // Z = 0
		  {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		  },
		  // Z = 0
		  {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		  },
		  // Z = 0
		  {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		  },
		  // Z = 0
		  {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		  },
		  // Z = 0
		  {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		  },
		  // Z = 0
		  {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		  },
		  // Z = 0
		  {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		  0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		  }
  }
};
/** PROPRIETARY DATA --- NOT TO BE COMMITTED INTo PUBLIC SVN
	END */

/** ----------------------------------------------- ATIMA DATA --------------------------------------------- */

/**
 * @struct AT_stopping_power_ATIMA_table
 */
typedef struct {
  const long       material_no;
  const long       number_of_data_points;
  const double     E_MeV_u_and_stopping_power_total_MeV_cm2_g[93][91];
} AT_stopping_power_ATIMA_table_struct;

/** PROPRIETARY DATA --- NOT TO BE COMMITTED INTo PUBLIC SVN
	BEGIN */
static const AT_stopping_power_ATIMA_table_struct AT_stopping_power_ATIMA_table = {
  LiF,
  91,
{{0.010000, 0.011659, 0.013594, 0.015849, 0.018478, 0.021544, 0.025119, 0.029286, 0.034145, 0.039811, 0.046416, 0.054117, 0.063096, 0.073564, 0.085770, 0.100000, 0.116591, 0.135936, 0.158489, 0.184785, 0.215443, 0.251189, 0.292864, 0.341455, 0.398107, 0.464159, 0.541170, 0.630957, 0.735642, 0.857696, 1.000000, 1.165914, 1.359356, 1.584893, 1.847850, 2.154435, 2.511886, 2.928645, 3.414549, 3.981072, 4.641589, 5.411695, 6.309573, 7.356423, 8.576959, 10.000000, 11.659144, 13.593564, 15.848932, 18.478498, 21.544347, 25.118864, 29.286446, 34.145489, 39.810717, 46.415888, 54.116953, 63.095734, 73.564225, 85.769590, 100.000000, 116.591440, 135.935639, 158.489319, 184.784980, 215.443469, 251.188643, 292.864456, 341.454887, 398.107171, 464.158883, 541.169527, 630.957344, 735.642254, 857.695899, 1000.000000, 1165.914401, 1359.356391, 1584.893192, 1847.849797, 2154.434690, 2511.886432, 2928.644565, 3414.548874, 3981.071706, 4641.588834, 5411.695265, 6309.573445, 7356.422545, 8576.958986, 10000.000000},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00},
{0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00}}  
};
/** PROPRIETARY DATA --- NOT TO BE COMMITTED INTO PUBLIC SVN
	END */


/**
 * Defines the total number of data entries in PSTAR structure below
 */
#define PSTAR_DATA_N_PER_MATERIAL       132
#define PSTAR_DATA_N                    PSTAR_DATA_N_PER_MATERIAL * MATERIAL_DATA_N


/**
 * Structure definition for PSTAR stopping power data
 */
typedef struct {
  const long     n;                                     /** number of items in local PSTAR data table */
  const double   kin_E_MeV[PSTAR_DATA_N];               /** vector of kinetic energies (of length 132) */
  const double   stp_pow_el_MeV_cm2_g[PSTAR_DATA_N];    /** Electronic (Collision) Stopping Power */
  const double   stp_pow_nuc_MeV_cm2_g[PSTAR_DATA_N];   /** Nuclear Stopping Power */
  const double   stp_pow_tot_MeV_cm2_g[PSTAR_DATA_N];   /** Total Stopping Power */
  const double   range_cdsa_g_cm2[PSTAR_DATA_N];        /** CSDA (continuous-slowing-down approximation) range */
  const double   range_proj_g_cm2[PSTAR_DATA_N];        /** Projected range */
  const double   detour_factor[PSTAR_DATA_N];           /** Detour factor (ratio of projected range to CSDA range) */
  const long     material_no[PSTAR_DATA_N];             /** Material number (see AT_DataMaterial.h) */
} AT_pstar_data_struct;

/**
 * Structure to hold PSTAR stopping power data
 * Tabulated data for protons with energies between 1keV and 10GeV (132 kinetic energies),
 * based on ICRU Report 49
 * See http://physics.nist.gov/PhysRefData/Star/Text/intro.html
 * TODO Should be moved to external file
*/
static const AT_pstar_data_struct AT_PSTAR_Data = {
    PSTAR_DATA_N,
    {
        1.0000e-3f,    1.5000e-3f,    2.0000e-3f,    2.5000e-3f,    3.0000e-3f,    4.0000e-3f,    5.0000e-3f,    6.0000e-3f,    7.0000e-3f,    8.0000e-3f,
        9.0000e-3f,    1.0000e-2f,    1.2500e-2f,    1.5000e-2f,    1.7500e-2f,    2.0000e-2f,    2.2500e-2f,    2.5000e-2f,    2.7500e-2f,    3.0000e-2f,
        3.5000e-2f,    4.0000e-2f,    4.5000e-2f,    5.0000e-2f,    5.5000e-2f,    6.0000e-2f,    6.5000e-2f,    7.0000e-2f,    7.5000e-2f,    8.0000e-2f,
        8.5000e-2f,    9.0000e-2f,    9.5000e-2f,    1.0000e-1f,    1.2500e-1f,    1.5000e-1f,    1.7500e-1f,    2.0000e-1f,    2.2500e-1f,    2.5000e-1f,
        2.7500e-1f,    3.0000e-1f,    3.5000e-1f,    4.0000e-1f,    4.5000e-1f,    5.0000e-1f,    5.5000e-1f,    6.0000e-1f,    6.5000e-1f,    7.0000e-1f,
        7.5000e-1f,    8.0000e-1f,    8.5000e-1f,    9.0000e-1f,    9.5000e-1f,    1.0000e0f,    1.2500e0f,    1.5000e0f,    1.7500e0f,    2.0000e0f,
        2.2500e0f,    2.5000e0f,    2.7500e0f,    3.0000e0f,    3.5000e0f,    4.0000e0f,    4.5000e0f,    5.0000e0f,    5.5000e0f,    6.0000e0f,
        6.5000e0f,    7.0000e0f,    7.5000e0f,    8.0000e0f,    8.5000e0f,    9.0000e0f,    9.5000e0f,    1.0000e1f,    1.2500e1f,    1.5000e1f,
        1.7500e1f,    2.0000e1f,    2.5000e1f,    2.7500e1f,    3.0000e1f,    3.5000e1f,    4.0000e1f,    4.5000e1f,    5.0000e1f,    5.5000e1f,
        6.0000e1f,    6.5000e1f,    7.0000e1f,    7.5000e1f,    8.0000e1f,    8.5000e1f,    9.0000e1f,    9.5000e1f,    1.0000e2f,    1.2500e2f,
        1.5000e2f,    1.7500e2f,    2.0000e2f,    2.2500e2f,    2.5000e2f,    2.7500e2f,    3.0000e2f,    3.5000e2f,    4.0000e2f,    4.5000e2f,
        5.0000e2f,    5.5000e2f,    6.0000e2f,    6.5000e2f,    7.0000e2f,    7.5000e2f,    8.0000e2f,    8.5000e2f,    9.0000e2f,    9.5000e2f,
        1.0000e3f,    1.5000e3f,    2.0000e3f,    2.5000e3f,    3.0000e3f,    4.0000e3f,    5.0000e3f,    6.0000e3f,    7.0000e3f,    8.0000e3f,
        9.0000e3f,    1.0000e4f,    1.0000e-3f,    1.5000e-3f,    2.0000e-3f,    2.5000e-3f,    3.0000e-3f,    4.0000e-3f,    5.0000e-3f,    6.0000e-3f,
        7.0000e-3f,    8.0000e-3f,    9.0000e-3f,    1.0000e-2f,    1.2500e-2f,    1.5000e-2f,    1.7500e-2f,    2.0000e-2f,    2.2500e-2f,    2.5000e-2f,
        2.7500e-2f,    3.0000e-2f,    3.5000e-2f,    4.0000e-2f,    4.5000e-2f,    5.0000e-2f,    5.5000e-2f,    6.0000e-2f,    6.5000e-2f,    7.0000e-2f,
        7.5000e-2f,    8.0000e-2f,    8.5000e-2f,    9.0000e-2f,    9.5000e-2f,    1.0000e-1f,    1.2500e-1f,    1.5000e-1f,    1.7500e-1f,    2.0000e-1f,
        2.2500e-1f,    2.5000e-1f,    2.7500e-1f,    3.0000e-1f,    3.5000e-1f,    4.0000e-1f,    4.5000e-1f,    5.0000e-1f,    5.5000e-1f,    6.0000e-1f,
        6.5000e-1f,    7.0000e-1f,    7.5000e-1f,    8.0000e-1f,    8.5000e-1f,    9.0000e-1f,    9.5000e-1f,    1.0000e0f,    1.2500e0f,    1.5000e0f,
        1.7500e0f,    2.0000e0f,    2.2500e0f,    2.5000e0f,    2.7500e0f,    3.0000e0f,    3.5000e0f,    4.0000e0f,    4.5000e0f,    5.0000e0f,
        5.5000e0f,    6.0000e0f,    6.5000e0f,    7.0000e0f,    7.5000e0f,    8.0000e0f,    8.5000e0f,    9.0000e0f,    9.5000e0f,    1.0000e1f,
        1.2500e1f,    1.5000e1f,    1.7500e1f,    2.0000e1f,    2.5000e1f,    2.7500e1f,    3.0000e1f,    3.5000e1f,    4.0000e1f,    4.5000e1f,
        5.0000e1f,    5.5000e1f,    6.0000e1f,    6.5000e1f,    7.0000e1f,    7.5000e1f,    8.0000e1f,    8.5000e1f,    9.0000e1f,    9.5000e1f,
        1.0000e2f,    1.2500e2f,    1.5000e2f,    1.7500e2f,    2.0000e2f,    2.2500e2f,    2.5000e2f,    2.7500e2f,    3.0000e2f,    3.5000e2f,
        4.0000e2f,    4.5000e2f,    5.0000e2f,    5.5000e2f,    6.0000e2f,    6.5000e2f,    7.0000e2f,    7.5000e2f,    8.0000e2f,    8.5000e2f,
        9.0000e2f,    9.5000e2f,    1.0000e3f,    1.5000e3f,    2.0000e3f,    2.5000e3f,    3.0000e3f,    4.0000e3f,    5.0000e3f,    6.0000e3f,
        7.0000e3f,    8.0000e3f,    9.0000e3f,    1.0000e4f,    1.0000e-3f,    1.5000e-3f,    2.0000e-3f,    2.5000e-3f,    3.0000e-3f,    4.0000e-3f,
        5.0000e-3f,    6.0000e-3f,    7.0000e-3f,    8.0000e-3f,    9.0000e-3f,    1.0000e-2f,    1.2500e-2f,    1.5000e-2f,    1.7500e-2f,    2.0000e-2f,
        2.2500e-2f,    2.5000e-2f,    2.7500e-2f,    3.0000e-2f,    3.5000e-2f,    4.0000e-2f,    4.5000e-2f,    5.0000e-2f,    5.5000e-2f,    6.0000e-2f,
        6.5000e-2f,    7.0000e-2f,    7.5000e-2f,    8.0000e-2f,    8.5000e-2f,    9.0000e-2f,    9.5000e-2f,    1.0000e-1f,    1.2500e-1f,    1.5000e-1f,
        1.7500e-1f,    2.0000e-1f,    2.2500e-1f,    2.5000e-1f,    2.7500e-1f,    3.0000e-1f,    3.5000e-1f,    4.0000e-1f,    4.5000e-1f,    5.0000e-1f,
        5.5000e-1f,    6.0000e-1f,    6.5000e-1f,    7.0000e-1f,    7.5000e-1f,    8.0000e-1f,    8.5000e-1f,    9.0000e-1f,    9.5000e-1f,    1.0000e0f,
        1.2500e0f,    1.5000e0f,    1.7500e0f,    2.0000e0f,    2.2500e0f,    2.5000e0f,    2.7500e0f,    3.0000e0f,    3.5000e0f,    4.0000e0f,
        4.5000e0f,    5.0000e0f,    5.5000e0f,    6.0000e0f,    6.5000e0f,    7.0000e0f,    7.5000e0f,    8.0000e0f,    8.5000e0f,    9.0000e0f,
        9.5000e0f,    1.0000e1f,    1.2500e1f,    1.5000e1f,    1.7500e1f,    2.0000e1f,    2.5000e1f,    2.7500e1f,    3.0000e1f,    3.5000e1f,
        4.0000e1f,    4.5000e1f,    5.0000e1f,    5.5000e1f,    6.0000e1f,    6.5000e1f,    7.0000e1f,    7.5000e1f,    8.0000e1f,    8.5000e1f,
        9.0000e1f,    9.5000e1f,    1.0000e2f,    1.2500e2f,    1.5000e2f,    1.7500e2f,    2.0000e2f,    2.2500e2f,    2.5000e2f,    2.7500e2f,
        3.0000e2f,    3.5000e2f,    4.0000e2f,    4.5000e2f,    5.0000e2f,    5.5000e2f,    6.0000e2f,    6.5000e2f,    7.0000e2f,    7.5000e2f,
        8.0000e2f,    8.5000e2f,    9.0000e2f,    9.5000e2f,    1.0000e3f,    1.5000e3f,    2.0000e3f,    2.5000e3f,    3.0000e3f,    4.0000e3f,
        5.0000e3f,    6.0000e3f,    7.0000e3f,    8.0000e3f,    9.0000e3f,    1.0000e4f,    1.0000e-3f,    1.5000e-3f,    2.0000e-3f,    2.5000e-3f,
        3.0000e-3f,    4.0000e-3f,    5.0000e-3f,    6.0000e-3f,    7.0000e-3f,    8.0000e-3f,    9.0000e-3f,    1.0000e-2f,    1.2500e-2f,    1.5000e-2f,
        1.7500e-2f,    2.0000e-2f,    2.2500e-2f,    2.5000e-2f,    2.7500e-2f,    3.0000e-2f,    3.5000e-2f,    4.0000e-2f,    4.5000e-2f,    5.0000e-2f,
        5.5000e-2f,    6.0000e-2f,    6.5000e-2f,    7.0000e-2f,    7.5000e-2f,    8.0000e-2f,    8.5000e-2f,    9.0000e-2f,    9.5000e-2f,    1.0000e-1f,
        1.2500e-1f,    1.5000e-1f,    1.7500e-1f,    2.0000e-1f,    2.2500e-1f,    2.5000e-1f,    2.7500e-1f,    3.0000e-1f,    3.5000e-1f,    4.0000e-1f,
        4.5000e-1f,    5.0000e-1f,    5.5000e-1f,    6.0000e-1f,    6.5000e-1f,    7.0000e-1f,    7.5000e-1f,    8.0000e-1f,    8.5000e-1f,    9.0000e-1f,
        9.5000e-1f,    1.0000e0f,    1.2500e0f,    1.5000e0f,    1.7500e0f,    2.0000e0f,    2.2500e0f,    2.5000e0f,    2.7500e0f,    3.0000e0f,
        3.5000e0f,    4.0000e0f,    4.5000e0f,    5.0000e0f,    5.5000e0f,    6.0000e0f,    6.5000e0f,    7.0000e0f,    7.5000e0f,    8.0000e0f,
        8.5000e0f,    9.0000e0f,    9.5000e0f,    1.0000e1f,    1.2500e1f,    1.5000e1f,    1.7500e1f,    2.0000e1f,    2.5000e1f,    2.7500e1f,
        3.0000e1f,    3.5000e1f,    4.0000e1f,    4.5000e1f,    5.0000e1f,    5.5000e1f,    6.0000e1f,    6.5000e1f,    7.0000e1f,    7.5000e1f,
        8.0000e1f,    8.5000e1f,    9.0000e1f,    9.5000e1f,    1.0000e2f,    1.2500e2f,    1.5000e2f,    1.7500e2f,    2.0000e2f,    2.2500e2f,
        2.5000e2f,    2.7500e2f,    3.0000e2f,    3.5000e2f,    4.0000e2f,    4.5000e2f,    5.0000e2f,    5.5000e2f,    6.0000e2f,    6.5000e2f,
        7.0000e2f,    7.5000e2f,    8.0000e2f,    8.5000e2f,    9.0000e2f,    9.5000e2f,    1.0000e3f,    1.5000e3f,    2.0000e3f,    2.5000e3f,
        3.0000e3f,    4.0000e3f,    5.0000e3f,    6.0000e3f,    7.0000e3f,    8.0000e3f,    9.0000e3f,    1.0000e4f,
        1.000e-03f,             1.500e-03f,             2.000e-03f,             2.500e-03f,             3.000e-03f,             4.000e-03f,             5.000e-03f,             6.000e-03f,             7.000e-03f,             8.000e-03f,
        9.000e-03f,             1.000e-02f,             1.250e-02f,             1.500e-02f,             1.750e-02f,             2.000e-02f,             2.250e-02f,             2.500e-02f,             2.750e-02f,             3.000e-02f,
        3.500e-02f,             4.000e-02f,             4.500e-02f,             5.000e-02f,             5.500e-02f,             6.000e-02f,             6.500e-02f,             7.000e-02f,             7.500e-02f,             8.000e-02f,
        8.500e-02f,             9.000e-02f,             9.500e-02f,             1.000e-01f,             1.250e-01f,             1.500e-01f,             1.750e-01f,             2.000e-01f,             2.250e-01f,             2.500e-01f,
        2.750e-01f,             3.000e-01f,             3.500e-01f,             4.000e-01f,             4.500e-01f,             5.000e-01f,             5.500e-01f,             6.000e-01f,             6.500e-01f,             7.000e-01f,
        7.500e-01f,             8.000e-01f,             8.500e-01f,             9.000e-01f,             9.500e-01f,             1.000e+00f,             1.250e+00f,             1.500e+00f,             1.750e+00f,             2.000e+00f,
        2.250e+00f,             2.500e+00f,             2.750e+00f,             3.000e+00f,             3.500e+00f,             4.000e+00f,             4.500e+00f,             5.000e+00f,             5.500e+00f,             6.000e+00f,
        6.500e+00f,             7.000e+00f,             7.500e+00f,             8.000e+00f,             8.500e+00f,             9.000e+00f,             9.500e+00f,             1.000e+01f,             1.250e+01f,             1.500e+01f,
        1.750e+01f,             2.000e+01f,             2.500e+01f,             2.750e+01f,             3.000e+01f,             3.500e+01f,             4.000e+01f,             4.500e+01f,             5.000e+01f,             5.500e+01f,
        6.000e+01f,             6.500e+01f,             7.000e+01f,             7.500e+01f,             8.000e+01f,             8.500e+01f,             9.000e+01f,             9.500e+01f,             1.000e+02f,             1.250e+02f,
        1.500e+02f,             1.750e+02f,             2.000e+02f,             2.250e+02f,             2.500e+02f,             2.750e+02f,             3.000e+02f,             3.500e+02f,             4.000e+02f,             4.500e+02f,
        5.000e+02f,             5.500e+02f,             6.000e+02f,             6.500e+02f,             7.000e+02f,             7.500e+02f,             8.000e+02f,             8.500e+02f,             9.000e+02f,             9.500e+02f,
        1.000e+03f,             1.500e+03f,             2.000e+03f,             2.500e+03f,             3.000e+03f,             4.000e+03f,             5.000e+03f,             6.000e+03f,             7.000e+03f,             8.000e+03f,
        9.000e+03f,             1.000e+04f,
        1.000e-03f,             1.500e-03f,             2.000e-03f,             2.500e-03f,             3.000e-03f,             4.000e-03f,             5.000e-03f,             6.000e-03f,             7.000e-03f,             8.000e-03f,
        9.000e-03f,             1.000e-02f,             1.250e-02f,             1.500e-02f,             1.750e-02f,             2.000e-02f,             2.250e-02f,             2.500e-02f,             2.750e-02f,             3.000e-02f,
        3.500e-02f,             4.000e-02f,             4.500e-02f,             5.000e-02f,             5.500e-02f,             6.000e-02f,             6.500e-02f,             7.000e-02f,             7.500e-02f,             8.000e-02f,
        8.500e-02f,             9.000e-02f,             9.500e-02f,             1.000e-01f,             1.250e-01f,             1.500e-01f,             1.750e-01f,             2.000e-01f,             2.250e-01f,             2.500e-01f,
        2.750e-01f,             3.000e-01f,             3.500e-01f,             4.000e-01f,             4.500e-01f,             5.000e-01f,             5.500e-01f,             6.000e-01f,             6.500e-01f,             7.000e-01f,
        7.500e-01f,             8.000e-01f,             8.500e-01f,             9.000e-01f,             9.500e-01f,             1.000e+00f,             1.250e+00f,             1.500e+00f,             1.750e+00f,             2.000e+00f,
        2.250e+00f,             2.500e+00f,             2.750e+00f,             3.000e+00f,             3.500e+00f,             4.000e+00f,             4.500e+00f,             5.000e+00f,             5.500e+00f,             6.000e+00f,
        6.500e+00f,             7.000e+00f,             7.500e+00f,             8.000e+00f,             8.500e+00f,             9.000e+00f,             9.500e+00f,             1.000e+01f,             1.250e+01f,             1.500e+01f,
        1.750e+01f,             2.000e+01f,             2.500e+01f,             2.750e+01f,             3.000e+01f,             3.500e+01f,             4.000e+01f,             4.500e+01f,             5.000e+01f,             5.500e+01f,
        6.000e+01f,             6.500e+01f,             7.000e+01f,             7.500e+01f,             8.000e+01f,             8.500e+01f,             9.000e+01f,             9.500e+01f,             1.000e+02f,             1.250e+02f,
        1.500e+02f,             1.750e+02f,             2.000e+02f,             2.250e+02f,             2.500e+02f,             2.750e+02f,             3.000e+02f,             3.500e+02f,             4.000e+02f,             4.500e+02f,
        5.000e+02f,             5.500e+02f,             6.000e+02f,             6.500e+02f,             7.000e+02f,             7.500e+02f,             8.000e+02f,             8.500e+02f,             9.000e+02f,             9.500e+02f,
        1.000e+03f,             1.500e+03f,             2.000e+03f,             2.500e+03f,             3.000e+03f,             4.000e+03f,             5.000e+03f,             6.000e+03f,             7.000e+03f,             8.000e+03f,
        9.000e+03f,             1.000e+04f,
		1.000e-03f,		1.500e-03f,		2.000e-03f,		2.500e-03f,		3.000e-03f,		4.000e-03f,		5.000e-03f,		6.000e-03f,		7.000e-03f,		8.000e-03f,
		9.000e-03f,		1.000e-02f,		1.250e-02f,		1.500e-02f,		1.750e-02f,		2.000e-02f,		2.250e-02f,		2.500e-02f,		2.750e-02f,		3.000e-02f,
		3.500e-02f,		4.000e-02f,		4.500e-02f,		5.000e-02f,		5.500e-02f,		6.000e-02f,		6.500e-02f,		7.000e-02f,		7.500e-02f,		8.000e-02f,
		8.500e-02f,		9.000e-02f,		9.500e-02f,		1.000e-01f,		1.250e-01f,		1.500e-01f,		1.750e-01f,		2.000e-01f,		2.250e-01f,		2.500e-01f,
		2.750e-01f,		3.000e-01f,		3.500e-01f,		4.000e-01f,		4.500e-01f,		5.000e-01f,		5.500e-01f,		6.000e-01f,		6.500e-01f,		7.000e-01f,
		7.500e-01f,		8.000e-01f,		8.500e-01f,		9.000e-01f,		9.500e-01f,		1.000e+00f,		1.250e+00f,		1.500e+00f,		1.750e+00f,		2.000e+00f,
		2.250e+00f,		2.500e+00f,		2.750e+00f,		3.000e+00f,		3.500e+00f,		4.000e+00f,		4.500e+00f,		5.000e+00f,		5.500e+00f,		6.000e+00f,
		6.500e+00f,		7.000e+00f,		7.500e+00f,		8.000e+00f,		8.500e+00f,		9.000e+00f,		9.500e+00f,		1.000e+01f,		1.250e+01f,		1.500e+01f,
		1.750e+01f,		2.000e+01f,		2.500e+01f,		2.750e+01f,		3.000e+01f,		3.500e+01f,		4.000e+01f,		4.500e+01f,		5.000e+01f,		5.500e+01f,
		6.000e+01f,		6.500e+01f,		7.000e+01f,		7.500e+01f,		8.000e+01f,		8.500e+01f,		9.000e+01f,		9.500e+01f,		1.000e+02f,		1.250e+02f,
		1.500e+02f,		1.750e+02f,		2.000e+02f,		2.250e+02f,		2.500e+02f,		2.750e+02f,		3.000e+02f,		3.500e+02f,		4.000e+02f,		4.500e+02f,
		5.000e+02f,		5.500e+02f,		6.000e+02f,		6.500e+02f,		7.000e+02f,		7.500e+02f,		8.000e+02f,		8.500e+02f,		9.000e+02f,		9.500e+02f,
		1.000e+03f,		1.500e+03f,		2.000e+03f,		2.500e+03f,		3.000e+03f,		4.000e+03f,		5.000e+03f,		6.000e+03f,		7.000e+03f,		8.000e+03f,
		9.000e+03f,		1.000e+04f

    },
    {
        1.3370e2f,    1.6380e2f,    1.8910e2f,    2.1140e2f,    2.3160e2f,    2.6750e2f,    2.9900e2f,    3.2760e2f,    3.5380e2f,    3.7820e2f,
        4.0120e2f,    4.2290e2f,    4.6600e2f,    5.0360e2f,    5.3720e2f,    5.6730e2f,    5.9460e2f,    6.1950e2f,    6.4210e2f,    6.6280e2f,
        6.9890e2f,    7.2900e2f,    7.5380e2f,    7.7400e2f,    7.9010e2f,    8.0260e2f,    8.1190e2f,    8.1830e2f,    8.2230e2f,    8.2410e2f,
        8.2390e2f,    8.2220e2f,    8.1900e2f,    8.1450e2f,    7.8010e2f,    7.3600e2f,    6.9590e2f,    6.6040e2f,    6.2860e2f,    5.9990e2f,
        5.7370e2f,    5.4970e2f,    5.0750e2f,    4.7140e2f,    4.4010e2f,    4.1280e2f,    3.8880e2f,    3.6760e2f,    3.4890e2f,    3.3220e2f,
        3.1720e2f,    3.0370e2f,    2.9140e2f,    2.8030e2f,    2.7000e2f,    2.6060e2f,    2.2280e2f,    1.9550e2f,    1.7480e2f,    1.5850e2f,
        1.4530e2f,    1.3430e2f,    1.2500e2f,    1.1710e2f,    1.0410e2f,    9.3980e1f,    8.5800e1f,    7.9060e1f,    7.3390e1f,    6.8540e1f,
        6.4340e1f,    6.0680e1f,    5.7440e1f,    5.4560e1f,    5.1990e1f,    4.9660e1f,    4.7560e1f,    4.5640e1f,    3.8130e1f,    3.2900e1f,
        2.9040e1f,    2.6050e1f,    2.1740e1f,    2.0120e1f,    1.8750e1f,    1.6560e1f,    1.4870e1f,    1.3530e1f,    1.2440e1f,    1.1540e1f,
        1.0780e1f,    1.0120e1f,    9.5550e0f,    9.0590e0f,    8.6220e0f,    8.2330e0f,    7.8840e0f,    7.5700e0f,    7.2860e0f,    6.1900e0f,
        5.4430e0f,    4.9010e0f,    4.4910e0f,    4.1690e0f,    3.9100e0f,    3.6970e0f,    3.5190e0f,    3.2400e0f,    3.0310e0f,    2.8700e0f,
        2.7430e0f,    2.6400e0f,    2.5550e0f,    2.4850e0f,    2.4260e0f,    2.3750e0f,    2.3330e0f,    2.2960e0f,    2.2640e0f,    2.2360e0f,
        2.2110e0f,    2.0700e0f,    2.0210e0f,    2.0040e0f,    2.0010e0f,    2.0120e0f,    2.0310e0f,    2.0520e0f,    2.0720e0f,    2.0910e0f,
        2.1090e0f,    2.1260e0f,    7.3510e1f,    9.0030e1f,    1.0400e2f,    1.1620e2f,    1.2730e2f,    1.4700e2f,    1.6440e2f,    1.8010e2f,
        1.9450e2f,    2.0790e2f,    2.2050e2f,    2.3250e2f,    2.5610e2f,    2.7680e2f,    2.9530e2f,    3.1190e2f,    3.2700e2f,    3.4080e2f,
        3.5340e2f,    3.6490e2f,    3.8530e2f,    4.0250e2f,    4.1700e2f,    4.2920e2f,    4.3930e2f,    4.4760e2f,    4.5430e2f,    4.5960e2f,
        4.6370e2f,    4.6660e2f,    4.6860e2f,    4.6980e2f,    4.7020e2f,    4.7000e2f,    4.6190e2f,    4.4720e2f,    4.3010e2f,    4.1270e2f,
        3.9650e2f,    3.8140e2f,    3.6750e2f,    3.5450e2f,    3.3110e2f,    3.1060e2f,    2.9260e2f,    2.7650e2f,    2.6220e2f,    2.4940e2f,
        2.3810e2f,    2.2790e2f,    2.1870e2f,    2.1020e2f,    2.0260e2f,    1.9570e2f,    1.8920e2f,    1.8320e2f,    1.5890e2f,    1.4090e2f,
        1.2700e2f,    1.1600e2f,    1.0690e2f,    9.9260e1f,    9.2790e1f,    8.7200e1f,    7.8030e1f,    7.0780e1f,    6.4900e1f,    6.0020e1f,
        5.5900e1f,    5.2360e1f,    4.9280e1f,    4.6580e1f,    4.4190e1f,    4.2060e1f,    4.0150e1f,    3.8420e1f,    3.6840e1f,    3.5410e1f,
        2.9760e1f,    2.5790e1f,    2.2840e1f,    2.0550e1f,    1.7220e1f,    1.5970e1f,    1.4910e1f,    1.3200e1f,    1.1880e1f,    1.0830e1f,
        9.9730e0f,    9.2610e0f,    8.6580e0f,    8.1410e0f,    7.6930e0f,    7.3000e0f,    6.9530e0f,    6.6440e0f,    6.3670e0f,    6.1180e0f,
        5.8910e0f,    5.0160e0f,    4.4190e0f,    3.9850e0f,    3.6560e0f,    3.3970e0f,    3.1890e0f,    3.0180e0f,    2.8750e0f,    2.6500e0f,
        2.4820e0f,    2.3520e0f,    2.2480e0f,    2.1640e0f,    2.0940e0f,    2.0360e0f,    1.9870e0f,    1.9460e0f,    1.9100e0f,    1.8790e0f,
        1.8520e0f,    1.8290e0f,    1.8090e0f,    1.6960e0f,    1.6600e0f,    1.6510e0f,    1.6530e0f,    1.6690e0f,    1.6900e0f,    1.7110e0f,
        1.7310e0f,    1.7490e0f,    1.7660e0f,    1.7820e0f,    9.2380e1f,    1.1310e2f,    1.3060e2f,    1.4610e2f,    1.6000e2f,    1.8480e2f,
        2.0660e2f,    2.2630e2f,    2.4440e2f,    2.6130e2f,    2.7710e2f,    2.9210e2f,    3.2060e2f,    3.4480e2f,    3.6570e2f,    3.8380e2f,
        3.9960e2f,    4.1320e2f,    4.2500e2f,    4.3510e2f,    4.5100e2f,    4.6200e2f,    4.6920e2f,    4.7340e2f,    4.7520e2f,    4.7510e2f,
        4.7370e2f,    4.7120e2f,    4.6800e2f,    4.6420e2f,    4.6010e2f,    4.5580e2f,    4.5130e2f,    4.4680e2f,    4.2450e2f,    4.0450e2f,
        3.8670e2f,    3.7100e2f,    3.5680e2f,    3.4400e2f,    3.3230e2f,    3.2150e2f,    3.0170e2f,    2.8420e2f,    2.6860e2f,    2.5480e2f,
        2.4250e2f,    2.3140e2f,    2.2150e2f,    2.1240e2f,    2.0420e2f,    1.9660e2f,    1.8970e2f,    1.8330e2f,    1.7740e2f,    1.7190e2f,
        1.4940e2f,    1.3270e2f,    1.1980e2f,    1.0940e2f,    1.0090e2f,    9.3770e1f,    8.7690e1f,    8.2450e1f,    7.3830e1f,    6.7030e1f,
        6.1510e1f,    5.6910e1f,    5.3030e1f,    4.9700e1f,    4.6810e1f,    4.4280e1f,    4.2030e1f,    4.0020e1f,    3.8220e1f,    3.6580e1f,
        3.5100e1f,    3.3750e1f,    2.8410e1f,    2.4650e1f,    2.1850e1f,    1.9680e1f,    1.6510e1f,    1.5320e1f,    1.4300e1f,    1.2670e1f,
        1.1410e1f,    1.0410e1f,    9.5900e0f,    8.9080e0f,    8.3300e0f,    7.8350e0f,    7.4050e0f,    7.0290e0f,    6.6960e0f,    6.3990e0f,
        6.1330e0f,    5.8930e0f,    5.6760e0f,    4.8350e0f,    4.2610e0f,    3.8430e0f,    3.5250e0f,    3.2760e0f,    3.0750e0f,    2.9100e0f,
        2.7720e0f,    2.5550e0f,    2.3920e0f,    2.2660e0f,    2.1660e0f,    2.0860e0f,    2.0190e0f,    1.9640e0f,    1.9180e0f,    1.8780e0f,
        1.8450e0f,    1.8160e0f,    1.7900e0f,    1.7680e0f,    1.7490e0f,    1.6470e0f,    1.6170e0f,    1.6130e0f,    1.6190e0f,    1.6420e0f,
        1.6680e0f,    1.6920e0f,    1.7140e0f,    1.7340e0f,    1.7520e0f,    1.7680e0f,    1.7490e2f,    2.1420e2f,    2.4740e2f,    2.7660e2f,
        3.0300e2f,    3.4590e2f,    3.8410e2f,    4.1870e2f,    4.5040e2f,    4.7900e2f,    5.0520e2f,    5.2980e2f,    5.7910e2f,    6.2080e2f,
        6.5710e2f,    6.8930e2f,    7.1790e2f,    7.4320e2f,    7.6570e2f,    7.8580e2f,    8.2080e2f,    8.5050e2f,    8.7500e2f,    8.9450e2f,
        9.0930e2f,    9.2010e2f,    9.2740e2f,    9.3190e2f,    9.3390e2f,    9.3370e2f,    9.3160e2f,    9.2800e2f,    9.2290e2f,    9.1680e2f,
        8.7440e2f,    8.2320e2f,    7.7130e2f,    7.2240e2f,    6.7620e2f,    6.3400e2f,    5.9620e2f,    5.6280e2f,    5.0720e2f,    4.6340e2f,
        4.2850e2f,    4.0020e2f,    3.7670e2f,    3.5610e2f,    3.3800e2f,    3.2190e2f,    3.0750e2f,    2.9450e2f,    2.8280e2f,    2.7200e2f,
        2.6210e2f,    2.5300e2f,    2.1660e2f,    1.9040e2f,    1.7030e2f,    1.5450e2f,    1.4170e2f,    1.3100e2f,    1.2200e2f,    1.1430e2f,
        1.0160e2f,    9.1730e1f,    8.3740e1f,    7.7140e1f,    7.1600e1f,    6.6860e1f,    6.2760e1f,    5.9180e1f,    5.6020e1f,    5.3210e1f,
        5.0700e1f,    4.8430e1f,    4.6370e1f,    4.4500e1f,    3.7170e1f,    3.2060e1f,    2.8290e1f,    2.5380e1f,    2.1170e1f,    1.9600e1f,
        1.8260e1f,    1.6120e1f,    1.4480e1f,    1.3180e1f,    1.2120e1f,    1.1230e1f,    1.0490e1f,    9.8540e0f,    9.3020e0f,    8.8190e0f,
        8.3930e0f,    8.0140e0f,    7.6750e0f,    7.3690e0f,    7.0930e0f,    6.0250e0f,    5.2980e0f,    4.7700e0f,    4.3710e0f,    4.0570e0f,
        3.8050e0f,    3.5980e0f,    3.4250e0f,    3.1530e0f,    2.9500e0f,    2.7930e0f,    2.6690e0f,    2.5690e0f,    2.4860e0f,    2.4180e0f,
        2.3600e0f,    2.3110e0f,    2.2690e0f,    2.2310e0f,    2.1990e0f,    2.1700e0f,    2.1450e0f,    2.0040e0f,    1.9540e0f,    1.9370e0f,
        1.9340e0f,    1.9450e0f,    1.9640e0f,    1.9850e0f,    2.0050e0f,    2.0240e0f,    2.0420e0f,    2.0590e0f,
        1.707e+02f,             2.091e+02f,             2.414e+02f,             2.699e+02f,             2.957e+02f,             3.387e+02f,             3.770e+02f,             4.116e+02f,             4.432e+02f,             4.722e+02f,
        4.990e+02f,             5.241e+02f,             5.739e+02f,             6.163e+02f,             6.534e+02f,             6.863e+02f,             7.156e+02f,             7.417e+02f,             7.648e+02f,             7.857e+02f,
        8.217e+02f,             8.517e+02f,             8.761e+02f,             8.955e+02f,             9.101e+02f,             9.208e+02f,             9.281e+02f,             9.325e+02f,             9.343e+02f,             9.340e+02f,
        9.319e+02f,             9.281e+02f,             9.230e+02f,             9.167e+02f,             8.740e+02f,             8.226e+02f,             7.705e+02f,             7.214e+02f,             6.766e+02f,             6.362e+02f,
        6.004e+02f,             5.689e+02f,             5.166e+02f,             4.750e+02f,             4.407e+02f,             4.119e+02f,             3.872e+02f,             3.655e+02f,             3.466e+02f,             3.298e+02f,
        3.147e+02f,             3.012e+02f,             2.890e+02f,             2.778e+02f,             2.676e+02f,             2.581e+02f,             2.206e+02f,             1.936e+02f,             1.730e+02f,             1.568e+02f,
        1.437e+02f,             1.328e+02f,             1.237e+02f,             1.158e+02f,             1.029e+02f,             9.284e+01f,             8.473e+01f,             7.803e+01f,             7.241e+01f,             6.760e+01f,
        6.345e+01f,             5.981e+01f,             5.661e+01f,             5.376e+01f,             5.122e+01f,             4.892e+01f,             4.684e+01f,             4.495e+01f,             3.753e+01f,             3.236e+01f,
        2.855e+01f,             2.561e+01f,             2.135e+01f,             1.976e+01f,             1.841e+01f,             1.625e+01f,             1.459e+01f,             1.328e+01f,             1.221e+01f,             1.132e+01f,
        1.057e+01f,             9.925e+00f,             9.369e+00f,             8.881e+00f,             8.452e+00f,             8.070e+00f,             7.727e+00f,             7.419e+00f,             7.140e+00f,             6.063e+00f,
        5.330e+00f,             4.798e+00f,             4.395e+00f,             4.078e+00f,             3.824e+00f,             3.615e+00f,             3.440e+00f,             3.166e+00f,             2.960e+00f,             2.802e+00f,
        2.675e+00f,             2.574e+00f,             2.491e+00f,             2.420e+00f,             2.362e+00f,             2.312e+00f,             2.269e+00f,             2.232e+00f,             2.200e+00f,             2.173e+00f,
        2.148e+00f,             2.015e+00f,             1.976e+00f,             1.967e+00f,             1.973e+00f,             2.000e+00f,             2.032e+00f,             2.064e+00f,             2.094e+00f,             2.122e+00f,
        2.148e+00f,             2.173e+00f,
        8.087e+01f,             9.904e+01f,             1.144e+02f,             1.279e+02f,             1.401e+02f,             1.617e+02f,             1.808e+02f,             1.981e+02f,             2.140e+02f,             2.287e+02f,
        2.426e+02f,             2.557e+02f,             2.819e+02f,             3.049e+02f,             3.254e+02f,             3.440e+02f,             3.609e+02f,             3.765e+02f,             3.908e+02f,             4.040e+02f,
        4.276e+02f,             4.478e+02f,             4.653e+02f,             4.803e+02f,             4.931e+02f,             5.041e+02f,             5.133e+02f,             5.211e+02f,             5.274e+02f,             5.325e+02f,
        5.365e+02f,             5.396e+02f,             5.417e+02f,             5.430e+02f,             5.403e+02f,             5.274e+02f,             5.092e+02f,             4.888e+02f,             4.674e+02f,             4.465e+02f,
        4.269e+02f,             4.087e+02f,             3.769e+02f,             3.504e+02f,             3.280e+02f,             3.091e+02f,             2.927e+02f,             2.781e+02f,             2.651e+02f,             2.533e+02f,
        2.425e+02f,             2.328e+02f,             2.239e+02f,             2.156e+02f,             2.082e+02f,             2.012e+02f,             1.731e+02f,             1.525e+02f,             1.368e+02f,             1.244e+02f,
        1.142e+02f,             1.057e+02f,             9.853e+01f,             9.238e+01f,             8.233e+01f,             7.443e+01f,             6.804e+01f,             6.276e+01f,             5.831e+01f,             5.451e+01f,
        5.122e+01f,             4.834e+01f,             4.579e+01f,             4.352e+01f,             4.149e+01f,             3.966e+01f,             3.799e+01f,             3.648e+01f,             3.053e+01f,             2.638e+01f,
        2.330e+01f,             2.093e+01f,             1.748e+01f,             1.619e+01f,             1.510e+01f,             1.334e+01f,             1.199e+01f,             1.092e+01f,             1.004e+01f,             9.317e+00f,
        8.704e+00f,             8.178e+00f,             7.723e+00f,             7.324e+00f,             6.972e+00f,             6.659e+00f,             6.379e+00f,             6.126e+00f,             5.897e+00f,             5.013e+00f,
        4.411e+00f,             3.974e+00f,             3.642e+00f,             3.382e+00f,             3.173e+00f,             3.001e+00f,             2.858e+00f,             2.632e+00f,             2.463e+00f,             2.332e+00f,
        2.228e+00f,             2.144e+00f,             2.075e+00f,             2.017e+00f,             1.968e+00f,             1.927e+00f,             1.891e+00f,             1.860e+00f,             1.833e+00f,             1.810e+00f,
        1.789e+00f,             1.674e+00f,             1.635e+00f,             1.622e+00f,             1.621e+00f,             1.631e+00f,             1.648e+00f,             1.665e+00f,             1.682e+00f,             1.697e+00f,
        1.712e+00f,             1.726e+00f,
		1.197e+02f,		1.467e+02f,		1.693e+02f,		1.893e+02f,		2.074e+02f,		2.395e+02f,		2.678e+02f,		2.933e+02f,		3.168e+02f,		3.387e+02f,
		3.592e+02f,		3.787e+02f,		4.170e+02f,		4.504e+02f,		4.801e+02f,		5.067e+02f,		5.307e+02f,		5.526e+02f,		5.724e+02f,		5.905e+02f,
		6.221e+02f,		6.483e+02f,		6.700e+02f,		6.877e+02f,		7.020e+02f,		7.132e+02f,		7.217e+02f,		7.278e+02f,		7.319e+02f,		7.341e+02f,
		7.348e+02f,		7.340e+02f,		7.320e+02f,		7.290e+02f,		7.029e+02f,		6.672e+02f,		6.291e+02f,		5.922e+02f,		5.583e+02f,		5.278e+02f,
		5.006e+02f,		4.763e+02f,		4.349e+02f,		4.012e+02f,		3.733e+02f,		3.498e+02f,		3.297e+02f,		3.121e+02f,		2.964e+02f,		2.824e+02f,
		2.699e+02f,		2.587e+02f,		2.485e+02f,		2.391e+02f,		2.306e+02f,		2.227e+02f,		1.911e+02f,		1.682e+02f,		1.508e+02f,		1.370e+02f,
		1.258e+02f,		1.164e+02f,		1.085e+02f,		1.017e+02f,		9.063e+01f,		8.192e+01f,		7.488e+01f,		6.905e+01f,		6.414e+01f,		5.994e+01f,
		5.630e+01f,		5.312e+01f,		5.031e+01f,		4.781e+01f,		4.557e+01f,		4.355e+01f,		4.171e+01f,		4.004e+01f,		3.349e+01f,		2.892e+01f,
		2.554e+01f,		2.293e+01f,		1.914e+01f,		1.773e+01f,		1.652e+01f,		1.460e+01f,		1.312e+01f,		1.194e+01f,		1.098e+01f,		1.019e+01f,
		9.514e+00f,		8.938e+00f,		8.440e+00f,		8.003e+00f,		7.618e+00f,		7.275e+00f,		6.968e+00f,		6.691e+00f,		6.441e+00f,		5.474e+00f,
		4.815e+00f,		4.337e+00f,		3.975e+00f,		3.690e+00f,		3.462e+00f,		3.274e+00f,		3.117e+00f,		2.870e+00f,		2.686e+00f,		2.544e+00f,
		2.431e+00f,		2.340e+00f,		2.265e+00f,		2.203e+00f,		2.151e+00f,		2.107e+00f,		2.069e+00f,		2.036e+00f,		2.008e+00f,		1.984e+00f,
		1.962e+00f,		1.850e+00f,		1.820e+00f,		1.818e+00f,		1.828e+00f,		1.861e+00f,		1.898e+00f,		1.933e+00f,		1.967e+00f,		1.998e+00f,
		2.026e+00f,		2.052e+00f
    },
    {
        4.3150e1f,    3.4600e1f,    2.9270e1f,    2.5570e1f,    2.2810e1f,    1.8940e1f,    1.6310e1f,    1.4390e1f,    1.2920e1f,    1.1750e1f,
        1.0800e1f,    1.0000e1f,    8.4850e0f,    7.4000e0f,    6.5810e0f,    5.9390e0f,    5.4210e0f,    4.9930e0f,    4.6330e0f,    4.3250e0f,
        3.8260e0f,    3.4370e0f,    3.1260e0f,    2.8700e0f,    2.6550e0f,    2.4730e0f,    2.3160e0f,    2.1780e0f,    2.0580e0f,    1.9510e0f,
        1.8550e0f,    1.7690e0f,    1.6910e0f,    1.6200e0f,    1.3430e0f,    1.1520e0f,    1.0100e0f,    9.0160e-1f,    8.1520e-1f,    7.4470e-1f,
        6.8550e-1f,    6.3510e-1f,    5.5450e-1f,    4.9280e-1f,    4.4390e-1f,    4.0430e-1f,    3.7150e-1f,    3.4380e-1f,    3.2010e-1f,    2.9960e-1f,
        2.8170e-1f,    2.6580e-1f,    2.5160e-1f,    2.3900e-1f,    2.2760e-1f,    2.1730e-1f,    1.7750e-1f,    1.5040e-1f,    1.3070e-1f,    1.1570e-1f,
        1.0380e-1f,    9.4300e-2f,    8.6400e-2f,    7.9700e-2f,    6.9200e-2f,    6.1100e-2f,    5.4800e-2f,    4.9700e-2f,    4.5500e-2f,    4.2000e-2f,
        3.8900e-2f,    3.6300e-2f,    3.4100e-2f,    3.2100e-2f,    3.0300e-2f,    2.8700e-2f,    2.7300e-2f,    2.6000e-2f,    2.1100e-2f,    1.7800e-2f,
        1.5400e-2f,    1.3600e-2f,    1.1000e-2f,    1.0000e-2f,    9.2000e-3f,    8.0000e-3f,    7.0000e-3f,    6.3000e-3f,    5.7000e-3f,    5.2000e-3f,
        4.8000e-3f,    4.4000e-3f,    4.1000e-3f,    3.9000e-3f,    3.6000e-3f,    3.4000e-3f,    3.3000e-3f,    3.1000e-3f,    2.9000e-3f,    2.4000e-3f,
        2.0000e-3f,    1.7000e-3f,    1.5000e-3f,    1.4000e-3f,    1.2000e-3f,    1.1000e-3f,    1.0000e-3f,    9.0000e-4f,    8.0000e-4f,    7.0000e-4f,
        6.0000e-4f,    6.0000e-4f,    5.0000e-4f,    5.0000e-4f,    5.0000e-4f,    4.0000e-4f,    4.0000e-4f,    4.0000e-4f,    4.0000e-4f,    3.0000e-4f,
        3.0000e-4f,    2.0000e-4f,    2.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    1.5780e1f,    1.3780e1f,    1.2310e1f,    1.1180e1f,    1.0270e1f,    8.9000e0f,    7.9000e0f,    7.1310e0f,
        6.5180e0f,    6.0160e0f,    5.5950e0f,    5.2370e0f,    4.5340e0f,    4.0150e0f,    3.6140e0f,    3.2940e0f,    3.0310e0f,    2.8110e0f,
        2.6240e0f,    2.4620e0f,    2.1980e0f,    1.9890e0f,    1.8190e0f,    1.6790e0f,    1.5600e0f,    1.4590e0f,    1.3710e0f,    1.2930e0f,
        1.2250e0f,    1.1640e0f,    1.1100e0f,    1.0610e0f,    1.0160e0f,    9.7510e-1f,    8.1480e-1f,    7.0260e-1f,    6.1930e-1f,    5.5480e-1f,
        5.0320e-1f,    4.6100e-1f,    4.2580e-1f,    3.9590e-1f,    3.4780e-1f,    3.1070e-1f,    2.8120e-1f,    2.5710e-1f,    2.3710e-1f,    2.2010e-1f,
        2.0550e-1f,    1.9280e-1f,    1.8170e-1f,    1.7190e-1f,    1.6310e-1f,    1.5520e-1f,    1.4800e-1f,    1.4160e-1f,    1.1650e-1f,    9.9200e-2f,
        8.6600e-2f,    7.6900e-2f,    6.9300e-2f,    6.3000e-2f,    5.7900e-2f,    5.3600e-2f,    4.6600e-2f,    4.1300e-2f,    3.7200e-2f,    3.3800e-2f,
        3.1000e-2f,    2.8600e-2f,    2.6600e-2f,    2.4800e-2f,    2.3300e-2f,    2.2000e-2f,    2.0800e-2f,    1.9700e-2f,    1.8800e-2f,    1.7900e-2f,
        1.4600e-2f,    1.2300e-2f,    1.0700e-2f,    9.4000e-3f,    7.6000e-3f,    7.0000e-3f,    6.4000e-3f,    5.6000e-3f,    4.9000e-3f,    4.4000e-3f,
        4.0000e-3f,    3.6000e-3f,    3.4000e-3f,    3.1000e-3f,    2.9000e-3f,    2.7000e-3f,    2.6000e-3f,    2.4000e-3f,    2.3000e-3f,    2.2000e-3f,
        2.1000e-3f,    1.7000e-3f,    1.4000e-3f,    1.2000e-3f,    1.1000e-3f,    1.0000e-3f,    9.0000e-4f,    8.0000e-4f,    7.0000e-4f,    6.0000e-4f,
        6.0000e-4f,    5.0000e-4f,    4.0000e-4f,    4.0000e-4f,    4.0000e-4f,    4.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,
        3.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    1.1970e1f,    1.0720e1f,    9.7490e0f,    8.9670e0f,    8.3240e0f,    7.3220e0f,
        6.5710e0f,    5.9820e0f,    5.5050e0f,    5.1100e0f,    4.7750e0f,    4.4880e0f,    3.9170e0f,    3.4910e0f,    3.1570e0f,    2.8890e0f,
        2.6670e0f,    2.4800e0f,    2.3210e0f,    2.1830e0f,    1.9550e0f,    1.7740e0f,    1.6270e0f,    1.5040e0f,    1.4010e0f,    1.3110e0f,
        1.2340e0f,    1.1660e0f,    1.1060e0f,    1.0520e0f,    1.0030e0f,    9.5970e-1f,    9.2000e-1f,    8.8370e-1f,    7.4060e-1f,    6.4000e-1f,
        5.6510e-1f,    5.0700e-1f,    4.6040e-1f,    4.2220e-1f,    3.9030e-1f,    3.6320e-1f,    3.1950e-1f,    2.8580e-1f,    2.5890e-1f,    2.3690e-1f,
        2.1850e-1f,    2.0300e-1f,    1.8970e-1f,    1.7810e-1f,    1.6790e-1f,    1.5890e-1f,    1.5090e-1f,    1.4360e-1f,    1.3710e-1f,    1.3120e-1f,
        1.0820e-1f,    9.2300e-2f,    8.0700e-2f,    7.1700e-2f,    6.4700e-2f,    5.8900e-2f,    5.4100e-2f,    5.0100e-2f,    4.3700e-2f,    3.8800e-2f,
        3.4900e-2f,    3.1700e-2f,    2.9100e-2f,    2.6900e-2f,    2.5000e-2f,    2.3400e-2f,    2.2000e-2f,    2.0800e-2f,    1.9600e-2f,    1.8600e-2f,
        1.7700e-2f,    1.6900e-2f,    1.3800e-2f,    1.1700e-2f,    1.0100e-2f,    8.9000e-3f,    7.3000e-3f,    6.6000e-3f,    6.1000e-3f,    5.3000e-3f,
        4.7000e-3f,    4.2000e-3f,    3.8000e-3f,    3.5000e-3f,    3.2000e-3f,    3.0000e-3f,    2.8000e-3f,    2.6000e-3f,    2.4000e-3f,    2.3000e-3f,
        2.2000e-3f,    2.1000e-3f,    2.0000e-3f,    1.6000e-3f,    1.3000e-3f,    1.2000e-3f,    1.0000e-3f,    9.0000e-4f,    8.0000e-4f,    8.0000e-4f,
        7.0000e-4f,    6.0000e-4f,    5.0000e-4f,    5.0000e-4f,    4.0000e-4f,    4.0000e-4f,    4.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,
        3.0000e-4f,    3.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    3.9730e1f,    3.2030e1f,    2.7180e1f,    2.3790e1f,
        2.1250e1f,    1.7660e1f,    1.5220e1f,    1.3440e1f,    1.2070e1f,    1.0980e1f,    1.0090e1f,    9.3420e0f,    7.9240e0f,    6.9100e0f,
        6.1450e0f,    5.5450e0f,    5.0610e0f,    4.6610e0f,    4.3240e0f,    4.0360e0f,    3.5700e0f,    3.2070e0f,    2.9160e0f,    2.6770e0f,
        2.4770e0f,    2.3070e0f,    2.1600e0f,    2.0320e0f,    1.9190e0f,    1.8190e0f,    1.7300e0f,    1.6500e0f,    1.5770e0f,    1.5110e0f,
        1.2530e0f,    1.0740e0f,    9.4190e-1f,    8.4050e-1f,    7.5990e-1f,    6.9420e-1f,    6.3920e-1f,    5.9230e-1f,    5.1740e-1f,    4.6000e-1f,
        4.1460e-1f,    3.7770e-1f,    3.4720e-1f,    3.2130e-1f,    2.9910e-1f,    2.7990e-1f,    2.6310e-1f,    2.4830e-1f,    2.3510e-1f,    2.2330e-1f,
        2.1270e-1f,    2.0300e-1f,    1.6590e-1f,    1.4050e-1f,    1.2210e-1f,    1.0810e-1f,    9.7000e-2f,    8.8000e-2f,    8.0600e-2f,    7.4400e-2f,
        6.4600e-2f,    5.7000e-2f,    5.1100e-2f,    4.6400e-2f,    4.2400e-2f,    3.9100e-2f,    3.6300e-2f,    3.3900e-2f,    3.1800e-2f,    2.9900e-2f,
        2.8200e-2f,    2.6800e-2f,    2.5400e-2f,    2.4200e-2f,    1.9700e-2f,    1.6600e-2f,    1.4300e-2f,    1.2600e-2f,    1.0200e-2f,    9.3000e-3f,
        8.6000e-3f,    7.4000e-3f,    6.5000e-3f,    5.8000e-3f,    5.3000e-3f,    4.8000e-3f,    4.5000e-3f,    4.1000e-3f,    3.8000e-3f,    3.6000e-3f,
        3.4000e-3f,    3.2000e-3f,    3.0000e-3f,    2.9000e-3f,    2.7000e-3f,    2.2000e-3f,    1.9000e-3f,    1.6000e-3f,    1.4000e-3f,    1.3000e-3f,
        1.1000e-3f,    1.0000e-3f,    1.0000e-3f,    8.0000e-4f,    7.0000e-4f,    7.0000e-4f,    6.0000e-4f,    5.0000e-4f,    5.0000e-4f,    5.0000e-4f,
        4.0000e-4f,    4.0000e-4f,    4.0000e-4f,    4.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    2.0000e-4f,    2.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        3.881e+01f,             3.135e+01f,             2.664e+01f,             2.332e+01f,             2.085e+01f,             1.735e+01f,             1.497e+01f,             1.322e+01f,             1.188e+01f,             1.081e+01f,
        9.935e+00f,             9.205e+00f,             7.813e+00f,             6.817e+00f,             6.064e+00f,             5.474e+00f,             4.997e+00f,             4.604e+00f,             4.271e+00f,             3.988e+00f,
        3.528e+00f,             3.171e+00f,             2.884e+00f,             2.648e+00f,             2.450e+00f,             2.282e+00f,             2.137e+00f,             2.011e+00f,             1.899e+00f,             1.800e+00f,
        1.712e+00f,             1.633e+00f,             1.561e+00f,             1.496e+00f,             1.240e+00f,             1.063e+00f,             9.331e-01f,             8.327e-01f,             7.529e-01f,             6.879e-01f,
        6.335e-01f,             5.870e-01f,             5.128e-01f,             4.560e-01f,             4.110e-01f,             3.746e-01f,             3.443e-01f,             3.186e-01f,             2.967e-01f,             2.776e-01f,
        2.611e-01f,             2.463e-01f,             2.333e-01f,             2.216e-01f,             2.110e-01f,             2.015e-01f,             1.647e-01f,             1.396e-01f,             1.213e-01f,             1.073e-01f,
        9.636e-02f,             8.748e-02f,             8.014e-02f,             7.397e-02f,             6.416e-02f,             5.670e-02f,             5.084e-02f,             4.610e-02f,             4.219e-02f,             3.891e-02f,
        3.611e-02f,             3.369e-02f,             3.159e-02f,             2.974e-02f,             2.810e-02f,             2.663e-02f,             2.532e-02f,             2.413e-02f,             1.956e-02f,             1.648e-02f,
        1.424e-02f,             1.256e-02f,             1.017e-02f,             9.292e-03f,             8.557e-03f,             7.393e-03f,             6.514e-03f,             5.825e-03f,             5.270e-03f,             4.814e-03f,
        4.431e-03f,             4.107e-03f,             3.827e-03f,             3.584e-03f,             3.371e-03f,             3.182e-03f,             3.013e-03f,             2.862e-03f,             2.726e-03f,             2.204e-03f,
        1.853e-03f,             1.600e-03f,             1.409e-03f,             1.259e-03f,             1.139e-03f,             1.041e-03f,             9.579e-04f,             8.273e-04f,             7.286e-04f,             6.514e-04f,
        5.894e-04f,             5.384e-04f,             4.957e-04f,             4.595e-04f,             4.283e-04f,             4.012e-04f,             3.773e-04f,             3.563e-04f,             3.375e-04f,             3.207e-04f,
        3.055e-04f,             2.083e-04f,             1.588e-04f,             1.287e-04f,             1.084e-04f,             8.277e-05f,             6.714e-05f,             5.660e-05f,             4.899e-05f,             4.323e-05f,
        3.871e-05f,             3.508e-05f,
        1.991e+01f,             1.669e+01f,             1.452e+01f,             1.293e+01f,             1.171e+01f,             9.926e+00f,             8.675e+00f,             7.739e+00f,             7.008e+00f,             6.418e+00f,
        5.931e+00f,             5.521e+00f,             4.728e+00f,             4.153e+00f,             3.714e+00f,             3.367e+00f,             3.085e+00f,             2.850e+00f,             2.652e+00f,             2.482e+00f,
        2.204e+00f,             1.987e+00f,             1.812e+00f,             1.667e+00f,             1.546e+00f,             1.442e+00f,             1.353e+00f,             1.275e+00f,             1.206e+00f,             1.144e+00f,
        1.089e+00f,             1.040e+00f,             9.948e-01f,             9.539e-01f,             7.939e-01f,             6.824e-01f,             6.001e-01f,             5.365e-01f,             4.858e-01f,             4.444e-01f,
        4.097e-01f,             3.804e-01f,             3.333e-01f,             2.971e-01f,             2.684e-01f,             2.449e-01f,             2.255e-01f,             2.090e-01f,             1.949e-01f,             1.827e-01f,
        1.720e-01f,             1.625e-01f,             1.540e-01f,             1.464e-01f,             1.396e-01f,             1.334e-01f,             1.094e-01f,             9.290e-02f,             8.089e-02f,             7.172e-02f,
        6.448e-02f,             5.861e-02f,             5.375e-02f,             4.966e-02f,             4.315e-02f,             3.819e-02f,             3.428e-02f,             3.112e-02f,             2.850e-02f,             2.630e-02f,
        2.443e-02f,             2.281e-02f,             2.140e-02f,             2.016e-02f,             1.905e-02f,             1.807e-02f,             1.718e-02f,             1.638e-02f,             1.330e-02f,             1.121e-02f,
        9.706e-03f,             8.561e-03f,             6.939e-03f,             6.343e-03f,             5.843e-03f,             5.051e-03f,             4.452e-03f,             3.983e-03f,             3.604e-03f,             3.293e-03f,
        3.032e-03f,             2.811e-03f,             2.620e-03f,             2.454e-03f,             2.308e-03f,             2.179e-03f,             2.064e-03f,             1.960e-03f,             1.867e-03f,             1.510e-03f,
        1.270e-03f,             1.097e-03f,             9.661e-04f,             8.637e-04f,             7.814e-04f,             7.137e-04f,             6.571e-04f,             5.675e-04f,             4.999e-04f,             4.470e-04f,
        4.045e-04f,             3.695e-04f,             3.403e-04f,             3.154e-04f,             2.940e-04f,             2.754e-04f,             2.591e-04f,             2.446e-04f,             2.318e-04f,             2.202e-04f,
        2.098e-04f,             1.431e-04f,             1.091e-04f,             8.845e-05f,             7.453e-05f,             5.690e-05f,             4.616e-05f,             3.892e-05f,             3.369e-05f,             2.973e-05f,
        2.663e-05f,             2.413e-05f,
		2.163e+01f,		1.840e+01f,		1.614e+01f,		1.446e+01f,		1.314e+01f,		1.120e+01f,		9.825e+00f,		8.786e+00f,		7.970e+00f,		7.310e+00f,
		6.762e+00f,		6.300e+00f,		5.404e+00f,		4.751e+00f,		4.253e+00f,		3.858e+00f,		3.536e+00f,		3.269e+00f,		3.042e+00f,		2.848e+00f,
		2.531e+00f,		2.282e+00f,		2.082e+00f,		1.917e+00f,		1.777e+00f,		1.659e+00f,		1.556e+00f,		1.466e+00f,		1.387e+00f,		1.316e+00f,
		1.253e+00f,		1.196e+00f,		1.145e+00f,		1.098e+00f,		9.142e-01f,		7.861e-01f,		6.914e-01f,		6.183e-01f,		5.600e-01f,		5.124e-01f,
		4.727e-01f,		4.390e-01f,		3.850e-01f,		3.435e-01f,		3.105e-01f,		2.836e-01f,		2.612e-01f,		2.423e-01f,		2.261e-01f,		2.119e-01f,
		1.995e-01f,		1.885e-01f,		1.787e-01f,		1.699e-01f,		1.620e-01f,		1.548e-01f,		1.270e-01f,		1.080e-01f,		9.404e-02f,		8.340e-02f,
		7.500e-02f,		6.818e-02f,		6.254e-02f,		5.778e-02f,		5.021e-02f,		4.444e-02f,		3.989e-02f,		3.621e-02f,		3.317e-02f,		3.061e-02f,
		2.843e-02f,		2.654e-02f,		2.490e-02f,		2.345e-02f,		2.217e-02f,		2.102e-02f,		1.999e-02f,		1.905e-02f,		1.547e-02f,		1.304e-02f,
		1.128e-02f,		9.953e-03f,		8.066e-03f,		7.372e-03f,		6.790e-03f,		5.870e-03f,		5.173e-03f,		4.627e-03f,		4.187e-03f,		3.826e-03f,
		3.523e-03f,		3.265e-03f,		3.043e-03f,		2.850e-03f,		2.681e-03f,		2.531e-03f,		2.397e-03f,		2.277e-03f,		2.169e-03f,		1.754e-03f,
		1.475e-03f,		1.274e-03f,		1.122e-03f,		1.003e-03f,		9.075e-04f,		8.289e-04f,		7.631e-04f,		6.591e-04f,		5.806e-04f,		5.191e-04f,
		4.697e-04f,		4.291e-04f,		3.951e-04f,		3.663e-04f,		3.414e-04f,		3.198e-04f,		3.009e-04f,		2.841e-04f,		2.691e-04f,		2.557e-04f,
		2.436e-04f,		1.661e-04f,		1.267e-04f,		1.027e-04f,		8.655e-05f,		6.608e-05f,		5.362e-05f,		4.520e-05f,		3.913e-05f,		3.453e-05f,
		3.093e-05f,		2.803e-05f
    },
    {
        1.7690e2f,    1.9840e2f,    2.1840e2f,    2.3700e2f,    2.5440e2f,    2.8640e2f,    3.1530e2f,    3.4200e2f,    3.6670e2f,    3.9000e2f,
        4.1200e2f,    4.3290e2f,    4.7450e2f,    5.1100e2f,    5.4370e2f,    5.7330e2f,    6.0010e2f,    6.2450e2f,    6.4670e2f,    6.6710e2f,
        7.0280e2f,    7.3240e2f,    7.5690e2f,    7.7680e2f,    7.9270e2f,    8.0500e2f,    8.1420e2f,    8.2050e2f,    8.2430e2f,    8.2600e2f,
        8.2580e2f,    8.2390e2f,    8.2060e2f,    8.1610e2f,    7.8140e2f,    7.3710e2f,    6.9690e2f,    6.6130e2f,    6.2940e2f,    6.0060e2f,
        5.7440e2f,    5.5040e2f,    5.0800e2f,    4.7190e2f,    4.4060e2f,    4.1320e2f,    3.8910e2f,    3.6800e2f,    3.4920e2f,    3.3250e2f,
        3.1750e2f,    3.0390e2f,    2.9170e2f,    2.8050e2f,    2.7020e2f,    2.6080e2f,    2.2290e2f,    1.9570e2f,    1.7490e2f,    1.5860e2f,
        1.4540e2f,    1.3440e2f,    1.2510e2f,    1.1720e2f,    1.0420e2f,    9.4040e1f,    8.5860e1f,    7.9110e1f,    7.3430e1f,    6.8580e1f,
        6.4380e1f,    6.0710e1f,    5.7470e1f,    5.4600e1f,    5.2020e1f,    4.9690e1f,    4.7590e1f,    4.5670e1f,    3.8150e1f,    3.2920e1f,
        2.9050e1f,    2.6070e1f,    2.1750e1f,    2.0130e1f,    1.8760e1f,    1.6560e1f,    1.4880e1f,    1.3540e1f,    1.2450e1f,    1.1540e1f,
        1.0780e1f,    1.0130e1f,    9.5590e0f,    9.0630e0f,    8.6250e0f,    8.2360e0f,    7.8880e0f,    7.5730e0f,    7.2890e0f,    6.1920e0f,
        5.4450e0f,    4.9030e0f,    4.4920e0f,    4.1700e0f,    3.9110e0f,    3.6980e0f,    3.5200e0f,    3.2410e0f,    3.0320e0f,    2.8710e0f,
        2.7430e0f,    2.6400e0f,    2.5560e0f,    2.4850e0f,    2.4260e0f,    2.3760e0f,    2.3330e0f,    2.2960e0f,    2.2640e0f,    2.2360e0f,
        2.2110e0f,    2.0700e0f,    2.0210e0f,    2.0040e0f,    2.0010e0f,    2.0120e0f,    2.0310e0f,    2.0520e0f,    2.0720e0f,    2.0910e0f,
        2.1090e0f,    2.1260e0f,    8.9290e1f,    1.0380e2f,    1.1630e2f,    1.2740e2f,    1.3760e2f,    1.5590e2f,    1.7230e2f,    1.8720e2f,
        2.0100e2f,    2.1390e2f,    2.2610e2f,    2.3770e2f,    2.6070e2f,    2.8090e2f,    2.9890e2f,    3.1520e2f,    3.3000e2f,    3.4360e2f,
        3.5600e2f,    3.6740e2f,    3.8750e2f,    4.0450e2f,    4.1880e2f,    4.3080e2f,    4.4080e2f,    4.4900e2f,    4.5570e2f,    4.6090e2f,
        4.6490e2f,    4.6780e2f,    4.6970e2f,    4.7080e2f,    4.7120e2f,    4.7090e2f,    4.6270e2f,    4.4790e2f,    4.3070e2f,    4.1330e2f,
        3.9700e2f,    3.8190e2f,    3.6790e2f,    3.5490e2f,    3.3140e2f,    3.1090e2f,    2.9280e2f,    2.7680e2f,    2.6250e2f,    2.4970e2f,
        2.3830e2f,    2.2810e2f,    2.1880e2f,    2.1040e2f,    2.0280e2f,    1.9580e2f,    1.8940e2f,    1.8340e2f,    1.5900e2f,    1.4100e2f,
        1.2710e2f,    1.1600e2f,    1.0690e2f,    9.9330e1f,    9.2850e1f,    8.7260e1f,    7.8070e1f,    7.0820e1f,    6.4940e1f,    6.0050e1f,
        5.5930e1f,    5.2380e1f,    4.9310e1f,    4.6610e1f,    4.4220e1f,    4.2080e1f,    4.0170e1f,    3.8440e1f,    3.6860e1f,    3.5430e1f,
        2.9770e1f,    2.5800e1f,    2.2850e1f,    2.0560e1f,    1.7230e1f,    1.5980e1f,    1.4910e1f,    1.3200e1f,    1.1880e1f,    1.0830e1f,
        9.9770e0f,    9.2650e0f,    8.6620e0f,    8.1450e0f,    7.6960e0f,    7.3030e0f,    6.9560e0f,    6.6470e0f,    6.3700e0f,    6.1200e0f,
        5.8930e0f,    5.0180e0f,    4.4210e0f,    3.9870e0f,    3.6570e0f,    3.3980e0f,    3.1900e0f,    3.0190e0f,    2.8760e0f,    2.6510e0f,
        2.4830e0f,    2.3530e0f,    2.2480e0f,    2.1640e0f,    2.0940e0f,    2.0360e0f,    1.9880e0f,    1.9460e0f,    1.9100e0f,    1.8790e0f,
        1.8530e0f,    1.8290e0f,    1.8090e0f,    1.6970e0f,    1.6600e0f,    1.6510e0f,    1.6530e0f,    1.6690e0f,    1.6900e0f,    1.7110e0f,
        1.7310e0f,    1.7490e0f,    1.7660e0f,    1.7820e0f,    1.0430e2f,    1.2390e2f,    1.4040e2f,    1.5500e2f,    1.6830e2f,    1.9210e2f,
        2.1310e2f,    2.3230e2f,    2.4990e2f,    2.6640e2f,    2.8190e2f,    2.9660e2f,    3.2450e2f,    3.4830e2f,    3.6890e2f,    3.8670e2f,
        4.0220e2f,    4.1570e2f,    4.2730e2f,    4.3730e2f,    4.5290e2f,    4.6380e2f,    4.7090e2f,    4.7490e2f,    4.7660e2f,    4.7640e2f,
        4.7490e2f,    4.7240e2f,    4.6910e2f,    4.6530e2f,    4.6110e2f,    4.5670e2f,    4.5220e2f,    4.4770e2f,    4.2530e2f,    4.0510e2f,
        3.8730e2f,    3.7150e2f,    3.5730e2f,    3.4440e2f,    3.3270e2f,    3.2180e2f,    3.0200e2f,    2.8440e2f,    2.6890e2f,    2.5500e2f,
        2.4270e2f,    2.3160e2f,    2.2160e2f,    2.1260e2f,    2.0430e2f,    1.9680e2f,    1.8990e2f,    1.8350e2f,    1.7750e2f,    1.7200e2f,
        1.4950e2f,    1.3280e2f,    1.1990e2f,    1.0950e2f,    1.0100e2f,    9.3830e1f,    8.7750e1f,    8.2500e1f,    7.3880e1f,    6.7070e1f,
        6.1540e1f,    5.6950e1f,    5.3060e1f,    4.9730e1f,    4.6840e1f,    4.4300e1f,    4.2050e1f,    4.0040e1f,    3.8240e1f,    3.6600e1f,
        3.5120e1f,    3.3760e1f,    2.8420e1f,    2.4660e1f,    2.1860e1f,    1.9690e1f,    1.6520e1f,    1.5320e1f,    1.4310e1f,    1.2680e1f,
        1.1420e1f,    1.0410e1f,    9.5940e0f,    8.9110e0f,    8.3340e0f,    7.8380e0f,    7.4080e0f,    7.0310e0f,    6.6980e0f,    6.4010e0f,
        6.1350e0f,    5.8950e0f,    5.6780e0f,    4.8370e0f,    4.2620e0f,    3.8440e0f,    3.5260e0f,    3.2770e0f,    3.0760e0f,    2.9110e0f,
        2.7730e0f,    2.5550e0f,    2.3930e0f,    2.2670e0f,    2.1670e0f,    2.0860e0f,    2.0200e0f,    1.9650e0f,    1.9180e0f,    1.8790e0f,
        1.8450e0f,    1.8160e0f,    1.7910e0f,    1.7690e0f,    1.7500e0f,    1.6470e0f,    1.6180e0f,    1.6130e0f,    1.6190e0f,    1.6420e0f,
        1.6680e0f,    1.6920e0f,    1.7140e0f,    1.7340e0f,    1.7520e0f,    1.7680e0f,    2.1470e2f,    2.4630e2f,    2.7460e2f,    3.0040e2f,
        3.2420e2f,    3.6350e2f,    3.9940e2f,    4.3220e2f,    4.6240e2f,    4.8990e2f,    5.1530e2f,    5.3910e2f,    5.8700e2f,    6.2770e2f,
        6.6320e2f,    6.9480e2f,    7.2290e2f,    7.4790e2f,    7.7000e2f,    7.8980e2f,    8.2440e2f,    8.5370e2f,    8.7790e2f,    8.9720e2f,
        9.1180e2f,    9.2240e2f,    9.2960e2f,    9.3400e2f,    9.3580e2f,    9.3550e2f,    9.3340e2f,    9.2960e2f,    9.2450e2f,    9.1830e2f,
        8.7570e2f,    8.2430e2f,    7.7230e2f,    7.2320e2f,    6.7690e2f,    6.3460e2f,    5.9690e2f,    5.6340e2f,    5.0770e2f,    4.6390e2f,
        4.2890e2f,    4.0060e2f,    3.7700e2f,    3.5640e2f,    3.3830e2f,    3.2220e2f,    3.0780e2f,    2.9480e2f,    2.8300e2f,    2.7220e2f,
        2.6230e2f,    2.5320e2f,    2.1680e2f,    1.9050e2f,    1.7050e2f,    1.5460e2f,    1.4180e2f,    1.3110e2f,    1.2210e2f,    1.1430e2f,
        1.0170e2f,    9.1790e1f,    8.3790e1f,    7.7190e1f,    7.1640e1f,    6.6900e1f,    6.2800e1f,    5.9210e1f,    5.6050e1f,    5.3240e1f,
        5.0730e1f,    4.8450e1f,    4.6400e1f,    4.4520e1f,    3.7190e1f,    3.2080e1f,    2.8310e1f,    2.5390e1f,    2.1180e1f,    1.9610e1f,
        1.8270e1f,    1.6130e1f,    1.4490e1f,    1.3180e1f,    1.2120e1f,    1.1240e1f,    1.0500e1f,    9.8580e0f,    9.3060e0f,    8.8230e0f,
        8.3970e0f,    8.0180e0f,    7.6780e0f,    7.3720e0f,    7.0950e0f,    6.0270e0f,    5.3000e0f,    4.7720e0f,    4.3720e0f,    4.0580e0f,
        3.8060e0f,    3.5990e0f,    3.4260e0f,    3.1540e0f,    2.9510e0f,    2.7940e0f,    2.6700e0f,    2.5690e0f,    2.4870e0f,    2.4180e0f,
        2.3610e0f,    2.3120e0f,    2.2690e0f,    2.2320e0f,    2.1990e0f,    2.1700e0f,    2.1450e0f,    2.0040e0f,    1.9540e0f,    1.9370e0f,
        1.9340e0f,    1.9450e0f,    1.9640e0f,    1.9850e0f,    2.0050e0f,    2.0240e0f,    2.0420e0f,    2.0590e0f,
        2.095e+02f,             2.404e+02f,             2.680e+02f,             2.932e+02f,             3.165e+02f,             3.561e+02f,             3.920e+02f,             4.247e+02f,             4.551e+02f,             4.830e+02f,
        5.088e+02f,             5.333e+02f,             5.817e+02f,             6.231e+02f,             6.595e+02f,             6.917e+02f,             7.206e+02f,             7.462e+02f,             7.691e+02f,             7.896e+02f,
        8.252e+02f,             8.548e+02f,             8.790e+02f,             8.981e+02f,             9.126e+02f,             9.231e+02f,             9.303e+02f,             9.345e+02f,             9.363e+02f,             9.359e+02f,
        9.336e+02f,             9.297e+02f,             9.246e+02f,             9.183e+02f,             8.753e+02f,             8.237e+02f,             7.715e+02f,             7.224e+02f,             6.773e+02f,             6.369e+02f,
        6.011e+02f,             5.694e+02f,             5.172e+02f,             4.754e+02f,             4.412e+02f,             4.123e+02f,             3.875e+02f,             3.659e+02f,             3.469e+02f,             3.301e+02f,
        3.150e+02f,             3.015e+02f,             2.892e+02f,             2.780e+02f,             2.677e+02f,             2.584e+02f,             2.207e+02f,             1.937e+02f,             1.731e+02f,             1.570e+02f,
        1.438e+02f,             1.330e+02f,             1.237e+02f,             1.158e+02f,             1.030e+02f,             9.289e+01f,             8.478e+01f,             7.808e+01f,             7.245e+01f,             6.764e+01f,
        6.348e+01f,             5.985e+01f,             5.664e+01f,             5.379e+01f,             5.125e+01f,             4.895e+01f,             4.687e+01f,             4.496e+01f,             3.755e+01f,             3.238e+01f,
        2.856e+01f,             2.562e+01f,             2.136e+01f,             1.977e+01f,             1.842e+01f,             1.626e+01f,             1.460e+01f,             1.328e+01f,             1.221e+01f,             1.132e+01f,
        1.057e+01f,             9.929e+00f,             9.372e+00f,             8.886e+00f,             8.455e+00f,             8.073e+00f,             7.730e+00f,             7.422e+00f,             7.142e+00f,             6.065e+00f,
        5.331e+00f,             4.800e+00f,             4.396e+00f,             4.080e+00f,             3.825e+00f,             3.616e+00f,             3.441e+00f,             3.167e+00f,             2.961e+00f,             2.802e+00f,
        2.677e+00f,             2.575e+00f,             2.491e+00f,             2.421e+00f,             2.362e+00f,             2.313e+00f,             2.270e+00f,             2.233e+00f,             2.200e+00f,             2.173e+00f,
        2.148e+00f,             2.016e+00f,             1.976e+00f,             1.968e+00f,             1.973e+00f,             2.000e+00f,             2.032e+00f,             2.064e+00f,             2.094e+00f,             2.122e+00f,
        2.149e+00f,             2.173e+00f,
        1.008e+02f,             1.157e+02f,             1.289e+02f,             1.408e+02f,             1.518e+02f,             1.717e+02f,             1.895e+02f,             2.058e+02f,             2.210e+02f,             2.352e+02f,
         2.485e+02f,             2.613e+02f,             2.866e+02f,             3.090e+02f,             3.291e+02f,             3.474e+02f,             3.640e+02f,             3.793e+02f,             3.934e+02f,             4.065e+02f,
         4.298e+02f,             4.498e+02f,             4.671e+02f,             4.820e+02f,             4.947e+02f,             5.055e+02f,             5.147e+02f,             5.223e+02f,             5.286e+02f,             5.337e+02f,
         5.376e+02f,             5.406e+02f,             5.427e+02f,             5.440e+02f,             5.411e+02f,             5.281e+02f,             5.098e+02f,             4.893e+02f,             4.679e+02f,             4.469e+02f,
         4.273e+02f,             4.091e+02f,             3.773e+02f,             3.507e+02f,             3.283e+02f,             3.093e+02f,             2.930e+02f,             2.783e+02f,             2.652e+02f,             2.534e+02f,
         2.427e+02f,             2.329e+02f,             2.240e+02f,             2.158e+02f,             2.083e+02f,             2.014e+02f,             1.732e+02f,             1.526e+02f,             1.369e+02f,             1.244e+02f,
         1.142e+02f,             1.058e+02f,             9.859e+01f,             9.243e+01f,             8.237e+01f,             7.447e+01f,             6.807e+01f,             6.279e+01f,             5.834e+01f,             5.453e+01f,
         5.124e+01f,             4.836e+01f,             4.581e+01f,             4.354e+01f,             4.151e+01f,             3.967e+01f,             3.801e+01f,             3.649e+01f,             3.054e+01f,             2.639e+01f,
         2.331e+01f,             2.093e+01f,             1.749e+01f,             1.620e+01f,             1.510e+01f,             1.334e+01f,             1.199e+01f,             1.092e+01f,             1.005e+01f,             9.320e+00f,
         8.707e+00f,             8.181e+00f,             7.725e+00f,             7.327e+00f,             6.975e+00f,             6.661e+00f,             6.381e+00f,             6.128e+00f,             5.899e+00f,             5.015e+00f,
         4.412e+00f,             3.975e+00f,             3.643e+00f,             3.383e+00f,             3.174e+00f,             3.002e+00f,             2.858e+00f,             2.633e+00f,             2.464e+00f,             2.333e+00f,
         2.229e+00f,             2.144e+00f,             2.075e+00f,             2.017e+00f,             1.968e+00f,             1.927e+00f,             1.891e+00f,             1.860e+00f,             1.833e+00f,             1.810e+00f,
         1.789e+00f,             1.675e+00f,             1.635e+00f,             1.622e+00f,             1.621e+00f,             1.631e+00f,             1.648e+00f,             1.665e+00f,             1.682e+00f,             1.697e+00f,
         1.712e+00f,             1.726e+00f,
		1.414e+02f,		1.651e+02f,		1.855e+02f,		2.038e+02f,		2.206e+02f,		2.507e+02f,		2.776e+02f,		3.021e+02f,		3.248e+02f,		3.460e+02f,
		3.660e+02f,		3.850e+02f,		4.224e+02f,		4.552e+02f,		4.843e+02f,		5.106e+02f,		5.343e+02f,		5.558e+02f,		5.755e+02f,		5.934e+02f,
		6.246e+02f,		6.506e+02f,		6.721e+02f,		6.897e+02f,		7.038e+02f,		7.149e+02f,		7.233e+02f,		7.293e+02f,		7.333e+02f,		7.355e+02f,
		7.360e+02f,		7.352e+02f,		7.332e+02f,		7.301e+02f,		7.038e+02f,		6.680e+02f,		6.298e+02f,		5.928e+02f,		5.589e+02f,		5.284e+02f,
		5.011e+02f,		4.767e+02f,		4.353e+02f,		4.015e+02f,		3.736e+02f,		3.501e+02f,		3.300e+02f,		3.123e+02f,		2.967e+02f,		2.826e+02f,
		2.701e+02f,		2.589e+02f,		2.486e+02f,		2.393e+02f,		2.308e+02f,		2.229e+02f,		1.912e+02f,		1.683e+02f,		1.509e+02f,		1.371e+02f,
		1.258e+02f,		1.165e+02f,		1.086e+02f,		1.018e+02f,		9.068e+01f,		8.197e+01f,		7.492e+01f,		6.909e+01f,		6.417e+01f,		5.997e+01f,
		5.633e+01f,		5.315e+01f,		5.033e+01f,		4.783e+01f,		4.559e+01f,		4.357e+01f,		4.173e+01f,		4.006e+01f,		3.351e+01f,		2.894e+01f,
		2.555e+01f,		2.294e+01f,		1.915e+01f,		1.773e+01f,		1.653e+01f,		1.460e+01f,		1.312e+01f,		1.194e+01f,		1.099e+01f,		1.019e+01f,
		9.517e+00f,		8.942e+00f,		8.443e+00f,		8.006e+00f,		7.620e+00f,		7.277e+00f,		6.970e+00f,		6.693e+00f,		6.443e+00f,		5.475e+00f,
		4.816e+00f,		4.338e+00f,		3.976e+00f,		3.691e+00f,		3.462e+00f,		3.275e+00f,		3.118e+00f,		2.871e+00f,		2.687e+00f,		2.544e+00f,
		2.431e+00f,		2.340e+00f,		2.266e+00f,		2.203e+00f,		2.151e+00f,		2.107e+00f,		2.069e+00f,		2.037e+00f,		2.008e+00f,		1.984e+00f,
		1.963e+00f,		1.850e+00f,		1.820e+00f,		1.818e+00f,		1.828e+00f,		1.861e+00f,		1.898e+00f,		1.934e+00f,		1.967e+00f,		1.998e+00f,
		2.026e+00f,		2.052e+00f
   },
    {
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    4.0000e-4f,
        4.0000e-4f,    5.0000e-4f,    6.0000e-4f,    7.0000e-4f,    8.0000e-4f,    9.0000e-4f,    1.0000e-3f,    1.1000e-3f,    1.3000e-3f,    1.4000e-3f,
        1.6000e-3f,    1.7000e-3f,    1.9000e-3f,    2.1000e-3f,    2.3000e-3f,    2.5000e-3f,    3.5000e-3f,    4.7000e-3f,    6.1000e-3f,    7.6000e-3f,
        9.2000e-3f,    1.1000e-2f,    1.2900e-2f,    1.5000e-2f,    1.9500e-2f,    2.4600e-2f,    3.0200e-2f,    3.6200e-2f,    4.2800e-2f,    4.9800e-2f,
        5.7400e-2f,    6.5400e-2f,    7.3800e-2f,    8.2800e-2f,    9.2200e-2f,    1.0200e-1f,    1.1230e-1f,    1.2300e-1f,    1.8320e-1f,    2.5390e-1f,
        3.3500e-1f,    4.2600e-1f,    6.3700e-1f,    7.5660e-1f,    8.8530e-1f,    1.1700e0f,    1.4890e0f,    1.8410e0f,    2.2270e0f,    2.6440e0f,
        3.0930e0f,    3.5720e0f,    4.0800e0f,    4.6180e0f,    5.1840e0f,    5.7770e0f,    6.3980e0f,    7.0450e0f,    7.7180e0f,    1.1460e1f,
        1.5770e1f,    2.0620e1f,    2.5960e1f,    3.1740e1f,    3.7940e1f,    4.4520e1f,    5.1450e1f,    6.6280e1f,    8.2250e1f,    9.9210e1f,
        1.1700e2f,    1.3560e2f,    1.5490e2f,    1.7470e2f,    1.9510e2f,    2.1590e2f,    2.3720e2f,    2.5880e2f,    2.8070e2f,    3.0290e2f,
        3.2540e2f,    5.6050e2f,    8.0540e2f,    1.0540e3f,    1.3040e3f,    1.8020e3f,    2.2970e3f,    2.7870e3f,    3.2720e3f,    3.7520e3f,
        4.2280e3f,    4.7000e3f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,
        2.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    4.0000e-4f,    5.0000e-4f,    5.0000e-4f,
        6.0000e-4f,    6.0000e-4f,    7.0000e-4f,    8.0000e-4f,    9.0000e-4f,    1.1000e-3f,    1.2000e-3f,    1.4000e-3f,    1.6000e-3f,    1.8000e-3f,
        2.0000e-3f,    2.2000e-3f,    2.4000e-3f,    2.7000e-3f,    2.9000e-3f,    3.2000e-3f,    3.4000e-3f,    3.7000e-3f,    5.2000e-3f,    6.8000e-3f,
        8.7000e-3f,    1.0800e-2f,    1.3000e-2f,    1.5400e-2f,    1.8000e-2f,    2.0800e-2f,    2.6900e-2f,    3.3600e-2f,    4.1000e-2f,    4.9000e-2f,
        5.7700e-2f,    6.6900e-2f,    7.6800e-2f,    8.7200e-2f,    9.8200e-2f,    1.0980e-1f,    1.2200e-1f,    1.3470e-1f,    1.4800e-1f,    1.6180e-1f,
        2.3910e-1f,    3.2960e-1f,    4.3280e-1f,    5.4830e-1f,    8.1510e-1f,    9.6590e-1f,    1.1280e0f,    1.4850e0f,    1.8850e0f,    2.3260e0f,
        2.8080e0f,    3.3280e0f,    3.8870e0f,    4.4820e0f,    5.1140e0f,    5.7810e0f,    6.4830e0f,    7.2190e0f,    7.9870e0f,    8.7890e0f,
        9.6210e0f,    1.4240e1f,    1.9560e1f,    2.5530e1f,    3.2090e1f,    3.9190e1f,    4.6790e1f,    5.4850e1f,    6.3340e1f,    8.1480e1f,
        1.0100e2f,    1.2170e2f,    1.4350e2f,    1.6610e2f,    1.8960e2f,    2.1390e2f,    2.3870e2f,    2.6420e2f,    2.9010e2f,    3.1650e2f,
        3.4330e2f,    3.7040e2f,    3.9790e2f,    6.8500e2f,    9.8350e2f,    1.2860e3f,    1.5890e3f,    2.1910e3f,    2.7860e3f,    3.3740e3f,
        3.9550e3f,    4.5300e3f,    5.0990e3f,    5.6630e3f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,
        2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    4.0000e-4f,
        4.0000e-4f,    5.0000e-4f,    6.0000e-4f,    6.0000e-4f,    7.0000e-4f,    8.0000e-4f,    1.0000e-3f,    1.1000e-3f,    1.3000e-3f,    1.5000e-3f,
        1.7000e-3f,    1.9000e-3f,    2.1000e-3f,    2.4000e-3f,    2.6000e-3f,    2.9000e-3f,    3.1000e-3f,    3.4000e-3f,    3.7000e-3f,    3.9000e-3f,
        5.5000e-3f,    7.3000e-3f,    9.3000e-3f,    1.1500e-2f,    1.3800e-2f,    1.6400e-2f,    1.9200e-2f,    2.2100e-2f,    2.8500e-2f,    3.5600e-2f,
        4.3400e-2f,    5.1900e-2f,    6.1000e-2f,    7.0700e-2f,    8.1100e-2f,    9.2100e-2f,    1.0370e-1f,    1.1580e-1f,    1.2860e-1f,    1.4200e-1f,
        1.5590e-1f,    1.7050e-1f,    2.5150e-1f,    3.4620e-1f,    4.5410e-1f,    5.7480e-1f,    8.5330e-1f,    1.0110e0f,    1.1800e0f,    1.5520e0f,
        1.9680e0f,    2.4270e0f,    2.9280e0f,    3.4690e0f,    4.0500e0f,    4.6690e0f,    5.3250e0f,    6.0180e0f,    6.7470e0f,    7.5110e0f,
        8.3090e0f,    9.1410e0f,    1.0010e1f,    1.4800e1f,    2.0320e1f,    2.6510e1f,    3.3310e1f,    4.0670e1f,    4.8550e1f,    5.6910e1f,
        6.5720e1f,    8.4540e1f,    1.0480e2f,    1.2630e2f,    1.4890e2f,    1.7240e2f,    1.9680e2f,    2.2190e2f,    2.4760e2f,    2.7400e2f,
        3.0080e2f,    3.2820e2f,    3.5590e2f,    3.8400e2f,    4.1240e2f,    7.0870e2f,    1.0160e3f,    1.3250e3f,    1.6350e3f,    2.2480e3f,
        2.8530e3f,    3.4480e3f,    4.0350e3f,    4.6150e3f,    5.1890e3f,    5.7570e3f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    3.0000e-4f,    3.0000e-4f,    3.0000e-4f,    4.0000e-4f,    4.0000e-4f,    5.0000e-4f,    6.0000e-4f,
        7.0000e-4f,    8.0000e-4f,    1.0000e-3f,    1.1000e-3f,    1.3000e-3f,    1.4000e-3f,    1.6000e-3f,    1.7000e-3f,    1.9000e-3f,    2.1000e-3f,
        2.3000e-3f,    2.5000e-3f,    3.5000e-3f,    4.8000e-3f,    6.2000e-3f,    7.7000e-3f,    9.4000e-3f,    1.1200e-2f,    1.3200e-2f,    1.5300e-2f,
        2.0000e-2f,    2.5200e-2f,    3.0900e-2f,    3.7100e-2f,    4.3800e-2f,    5.1000e-2f,    5.8800e-2f,    6.7000e-2f,    7.5600e-2f,    8.4800e-2f,
        9.4400e-2f,    1.0450e-1f,    1.1510e-1f,    1.2610e-1f,    1.8780e-1f,    2.6040e-1f,    3.4350e-1f,    4.3690e-1f,    6.5350e-1f,    7.7630e-1f,
        9.0850e-1f,    1.2000e0f,    1.5280e0f,    1.8900e0f,    2.2860e0f,    2.7150e0f,    3.1760e0f,    3.6680e0f,    4.1900e0f,    4.7420e0f,
        5.3230e0f,    5.9330e0f,    6.5700e0f,    7.2350e0f,    7.9270e0f,    1.1770e1f,    1.6200e1f,    2.1190e1f,    2.6670e1f,    3.2610e1f,
        3.8980e1f,    4.5740e1f,    5.2860e1f,    6.8100e1f,    8.4510e1f,    1.0190e2f,    1.2030e2f,    1.3940e2f,    1.5920e2f,    1.7960e2f,
        2.0050e2f,    2.2190e2f,    2.4370e2f,    2.6600e2f,    2.8850e2f,    3.1140e2f,    3.3460e2f,    5.7720e2f,    8.3030e2f,    1.0880e3f,
        1.3460e3f,    1.8620e3f,    2.3740e3f,    2.8800e3f,    3.3810e3f,    3.8770e3f,    4.3690e3f,    4.8570e3f,
        9.096e-06f,             1.213e-05f,             1.481e-05f,             1.722e-05f,             1.944e-05f,             2.345e-05f,             2.706e-05f,             3.037e-05f,             3.343e-05f,             3.631e-05f,
        3.903e-05f,             4.163e-05f,             4.769e-05f,             5.330e-05f,             5.856e-05f,             6.354e-05f,             6.831e-05f,             7.288e-05f,             7.730e-05f,             8.161e-05f,
        8.983e-05f,             9.772e-05f,             1.053e-04f,             1.127e-04f,             1.199e-04f,             1.269e-04f,             1.339e-04f,             1.407e-04f,             1.475e-04f,             1.543e-04f,
        1.610e-04f,             1.677e-04f,             1.744e-04f,             1.812e-04f,             2.153e-04f,             2.509e-04f,             2.884e-04f,             3.280e-04f,             3.699e-04f,             4.142e-04f,
        4.611e-04f,             5.103e-04f,             6.157e-04f,             7.303e-04f,             8.534e-04f,             9.853e-04f,             1.125e-03f,             1.273e-03f,             1.429e-03f,             1.594e-03f,
        1.766e-03f,             1.946e-03f,             2.132e-03f,             2.327e-03f,             2.529e-03f,             2.739e-03f,             3.892e-03f,             5.218e-03f,             6.710e-03f,             8.362e-03f,
        1.017e-02f,             1.213e-02f,             1.425e-02f,             1.651e-02f,             2.146e-02f,             2.698e-02f,             3.305e-02f,             3.966e-02f,             4.681e-02f,             5.448e-02f,
        6.267e-02f,             7.137e-02f,             8.057e-02f,             9.027e-02f,             1.005e-01f,             1.112e-01f,             1.223e-01f,             1.340e-01f,             1.992e-01f,             2.759e-01f,
        3.637e-01f,             4.623e-01f,             6.907e-01f,             8.201e-01f,             9.595e-01f,             1.267e+00f,             1.612e+00f,             1.993e+00f,             2.410e+00f,             2.861e+00f,
        3.346e+00f,             3.863e+00f,             4.413e+00f,             4.993e+00f,             5.605e+00f,             6.245e+00f,             6.916e+00f,             7.614e+00f,             8.341e+00f,             1.238e+01f,
        1.704e+01f,             2.228e+01f,             2.803e+01f,             3.428e+01f,             4.097e+01f,             4.807e+01f,             5.556e+01f,             7.158e+01f,             8.884e+01f,             1.072e+02f,
        1.265e+02f,             1.466e+02f,             1.674e+02f,             1.889e+02f,             2.109e+02f,             2.335e+02f,             2.565e+02f,             2.799e+02f,             3.037e+02f,             3.278e+02f,
        3.522e+02f,             6.070e+02f,             8.716e+02f,             1.139e+03f,             1.407e+03f,             1.938e+03f,             2.461e+03f,             2.976e+03f,             3.484e+03f,             3.984e+03f,
        4.479e+03f,             4.967e+03f,
        1.271e-05f,             1.730e-05f,             2.139e-05f,             2.510e-05f,             2.852e-05f,             3.470e-05f,             4.024e-05f,             4.530e-05f,             4.999e-05f,             5.437e-05f,
        5.851e-05f,             6.243e-05f,             7.154e-05f,             7.994e-05f,             8.777e-05f,             9.516e-05f,             1.022e-04f,             1.089e-04f,             1.154e-04f,             1.216e-04f,
        1.336e-04f,             1.450e-04f,             1.559e-04f,             1.664e-04f,             1.766e-04f,             1.866e-04f,             1.964e-04f,             2.061e-04f,             2.156e-04f,             2.250e-04f,
        2.343e-04f,             2.436e-04f,             2.528e-04f,             2.620e-04f,             3.080e-04f,             3.547e-04f,             4.029e-04f,             4.529e-04f,             5.052e-04f,             5.598e-04f,
        6.171e-04f,             6.769e-04f,             8.043e-04f,             9.418e-04f,             1.089e-03f,             1.246e-03f,             1.412e-03f,             1.588e-03f,             1.772e-03f,             1.965e-03f,
        2.166e-03f,             2.377e-03f,             2.596e-03f,             2.823e-03f,             3.059e-03f,             3.303e-03f,             4.646e-03f,             6.187e-03f,             7.920e-03f,             9.838e-03f,
        1.194e-02f,             1.421e-02f,             1.666e-02f,             1.928e-02f,             2.503e-02f,             3.142e-02f,             3.845e-02f,             4.611e-02f,             5.437e-02f,             6.324e-02f,
        7.271e-02f,             8.276e-02f,             9.338e-02f,             1.046e-01f,             1.163e-01f,             1.287e-01f,             1.416e-01f,             1.550e-01f,             2.302e-01f,             3.185e-01f,
        4.196e-01f,             5.329e-01f,             7.955e-01f,             9.441e-01f,             1.104e+00f,             1.457e+00f,             1.853e+00f,             2.290e+00f,             2.768e+00f,             3.285e+00f,
        3.841e+00f,             4.434e+00f,             5.063e+00f,             5.728e+00f,             6.428e+00f,             7.161e+00f,             7.929e+00f,             8.728e+00f,             9.560e+00f,             1.418e+01f,
        1.951e+01f,             2.549e+01f,             3.207e+01f,             3.920e+01f,             4.684e+01f,             5.494e+01f,             6.348e+01f,             8.174e+01f,             1.014e+02f,             1.223e+02f,
        1.442e+02f,             1.671e+02f,             1.908e+02f,             2.153e+02f,             2.404e+02f,             2.660e+02f,             2.922e+02f,             3.189e+02f,             3.460e+02f,             3.734e+02f,
        4.012e+02f,             6.918e+02f,             9.945e+02f,             1.302e+03f,             1.610e+03f,             2.226e+03f,             2.836e+03f,             3.439e+03f,             4.037e+03f,             4.629e+03f,
        5.215e+03f,             5.797e+03f,
		9.857e-06f,		1.310e-05f,		1.595e-05f,		1.852e-05f,		2.088e-05f,		2.512e-05f,		2.891e-05f,		3.236e-05f,		3.555e-05f,		3.853e-05f,
		4.134e-05f,		4.400e-05f,		5.019e-05f,		5.588e-05f,		6.120e-05f,		6.623e-05f,		7.101e-05f,		7.560e-05f,		8.002e-05f,		8.430e-05f,
		9.250e-05f,		1.003e-04f,		1.079e-04f,		1.152e-04f,		1.224e-04f,		1.295e-04f,		1.364e-04f,		1.433e-04f,		1.501e-04f,		1.569e-04f,
		1.637e-04f,		1.705e-04f,		1.773e-04f,		1.842e-04f,		2.190e-04f,		2.554e-04f,		2.940e-04f,		3.349e-04f,		3.783e-04f,		4.244e-04f,
		4.730e-04f,		5.241e-04f,		6.340e-04f,		7.538e-04f,		8.830e-04f,		1.021e-03f,		1.169e-03f,		1.324e-03f,		1.489e-03f,		1.661e-03f,
		1.842e-03f,		2.032e-03f,		2.229e-03f,		2.434e-03f,		2.646e-03f,		2.867e-03f,		4.082e-03f,		5.479e-03f,		7.051e-03f,		8.792e-03f,
		1.070e-02f,		1.276e-02f,		1.499e-02f,		1.737e-02f,		2.258e-02f,		2.839e-02f,		3.478e-02f,		4.173e-02f,		4.925e-02f,		5.731e-02f,
		6.592e-02f,		7.506e-02f,		8.474e-02f,		9.493e-02f,		1.056e-01f,		1.169e-01f,		1.286e-01f,		1.408e-01f,		2.094e-01f,		2.899e-01f,
		3.820e-01f,		4.855e-01f,		7.252e-01f,		8.609e-01f,		1.007e+00f,		1.330e+00f,		1.691e+00f,		2.091e+00f,		2.528e+00f,		3.001e+00f,
		3.509e+00f,		4.052e+00f,		4.628e+00f,		5.236e+00f,		5.876e+00f,		6.548e+00f,		7.250e+00f,		7.983e+00f,		8.744e+00f,		1.297e+01f,
		1.786e+01f,		2.334e+01f,		2.937e+01f,		3.590e+01f,		4.290e+01f,		5.033e+01f,		5.816e+01f,		7.490e+01f,		9.293e+01f,		1.121e+02f,
		1.322e+02f,		1.532e+02f,		1.749e+02f,		1.973e+02f,		2.203e+02f,		2.437e+02f,		2.677e+02f,		2.921e+02f,		3.168e+02f,		3.418e+02f,
		3.672e+02f,		6.311e+02f,		9.041e+02f,		1.179e+03f,		1.454e+03f,		1.996e+03f,		2.528e+03f,		3.050e+03f,		3.563e+03f,		4.067e+03f,
		4.564e+03f,		5.054e+03f
    },
    {
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    3.0000e-4f,    3.0000e-4f,    4.0000e-4f,
        4.0000e-4f,    4.0000e-4f,    5.0000e-4f,    6.0000e-4f,    8.0000e-4f,    9.0000e-4f,    1.0000e-3f,    1.1000e-3f,    1.3000e-3f,    1.4000e-3f,
        1.6000e-3f,    1.7000e-3f,    1.9000e-3f,    2.1000e-3f,    2.2000e-3f,    2.4000e-3f,    3.5000e-3f,    4.7000e-3f,    6.0000e-3f,    7.5000e-3f,
        9.2000e-3f,    1.1000e-2f,    1.2900e-2f,    1.4900e-2f,    1.9500e-2f,    2.4500e-2f,    3.0100e-2f,    3.6100e-2f,    4.2700e-2f,    4.9700e-2f,
        5.7200e-2f,    6.5200e-2f,    7.3700e-2f,    8.2600e-2f,    9.2000e-2f,    1.0180e-1f,    1.1200e-1f,    1.2280e-1f,    1.8280e-1f,    2.5350e-1f,
        3.3440e-1f,    4.2520e-1f,    6.3590e-1f,    7.5530e-1f,    8.8390e-1f,    1.1680e0f,    1.4860e0f,    1.8390e0f,    2.2240e0f,    2.6410e0f,
        3.0890e0f,    3.5670e0f,    4.0750e0f,    4.6110e0f,    5.1760e0f,    5.7690e0f,    6.3890e0f,    7.0350e0f,    7.7070e0f,    1.1440e1f,
        1.5760e1f,    2.0600e1f,    2.5930e1f,    3.1710e1f,    3.7900e1f,    4.4470e1f,    5.1390e1f,    6.6210e1f,    8.2170e1f,    9.9120e1f,
        1.1690e2f,    1.3550e2f,    1.5470e2f,    1.7460e2f,    1.9490e2f,    2.1580e2f,    2.3700e2f,    2.5860e2f,    2.8050e2f,    3.0270e2f,
        3.2520e2f,    5.6010e2f,    8.0490e2f,    1.0530e3f,    1.3030e3f,    1.8020e3f,    2.2960e3f,    2.7860e3f,    3.2710e3f,    3.7510e3f,
        4.2270e3f,    4.6990e3f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,
        2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    3.0000e-4f,    4.0000e-4f,    4.0000e-4f,    5.0000e-4f,
        5.0000e-4f,    6.0000e-4f,    7.0000e-4f,    7.0000e-4f,    9.0000e-4f,    1.0000e-3f,    1.2000e-3f,    1.4000e-3f,    1.5000e-3f,    1.7000e-3f,
        1.9000e-3f,    2.2000e-3f,    2.4000e-3f,    2.6000e-3f,    2.9000e-3f,    3.1000e-3f,    3.4000e-3f,    3.6000e-3f,    5.1000e-3f,    6.8000e-3f,
        8.6000e-3f,    1.0700e-2f,    1.2900e-2f,    1.5300e-2f,    1.7900e-2f,    2.0700e-2f,    2.6700e-2f,    3.3400e-2f,    4.0800e-2f,    4.8800e-2f,
        5.7400e-2f,    6.6600e-2f,    7.6400e-2f,    8.6800e-2f,    9.7800e-2f,    1.0930e-1f,    1.2150e-1f,    1.3420e-1f,    1.4740e-1f,    1.6120e-1f,
        2.3830e-1f,    3.2850e-1f,    4.3130e-1f,    5.4650e-1f,    8.1270e-1f,    9.6310e-1f,    1.1250e0f,    1.4810e0f,    1.8800e0f,    2.3200e0f,
        2.8000e0f,    3.3190e0f,    3.8770e0f,    4.4710e0f,    5.1010e0f,    5.7670e0f,    6.4670e0f,    7.2010e0f,    7.9680e0f,    8.7680e0f,
        9.5990e0f,    1.4210e1f,    1.9520e1f,    2.5480e1f,    3.2020e1f,    3.9110e1f,    4.6700e1f,    5.4740e1f,    6.3220e1f,    8.1330e1f,
        1.0080e2f,    1.2150e2f,    1.4320e2f,    1.6590e2f,    1.8930e2f,    2.1350e2f,    2.3840e2f,    2.6380e2f,    2.8970e2f,    3.1600e2f,
        3.4280e2f,    3.6990e2f,    3.9740e2f,    6.8430e2f,    9.8250e2f,    1.2850e3f,    1.5870e3f,    2.1890e3f,    2.7850e3f,    3.3720e3f,
        3.9530e3f,    4.5280e3f,    5.0970e3f,    5.6600e3f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    3.0000e-4f,    3.0000e-4f,
        4.0000e-4f,    5.0000e-4f,    5.0000e-4f,    6.0000e-4f,    7.0000e-4f,    8.0000e-4f,    9.0000e-4f,    1.1000e-3f,    1.3000e-3f,    1.4000e-3f,
        1.6000e-3f,    1.9000e-3f,    2.1000e-3f,    2.3000e-3f,    2.5000e-3f,    2.8000e-3f,    3.0000e-3f,    3.3000e-3f,    3.6000e-3f,    3.9000e-3f,
        5.4000e-3f,    7.2000e-3f,    9.2000e-3f,    1.1300e-2f,    1.3700e-2f,    1.6300e-2f,    1.9000e-2f,    2.1900e-2f,    2.8300e-2f,    3.5400e-2f,
        4.3200e-2f,    5.1600e-2f,    6.0600e-2f,    7.0300e-2f,    8.0600e-2f,    9.1600e-2f,    1.0310e-1f,    1.1530e-1f,    1.2800e-1f,    1.4130e-1f,
        1.5520e-1f,    1.6970e-1f,    2.5040e-1f,    3.4480e-1f,    4.5230e-1f,    5.7260e-1f,    8.5010e-1f,    1.0070e0f,    1.1750e0f,    1.5460e0f,
        1.9610e0f,    2.4190e0f,    2.9180e0f,    3.4580e0f,    4.0370e0f,    4.6540e0f,    5.3090e0f,    6.0000e0f,    6.7270e0f,    7.4880e0f,
        8.2840e0f,    9.1140e0f,    9.9760e0f,    1.4760e1f,    2.0260e1f,    2.6440e1f,    3.3220e1f,    4.0570e1f,    4.8430e1f,    5.6780e1f,
        6.5560e1f,    8.4340e1f,    1.0460e2f,    1.2600e2f,    1.4850e2f,    1.7200e2f,    1.9640e2f,    2.2140e2f,    2.4720e2f,    2.7350e2f,
        3.0030e2f,    3.2760e2f,    3.5530e2f,    3.8330e2f,    4.1170e2f,    7.0770e2f,    1.0140e3f,    1.3240e3f,    1.6330e3f,    2.2460e3f,
        2.8500e3f,    3.4450e3f,    4.0320e3f,    4.6120e3f,    5.1860e3f,    5.7540e3f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,
        0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    0.0000e0f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,    1.0000e-4f,
        2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    2.0000e-4f,    3.0000e-4f,    3.0000e-4f,    4.0000e-4f,    4.0000e-4f,    5.0000e-4f,    6.0000e-4f,
        7.0000e-4f,    8.0000e-4f,    1.0000e-3f,    1.1000e-3f,    1.2000e-3f,    1.4000e-3f,    1.6000e-3f,    1.7000e-3f,    1.9000e-3f,    2.1000e-3f,
        2.3000e-3f,    2.4000e-3f,    3.5000e-3f,    4.7000e-3f,    6.1000e-3f,    7.7000e-3f,    9.4000e-3f,    1.1200e-2f,    1.3200e-2f,    1.5300e-2f,
        1.9900e-2f,    2.5100e-2f,    3.0800e-2f,    3.7000e-2f,    4.3700e-2f,    5.0900e-2f,    5.8600e-2f,    6.6800e-2f,    7.5500e-2f,    8.4600e-2f,
        9.4200e-2f,    1.0430e-1f,    1.1490e-1f,    1.2580e-1f,    1.8750e-1f,    2.6000e-1f,    3.4300e-1f,    4.3630e-1f,    6.5260e-1f,    7.7520e-1f,
        9.0720e-1f,    1.1990e0f,    1.5260e0f,    1.8880e0f,    2.2830e0f,    2.7120e0f,    3.1720e0f,    3.6630e0f,    4.1850e0f,    4.7360e0f,
        5.3170e0f,    5.9260e0f,    6.5620e0f,    7.2260e0f,    7.9170e0f,    1.1750e1f,    1.6190e1f,    2.1160e1f,    2.6640e1f,    3.2580e1f,
        3.8940e1f,    4.5690e1f,    5.2810e1f,    6.8040e1f,    8.4430e1f,    1.0190e2f,    1.2020e2f,    1.3930e2f,    1.5900e2f,    1.7940e2f,
        2.0030e2f,    2.2170e2f,    2.4350e2f,    2.6580e2f,    2.8830e2f,    3.1120e2f,    3.3440e2f,    5.7680e2f,    8.2990e2f,    1.0870e3f,
        1.3450e3f,    1.8610e3f,    2.3730e3f,    2.8790e3f,    3.3800e3f,    3.8760e3f,    4.3680e3f,    4.8560e3f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,
        3.947e-06f,             6.091e-06f,             8.255e-06f,             1.039e-05f,             1.250e-05f,             1.660e-05f,             2.052e-05f,             2.430e-05f,             2.792e-05f,             3.142e-05f,
        3.479e-05f,             3.805e-05f,             4.582e-05f,             5.318e-05f,             6.018e-05f,             6.689e-05f,             7.333e-05f,             7.956e-05f,             8.559e-05f,             9.146e-05f,
        1.028e-04f,             1.136e-04f,             1.240e-04f,             1.342e-04f,             1.441e-04f,             1.538e-04f,             1.633e-04f,             1.727e-04f,             1.820e-04f,             1.912e-04f,
        2.003e-04f,             2.094e-04f,             2.185e-04f,             2.275e-04f,             2.728e-04f,             3.189e-04f,             3.666e-04f,             4.162e-04f,             4.680e-04f,             5.222e-04f,
        5.790e-04f,             6.384e-04f,             7.651e-04f,             9.019e-04f,             1.049e-03f,             1.205e-03f,             1.370e-03f,             1.545e-03f,             1.728e-03f,             1.920e-03f,
        2.121e-03f,             2.330e-03f,             2.549e-03f,             2.775e-03f,             3.010e-03f,             3.254e-03f,             4.592e-03f,             6.129e-03f,             7.856e-03f,             9.769e-03f,
        1.186e-02f,             1.413e-02f,             1.658e-02f,             1.919e-02f,             2.492e-02f,             3.130e-02f,             3.831e-02f,             4.595e-02f,             5.420e-02f,             6.305e-02f,
        7.249e-02f,             8.252e-02f,             9.312e-02f,             1.043e-01f,             1.160e-01f,             1.283e-01f,             1.412e-01f,             1.546e-01f,             2.297e-01f,             3.178e-01f,
        4.187e-01f,             5.318e-01f,             7.939e-01f,             9.423e-01f,             1.102e+00f,             1.454e+00f,             1.850e+00f,             2.286e+00f,             2.763e+00f,             3.280e+00f,
        3.834e+00f,             4.426e+00f,             5.055e+00f,             5.718e+00f,             6.417e+00f,             7.150e+00f,             7.916e+00f,             8.715e+00f,             9.545e+00f,             1.416e+01f,
        1.948e+01f,             2.545e+01f,             3.203e+01f,             3.915e+01f,             4.677e+01f,             5.487e+01f,             6.340e+01f,             8.164e+01f,             1.013e+02f,             1.221e+02f,
        1.441e+02f,             1.669e+02f,             1.906e+02f,             2.150e+02f,             2.401e+02f,             2.658e+02f,             2.920e+02f,             3.186e+02f,             3.457e+02f,             3.731e+02f,
        4.009e+02f,             6.912e+02f,             9.939e+02f,             1.301e+03f,             1.609e+03f,             2.224e+03f,             2.834e+03f,             3.438e+03f,             4.035e+03f,             4.627e+03f,
        5.214e+03f,             5.795e+03f,
		3.257e-06f,		4.925e-06f,		6.577e-06f,		8.189e-06f,		9.759e-06f,		1.277e-05f,		1.563e-05f,		1.834e-05f,		2.094e-05f,		2.342e-05f,
		2.580e-05f,		2.810e-05f,		3.356e-05f,		3.869e-05f,		4.356e-05f,		4.822e-05f,		5.269e-05f,		5.701e-05f,		6.120e-05f,		6.527e-05f,
		7.313e-05f,		8.069e-05f,		8.800e-05f,		9.513e-05f,		1.021e-04f,		1.090e-04f,		1.158e-04f,		1.226e-04f,		1.293e-04f,		1.360e-04f,
		1.426e-04f,		1.493e-04f,		1.561e-04f,		1.628e-04f,		1.972e-04f,		2.333e-04f,		2.715e-04f,		3.121e-04f,		3.553e-04f,		4.010e-04f,
		4.493e-04f,		5.002e-04f,		6.095e-04f,		7.287e-04f,		8.573e-04f,		9.951e-04f,		1.142e-03f,		1.297e-03f,		1.461e-03f,		1.633e-03f,
		1.813e-03f,		2.002e-03f,		2.198e-03f,		2.402e-03f,		2.614e-03f,		2.834e-03f,		4.046e-03f,		5.438e-03f,		7.006e-03f,		8.742e-03f,
		1.064e-02f,		1.270e-02f,		1.492e-02f,		1.730e-02f,		2.250e-02f,		2.829e-02f,		3.466e-02f,		4.161e-02f,		4.910e-02f,		5.715e-02f,
		6.574e-02f,		7.486e-02f,		8.452e-02f,		9.469e-02f,		1.054e-01f,		1.166e-01f,		1.283e-01f,		1.405e-01f,		2.089e-01f,		2.893e-01f,
		3.813e-01f,		4.845e-01f,		7.238e-01f,		8.593e-01f,		1.005e+00f,		1.327e+00f,		1.688e+00f,		2.088e+00f,		2.524e+00f,		2.996e+00f,
		3.504e+00f,		4.045e+00f,		4.620e+00f,		5.228e+00f,		5.867e+00f,		6.538e+00f,		7.239e+00f,		7.970e+00f,		8.731e+00f,		1.295e+01f,
		1.783e+01f,		2.330e+01f,		2.933e+01f,		3.585e+01f,		4.284e+01f,		5.026e+01f,		5.809e+01f,		7.481e+01f,		9.282e+01f,		1.119e+02f,
		1.320e+02f,		1.530e+02f,		1.747e+02f,		1.971e+02f,		2.200e+02f,		2.435e+02f,		2.674e+02f,		2.918e+02f,		3.165e+02f,		3.415e+02f,
		3.668e+02f,		6.306e+02f,		9.035e+02f,		1.178e+03f,		1.453e+03f,		1.995e+03f,		2.527e+03f,		3.049e+03f,		3.561e+03f,		4.066e+03f,
		4.563e+03f,		5.053e+03f
    },
    {
        4.5550e-1f,    4.9060e-1f,    5.1970e-1f,    5.4400e-1f,    5.6470e-1f,    5.9860e-1f,    6.2540e-1f,    6.4730e-1f,    6.6560e-1f,    6.8130e-1f,
        6.9500e-1f,    7.0700e-1f,    7.3180e-1f,    7.5140e-1f,    7.6740e-1f,    7.8080e-1f,    7.9230e-1f,    8.0220e-1f,    8.1090e-1f,    8.1870e-1f,
        8.3190e-1f,    8.4290e-1f,    8.5220e-1f,    8.6020e-1f,    8.6730e-1f,    8.7350e-1f,    8.7910e-1f,    8.8420e-1f,    8.8890e-1f,    8.9310e-1f,
        8.9710e-1f,    9.0070e-1f,    9.0410e-1f,    9.0730e-1f,    9.2070e-1f,    9.3100e-1f,    9.3930e-1f,    9.4600e-1f,    9.5150e-1f,    9.5620e-1f,
        9.6010e-1f,    9.6350e-1f,    9.6890e-1f,    9.7310e-1f,    9.7640e-1f,    9.7900e-1f,    9.8110e-1f,    9.8290e-1f,    9.8440e-1f,    9.8570e-1f,
        9.8680e-1f,    9.8770e-1f,    9.8860e-1f,    9.8930e-1f,    9.8990e-1f,    9.9050e-1f,    9.9250e-1f,    9.9380e-1f,    9.9460e-1f,    9.9520e-1f,
        9.9570e-1f,    9.9600e-1f,    9.9630e-1f,    9.9650e-1f,    9.9680e-1f,    9.9710e-1f,    9.9730e-1f,    9.9740e-1f,    9.9750e-1f,    9.9760e-1f,
        9.9770e-1f,    9.9770e-1f,    9.9780e-1f,    9.9780e-1f,    9.9790e-1f,    9.9790e-1f,    9.9790e-1f,    9.9800e-1f,    9.9810e-1f,    9.9820e-1f,
        9.9820e-1f,    9.9830e-1f,    9.9830e-1f,    9.9840e-1f,    9.9840e-1f,    9.9840e-1f,    9.9850e-1f,    9.9850e-1f,    9.9850e-1f,    9.9850e-1f,
        9.9860e-1f,    9.9860e-1f,    9.9860e-1f,    9.9860e-1f,    9.9860e-1f,    9.9860e-1f,    9.9860e-1f,    9.9860e-1f,    9.9870e-1f,    9.9870e-1f,
        9.9870e-1f,    9.9880e-1f,    9.9880e-1f,    9.9880e-1f,    9.9890e-1f,    9.9890e-1f,    9.9890e-1f,    9.9890e-1f,    9.9900e-1f,    9.9900e-1f,
        9.9900e-1f,    9.9910e-1f,    9.9910e-1f,    9.9910e-1f,    9.9910e-1f,    9.9910e-1f,    9.9920e-1f,    9.9920e-1f,    9.9920e-1f,    9.9920e-1f,
        9.9920e-1f,    9.9930e-1f,    9.9940e-1f,    9.9950e-1f,    9.9950e-1f,    9.9960e-1f,    9.9960e-1f,    9.9970e-1f,    9.9970e-1f,    9.9970e-1f,
        9.9970e-1f,    9.9980e-1f,    2.4840e-1f,    2.8570e-1f,    3.1680e-1f,    3.4320e-1f,    3.6630e-1f,    4.0500e-1f,    4.3660e-1f,    4.6310e-1f,
        4.8600e-1f,    5.0590e-1f,    5.2350e-1f,    5.3930e-1f,    5.7260e-1f,    5.9970e-1f,    6.2220e-1f,    6.4150e-1f,    6.5810e-1f,    6.7270e-1f,
        6.8570e-1f,    6.9730e-1f,    7.1740e-1f,    7.3420e-1f,    7.4870e-1f,    7.6130e-1f,    7.7240e-1f,    7.8230e-1f,    7.9120e-1f,    7.9930e-1f,
        8.0670e-1f,    8.1360e-1f,    8.1990e-1f,    8.2580e-1f,    8.3130e-1f,    8.3640e-1f,    8.5800e-1f,    8.7480e-1f,    8.8830e-1f,    8.9940e-1f,
        9.0870e-1f,    9.1660e-1f,    9.2330e-1f,    9.2920e-1f,    9.3880e-1f,    9.4640e-1f,    9.5240e-1f,    9.5730e-1f,    9.6140e-1f,    9.6480e-1f,
        9.6770e-1f,    9.7020e-1f,    9.7240e-1f,    9.7430e-1f,    9.7590e-1f,    9.7730e-1f,    9.7860e-1f,    9.7980e-1f,    9.8400e-1f,    9.8670e-1f,
        9.8860e-1f,    9.8990e-1f,    9.9090e-1f,    9.9170e-1f,    9.9230e-1f,    9.9280e-1f,    9.9350e-1f,    9.9410e-1f,    9.9450e-1f,    9.9480e-1f,
        9.9510e-1f,    9.9530e-1f,    9.9540e-1f,    9.9560e-1f,    9.9570e-1f,    9.9580e-1f,    9.9590e-1f,    9.9600e-1f,    9.9610e-1f,    9.9620e-1f,
        9.9640e-1f,    9.9660e-1f,    9.9670e-1f,    9.9680e-1f,    9.9700e-1f,    9.9700e-1f,    9.9710e-1f,    9.9720e-1f,    9.9720e-1f,    9.9730e-1f,
        9.9730e-1f,    9.9740e-1f,    9.9740e-1f,    9.9750e-1f,    9.9750e-1f,    9.9750e-1f,    9.9760e-1f,    9.9760e-1f,    9.9760e-1f,    9.9760e-1f,
        9.9760e-1f,    9.9770e-1f,    9.9780e-1f,    9.9790e-1f,    9.9790e-1f,    9.9800e-1f,    9.9800e-1f,    9.9810e-1f,    9.9810e-1f,    9.9820e-1f,
        9.9820e-1f,    9.9830e-1f,    9.9830e-1f,    9.9840e-1f,    9.9840e-1f,    9.9840e-1f,    9.9850e-1f,    9.9850e-1f,    9.9850e-1f,    9.9860e-1f,
        9.9860e-1f,    9.9860e-1f,    9.9870e-1f,    9.9880e-1f,    9.9900e-1f,    9.9910e-1f,    9.9920e-1f,    9.9930e-1f,    9.9940e-1f,    9.9940e-1f,
        9.9950e-1f,    9.9950e-1f,    9.9950e-1f,    9.9960e-1f,    2.5550e-1f,    2.9330e-1f,    3.2450e-1f,    3.5090e-1f,    3.7380e-1f,    4.1220e-1f,
        4.4340e-1f,    4.6960e-1f,    4.9210e-1f,    5.1170e-1f,    5.2910e-1f,    5.4450e-1f,    5.7730e-1f,    6.0400e-1f,    6.2630e-1f,    6.4540e-1f,
        6.6200e-1f,    6.7670e-1f,    6.8980e-1f,    7.0150e-1f,    7.2200e-1f,    7.3940e-1f,    7.5440e-1f,    7.6760e-1f,    7.7930e-1f,    7.8980e-1f,
        7.9940e-1f,    8.0810e-1f,    8.1610e-1f,    8.2340e-1f,    8.3020e-1f,    8.3660e-1f,    8.4250e-1f,    8.4800e-1f,    8.7080e-1f,    8.8800e-1f,
        9.0130e-1f,    9.1190e-1f,    9.2050e-1f,    9.2770e-1f,    9.3370e-1f,    9.3890e-1f,    9.4710e-1f,    9.5350e-1f,    9.5850e-1f,    9.6260e-1f,
        9.6590e-1f,    9.6870e-1f,    9.7110e-1f,    9.7310e-1f,    9.7490e-1f,    9.7640e-1f,    9.7780e-1f,    9.7900e-1f,    9.8000e-1f,    9.8100e-1f,
        9.8450e-1f,    9.8680e-1f,    9.8840e-1f,    9.8950e-1f,    9.9040e-1f,    9.9110e-1f,    9.9160e-1f,    9.9210e-1f,    9.9280e-1f,    9.9330e-1f,
        9.9370e-1f,    9.9400e-1f,    9.9420e-1f,    9.9440e-1f,    9.9460e-1f,    9.9480e-1f,    9.9490e-1f,    9.9500e-1f,    9.9510e-1f,    9.9520e-1f,
        9.9530e-1f,    9.9530e-1f,    9.9560e-1f,    9.9580e-1f,    9.9600e-1f,    9.9610e-1f,    9.9630e-1f,    9.9630e-1f,    9.9640e-1f,    9.9650e-1f,
        9.9660e-1f,    9.9660e-1f,    9.9670e-1f,    9.9680e-1f,    9.9680e-1f,    9.9680e-1f,    9.9690e-1f,    9.9690e-1f,    9.9700e-1f,    9.9700e-1f,
        9.9700e-1f,    9.9700e-1f,    9.9710e-1f,    9.9720e-1f,    9.9730e-1f,    9.9730e-1f,    9.9740e-1f,    9.9750e-1f,    9.9750e-1f,    9.9760e-1f,
        9.9760e-1f,    9.9770e-1f,    9.9780e-1f,    9.9790e-1f,    9.9790e-1f,    9.9800e-1f,    9.9800e-1f,    9.9810e-1f,    9.9810e-1f,    9.9820e-1f,
        9.9820e-1f,    9.9820e-1f,    9.9830e-1f,    9.9830e-1f,    9.9830e-1f,    9.9860e-1f,    9.9870e-1f,    9.9890e-1f,    9.9900e-1f,    9.9910e-1f,
        9.9920e-1f,    9.9930e-1f,    9.9930e-1f,    9.9940e-1f,    9.9940e-1f,    9.9950e-1f,    4.7130e-1f,    5.1250e-1f,    5.4510e-1f,    5.7140e-1f,
        5.9340e-1f,    6.2870e-1f,    6.5600e-1f,    6.7790e-1f,    6.9610e-1f,    7.1160e-1f,    7.2490e-1f,    7.3650e-1f,    7.6040e-1f,    7.7910e-1f,
        7.9430e-1f,    8.0690e-1f,    8.1770e-1f,    8.2700e-1f,    8.3510e-1f,    8.4240e-1f,    8.5460e-1f,    8.6480e-1f,    8.7330e-1f,    8.8070e-1f,
        8.8710e-1f,    8.9280e-1f,    8.9780e-1f,    9.0240e-1f,    9.0660e-1f,    9.1040e-1f,    9.1390e-1f,    9.1710e-1f,    9.2020e-1f,    9.2300e-1f,
        9.3470e-1f,    9.4360e-1f,    9.5060e-1f,    9.5630e-1f,    9.6110e-1f,    9.6500e-1f,    9.6840e-1f,    9.7130e-1f,    9.7590e-1f,    9.7930e-1f,
        9.8200e-1f,    9.8420e-1f,    9.8590e-1f,    9.8720e-1f,    9.8840e-1f,    9.8930e-1f,    9.9010e-1f,    9.9080e-1f,    9.9140e-1f,    9.9200e-1f,
        9.9240e-1f,    9.9280e-1f,    9.9430e-1f,    9.9520e-1f,    9.9580e-1f,    9.9620e-1f,    9.9650e-1f,    9.9680e-1f,    9.9700e-1f,    9.9710e-1f,
        9.9740e-1f,    9.9750e-1f,    9.9770e-1f,    9.9780e-1f,    9.9780e-1f,    9.9790e-1f,    9.9800e-1f,    9.9800e-1f,    9.9810e-1f,    9.9810e-1f,
        9.9810e-1f,    9.9820e-1f,    9.9820e-1f,    9.9820e-1f,    9.9830e-1f,    9.9840e-1f,    9.9840e-1f,    9.9850e-1f,    9.9850e-1f,    9.9850e-1f,
        9.9860e-1f,    9.9860e-1f,    9.9860e-1f,    9.9870e-1f,    9.9870e-1f,    9.9870e-1f,    9.9870e-1f,    9.9870e-1f,    9.9870e-1f,    9.9880e-1f,
        9.9880e-1f,    9.9880e-1f,    9.9880e-1f,    9.9880e-1f,    9.9880e-1f,    9.9880e-1f,    9.9890e-1f,    9.9890e-1f,    9.9890e-1f,    9.9900e-1f,
        9.9900e-1f,    9.9900e-1f,    9.9900e-1f,    9.9910e-1f,    9.9910e-1f,    9.9910e-1f,    9.9910e-1f,    9.9920e-1f,    9.9920e-1f,    9.9920e-1f,
        9.9920e-1f,    9.9920e-1f,    9.9920e-1f,    9.9930e-1f,    9.9930e-1f,    9.9930e-1f,    9.9930e-1f,    9.9940e-1f,    9.9950e-1f,    9.9950e-1f,
        9.9960e-1f,    9.9960e-1f,    9.9970e-1f,    9.9970e-1f,    9.9970e-1f,    9.9970e-1f,    9.9980e-1f,    9.9980e-1f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,            -1.000e+04f,
        -1.000e+04f,            -1.000e+04f,
        3.105e-01f,             3.520e-01f,             3.859e-01f,             4.141e-01f,             4.384e-01f,             4.782e-01f,             5.101e-01f,             5.364e-01f,             5.586e-01f,             5.778e-01f,
        5.946e-01f,             6.095e-01f,             6.405e-01f,             6.653e-01f,             6.857e-01f,             7.029e-01f,             7.176e-01f,             7.305e-01f,             7.418e-01f,             7.519e-01f,
        7.692e-01f,             7.836e-01f,             7.959e-01f,             8.064e-01f,             8.157e-01f,             8.240e-01f,             8.314e-01f,             8.381e-01f,             8.442e-01f,             8.498e-01f,
        8.549e-01f,             8.597e-01f,             8.641e-01f,             8.683e-01f,             8.857e-01f,             8.991e-01f,             9.099e-01f,             9.188e-01f,             9.264e-01f,             9.328e-01f,
        9.384e-01f,             9.432e-01f,             9.513e-01f,             9.576e-01f,             9.626e-01f,             9.667e-01f,             9.701e-01f,             9.729e-01f,             9.753e-01f,             9.773e-01f,
        9.791e-01f,             9.806e-01f,             9.819e-01f,             9.831e-01f,             9.841e-01f,             9.850e-01f,             9.884e-01f,             9.905e-01f,             9.919e-01f,             9.930e-01f,
        9.937e-01f,             9.943e-01f,             9.948e-01f,             9.951e-01f,             9.957e-01f,             9.961e-01f,             9.963e-01f,             9.966e-01f,             9.967e-01f,             9.969e-01f,
        9.970e-01f,             9.971e-01f,             9.972e-01f,             9.973e-01f,             9.973e-01f,             9.974e-01f,             9.974e-01f,             9.975e-01f,             9.977e-01f,             9.978e-01f,
        9.979e-01f,             9.979e-01f,             9.980e-01f,             9.980e-01f,             9.981e-01f,             9.981e-01f,             9.982e-01f,             9.982e-01f,             9.982e-01f,             9.983e-01f,
        9.983e-01f,             9.983e-01f,             9.983e-01f,             9.984e-01f,             9.984e-01f,             9.984e-01f,             9.984e-01f,             9.984e-01f,             9.984e-01f,             9.985e-01f,
        9.985e-01f,             9.986e-01f,             9.986e-01f,             9.986e-01f,             9.987e-01f,             9.987e-01f,             9.987e-01f,             9.988e-01f,             9.988e-01f,             9.988e-01f,
        9.989e-01f,             9.989e-01f,             9.989e-01f,             9.989e-01f,             9.990e-01f,             9.990e-01f,             9.990e-01f,             9.990e-01f,             9.991e-01f,             9.991e-01f,
        9.991e-01f,             9.992e-01f,             9.993e-01f,             9.994e-01f,             9.994e-01f,             9.995e-01f,             9.996e-01f,             9.996e-01f,             9.996e-01f,             9.997e-01f,
        9.997e-01f,             9.997e-01f,
		3.304e-01f,		3.760e-01f,		4.123e-01f,		4.421e-01f,		4.674e-01f,		5.084e-01f,		5.406e-01f,		5.669e-01f,		5.889e-01f,		6.078e-01f,
		6.242e-01f,		6.387e-01f,		6.686e-01f,		6.923e-01f,		7.118e-01f,		7.281e-01f,		7.420e-01f,		7.542e-01f,		7.648e-01f,		7.743e-01f,
		7.906e-01f,		8.041e-01f,		8.156e-01f,		8.256e-01f,		8.343e-01f,		8.420e-01f,		8.490e-01f,		8.553e-01f,		8.611e-01f,		8.664e-01f,
		8.712e-01f,		8.758e-01f,		8.800e-01f,		8.839e-01f,		9.005e-01f,		9.133e-01f,		9.236e-01f,		9.320e-01f,		9.390e-01f,		9.450e-01f,
		9.500e-01f,		9.544e-01f,		9.614e-01f,		9.668e-01f,		9.710e-01f,		9.743e-01f,		9.770e-01f,		9.793e-01f,		9.812e-01f,		9.827e-01f,
		9.841e-01f,		9.852e-01f,		9.862e-01f,		9.871e-01f,		9.879e-01f,		9.886e-01f,		9.910e-01f,		9.926e-01f,		9.936e-01f,		9.943e-01f,
		9.949e-01f,		9.953e-01f,		9.956e-01f,		9.959e-01f,		9.963e-01f,		9.966e-01f,		9.968e-01f,		9.969e-01f,		9.971e-01f,		9.972e-01f,
		9.973e-01f,		9.973e-01f,		9.974e-01f,		9.975e-01f,		9.975e-01f,		9.976e-01f,		9.976e-01f,		9.976e-01f,		9.978e-01f,		9.979e-01f,
		9.979e-01f,		9.980e-01f,		9.981e-01f,		9.981e-01f,		9.981e-01f,		9.982e-01f,		9.982e-01f,		9.983e-01f,		9.983e-01f,		9.983e-01f,
		9.983e-01f,		9.984e-01f,		9.984e-01f,		9.984e-01f,		9.984e-01f,		9.984e-01f,		9.984e-01f,		9.985e-01f,		9.985e-01f,		9.985e-01f,
		9.986e-01f,		9.986e-01f,		9.986e-01f,		9.987e-01f,		9.987e-01f,		9.987e-01f,		9.987e-01f,		9.988e-01f,		9.988e-01f,		9.989e-01f,
		9.989e-01f,		9.989e-01f,		9.990e-01f,		9.990e-01f,		9.990e-01f,		9.990e-01f,		9.990e-01f,		9.991e-01f,		9.991e-01f,		9.991e-01f,
		9.991e-01f,		9.992e-01f,		9.993e-01f,		9.994e-01f,		9.994e-01f,		9.995e-01f,		9.996e-01f,		9.996e-01f,		9.997e-01f,		9.997e-01f,
		9.997e-01f,		9.997e-01f
    },
    {
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,    Water_Liquid,
        Water_Liquid,    Water_Liquid,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,
        Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum_Oxide,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,
        Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    Aluminum,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,    PMMA,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,    Alanine,
        Alanine,    Alanine,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,    LiF,
        LiF,    LiF,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,    Air,
        Air,    Air
    }
};

/** ----------------------------------------------- Static objects  --------------------------------------- */


static const AT_stopping_power_sources_struct AT_stopping_power_sources = {
	STOPPING_POWER_SOURCE_N,
    {  PSTAR,         Bethe, 			ShieldHit, 				ICRU,			FLUKA, ATIMA}, // source_no
    {  "PSTAR data",  "Bethe formula",  "ShieldHit (Bethe)",    "ICRU 49&73",	"FLUKA dedx I=77.3eV", "ATIMA"}    // source_name
};


static const AT_stopping_power_tabulated_source_group_for_all_materials_struct AT_data_PSTAR_source = {
  MATERIAL_DATA_N,
  PSTAR,
  { User_Defined_Material,  Water_Liquid,                         Aluminum_Oxide,                                 Aluminum,                                 PMMA,                                Alanine,                                LiF,                                Air,								Silicon,								Copper },
  { NULL,                   &AT_stopping_power_data_PSTAR_Water,  &AT_stopping_power_data_PSTAR_Aluminum_Oxide,   &AT_stopping_power_data_PSTAR_Aluminum,   &AT_stopping_power_data_PSTAR_PMMA,  &AT_stopping_power_data_PSTAR_Alanine,  &AT_stopping_power_data_PSTAR_LiF,  &AT_stopping_power_data_PSTAR_Air,	&AT_stopping_power_data_PSTAR_Silicon,	&AT_stopping_power_data_PSTAR_Copper}
};


static const AT_stopping_power_tabulated_source_group_struct AT_stopping_power_tabulated_source = {
  STOPPING_POWER_SOURCE_N,
  { PSTAR,                  Bethe, 	ShieldHit,	ICRU,	FLUKA, ATIMA},
  { &AT_data_PSTAR_source,  NULL,	NULL,       NULL,	NULL, NULL}
};


static const AT_stopping_power_analytical_sources_struct AT_stopping_power_analytical_source = {
  STOPPING_POWER_SOURCE_N,
  { PSTAR,   Bethe,				ShieldHit,             ICRU,				FLUKA, ATIMA},
  { NULL,    &AT_Bethe_wrapper, &AT_ShieldHit_wrapper, &AT_ICRU_wrapper,		&AT_FLUKA_wrapper, &AT_ATIMA_wrapper}
};

//////////////////////////////////////
// FORMER AT_DataLET, integrate here
//////////////////////////////////////

/**
 * TODO find better name for this function
 * it just finds a vector of interpolated values
 * and is not related only to PSTAR tables
 *
 * @param[in]  n             number of points to interpolate
 * @param[in]  x             array of x values for which interpolation is done
 * @param[in]  subset_no     TODO
 * @param[in]  x_table       x part of data table
 * @param[in]  y_table       y part of data table
 * @param[out] y             array of interpolated y values
 */
void get_table_value( const long    n,
    const double  x[],
    const long    subset_no,
    const double  x_table[],
    const double  y_table[],
    double        y[]);


/**
 * Returns CSDA range (in g/cm2) from pstar tables for given energy.
 * In case of ions a simple scaling procedure (A/Z^2) will be used (even effective charge will be neglected)
 * @param[in]   E_MeV_u                  energy of particle
 * @param[in]   particle_no              type of the particle
 * @see          AT_DataParticle.h for definition
 * @param[in]   material_no              material index
 * @see          AT_DataMaterial.h for definition
 * @return      CSDA_range_g_cm2         CSDA range
  */
double AT_CSDA_range_g_cm2_single(   const double  E_MeV_u,
    const long    particle_no,
    const long    material_no);


/**
 * Returns CSDA range (in g/cm2) from pstar tables for given energy.
 * In case of ions a simple scaling procedure (A/Z^2) will be used (even effective charge will be neglected)
 * @param[in]   number_of_particles      number of particle types in the mixed particle field
 * @param[in]   E_MeV_u                  energy of particles in the mixed particle field (array of size number_of_particles)
 * @param[in]   particle_no              type of the particles in the mixed particle field (array of size number_of_particles)
 * @see          AT_DataParticle.h for definition
 * @param[in]   material_no              material index
 * @see          AT_DataMaterial.h for definition
 * @param[out]  CSDA_range_g_cm2         (array of size number_of_particles) to be allocated by the user which will be used to return the results
 */
void AT_CSDA_range_g_cm2(  const long  number_of_particles,
    const double   E_MeV_u[],
    const long     particle_no[],
    const long     material_no,
    double         CSDA_range_g_cm2[]);


/**
 * Returns CSDA range (in m) from pstar tables for given energy.
 * In case of ions a simple scaling procedure (A/Z^2) will be used (even effective charge will be neglected)
 * @param[in]   E_MeV_u                  energy of particle
 * @param[in]   particle_no              type of the particle
 * @see          AT_DataParticle.h for definition
 * @param[in]   material_no              material index
 * @see          AT_DataMaterial.h for definition
 * @return      CSDA_range_m             CSDA range
 */
double AT_CSDA_range_m_single(  const double  E_MeV_u,
    const long    particle_no,
    const long    material_no);

/**
 * Returns CSDA range (in m) from pstar tables for given energy.
 * In case of ions a simple scaling procedure (A/Z^2) will be used (even effective charge will be neglected)
 * @param[in]   number_of_particles      number of particle types in the mixed particle field
 * @param[in]   E_MeV_u                  energy of particles in the mixed particle field (array of size number_of_particles)
 * @param[in]   particle_no              type of the particles in the mixed particle field (array of size number_of_particles)
 * @see          AT_DataParticle.h for definition
 * @param[in]   material_no              material index
 * @see          AT_DataMaterial.h for definition
 * @param[out]  CSDA_range_m            (array of size number_of_particles) to be allocated by the user which will be used to return the results
 */
void AT_CSDA_range_m(  const long  number_of_particles,
	    const double  E_MeV_u[],
	    const long    particle_no[],
	    const long    material_no,
	    double        CSDA_range_m[]);

#endif /* AT_DATASTOPPINGPOWER_H_ */
