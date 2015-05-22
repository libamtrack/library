#ifndef AT_DATAMATERIAL_H_
#define AT_DATAMATERIAL_H_

/**
 * @brief Material properties
 */

/*
 *    AT_DataMaterial.h
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
#include <math.h>

#include "AT_Error.h"
#include "AT_NumericalRoutines.h"

/**
 * @enum material_no Materials code numbers
 */
enum material_no{
  User_Defined_Material     		= 0,  /**< To be defined by the user during runtime >**/
  Water_Liquid             			= 1,  /**< Liquid water */
  Aluminum_Oxide            		= 2,  /**< Aluminium oxide */
  Aluminum                  		= 3,  /**< Aluminium */
  PMMA                      		= 4,  /**< PMMA */
  Alanine                   		= 5,  /**< Alanine */
  LiF                       		= 6,  /**< Lithium Fluoride */
  Air				        		= 7,  /**< Air dry (at sea level) */
  Silicon                    		= 8,  /**< Silicon */
  Copper                    		= 9,  /**< Copper */
  Tungsten                   		= 10, /**< Tungsten */
  Gammex_Lung_LN450         		= 11, /**< Gammex tissue surrogate "Lung (LN450)" */
  Gammex_AP6_Adipose_RMI453 		= 12, /**< Gammex tissue surrogate "AP6 Adipose RMI 453" */
  Gammex_BR12_Breast_RMI454 		= 13, /**< Gammex tissue surrogate "BR12 Breast RMI454" */
  Gammex_CT_Solid_Water_RMI451		= 14, /**< Gammex tissue surrogate "CT Solid Water RMI451" */
  Gammex_Water				 		= 15, /**< Gammex tissue surrogate "Gammex Water" */
  Gammex_Muscle_RMI452         		= 16, /**< Gammex tissue surrogate "Muscle RMI452" */
  Gammex_LV1_RMI 					= 17, /**< Gammex tissue surrogate "LV1 Liver RMI" */
  Gammex_SR2_Brain 			 		= 18, /**< Gammex tissue surrogate "SR2 Brain" */
  Gammex_IB3_Inner_Bone_RMI456		= 19, /**< Gammex tissue surrogate "IB3 Inner Bone RMI 456" */
  Gammex_B200_Bone_Mineral			= 20, /**< Gammex tissue surrogate "B200 Bone Mineral" */
  Gammex_CB2_30_CaCO3				= 21, /**< Gammex tissue surrogate "CB2 30%  CaCO3" */
  Gammex_CB2_50_CaCO3				= 22, /**< Gammex tissue surrogate "CB2 50%  CaCO3" */
  Gammex_SB3_Cortical_Bone_RMI450	= 23, /**< Gammex tissue surrogate "SB3 Cortical Bone RMI 450" */
  Lead                              = 24  /**< Lead */
};

enum material_phase{
  phase_undefined		= 0,
  phase_condensed     	= 1,
  phase_gaseous			= 2
};

#define MATERIAL_DATA_N    25

// TODO The next two LET-related structures must be declared here rather than in AT_DataLET.h to avoid circular dependencies

/**
 * @struct AT_LET_data_single
 * Stopping power data for a material and a particle type
 */
typedef struct {
  long     n;                                     /**< number of items in the data table */
  long     particle_no;                           /**< particle number this data is for, see AT_DataParticle.h for definition */
  double*  kin_E_MeV;                             /**< Kinetic energy in MeV, pointer to array of size n */
  double*  stp_pow_el_MeV_cm2_g;                  /**< Electronic (Collision) Stopping Power, pointer to array of size n */
  double*  range_cdsa_g_cm2;                      /**< CSDA (continuous-slowing-down approximation) range, pointer to array of size n */
} AT_LET_data_single;

/**
 * @struct AT_LET_data
 * Stopping power data for a material
 */
typedef struct {
  long            n;                              /**< number of data tables in the structure */
  AT_LET_data_single*   LET_data_single;          /**< pointer to LET data tables for a particle type */
} AT_LET_data;

#define MATERIAL_NAME_LENGTH    255


/**
 * @struct AT_single_material_data_struct
 * TODO
 */
typedef struct {
  long            material_no;                            /**< material number - if 0 the user has to specify the material properties by at least handing over the data marked 'essential', all other data will be computed if not overridden by the user. If a number > 0 is used, a predefined material will be loaded */
  bool            material_established;                   /**< if true, the material has been established, MUST BE SET TO "FALSE" BY USER */

  double          density_g_cm3;                          /**< physical density in g/cm3 [ESSENTIAL] */
  double          I_eV;                                   /**< Ionization potential in eV [ESSENTIAL] */ //TODO: check if this could not be calculated

  long            n_elements;                             /**< number of elements constituting the material [ESSENTIAL] */
  long*           elements_Z;                             /**< Atomic numbers Z for the constituting elements, pointer to array of length n_elements [ESSENTIAL] */
  long*           elements_A;                             /**< Mass numbers A for the constituting elements, pointer to array of length n_elements [ESSENTIAL] */
  double*         elements_weight_fraction;               /**< Weight fractions for the constituting elements, pointer to array of length n_elements, have to add up to 1.0 [ESSENTIAL] */

  char*           material_name;                          /**< material name */
  long            ICRU_ID;                                /**< ICRU ID of the material, might serve for automatic look-up later */

  double          average_Z;                              /**< Average atomic number, will be calculated if not given */
  double          average_A;                              /**< Average mass number,  will be calculated if not given */

  long            LET_data_source;                        /**< Defines the source for stopping power data, see enum AT_LET_data_source [ESSENTIAL]. The user has to make sure that he provides the necessary data. */

  double          p_MeV;                                  /**< Prefactor for the power-law description of stopping power: S = p*E^alpha. In MeV^(1/alpha) */
  double          alpha_g_cm2_MeV;                        /**< Exponent for the power-law description of stopping */

  AT_LET_data     LET_data;                               /**< Stopping power data for the material, see AT_DataLET.h for definition. Will be calculated in case of power-law or given by the user (or read in from PSTAR data?). */
                                                          /**<  LET_data has to hold at_least the proton stopping powers. For any particle where no explicit stopping power data is given, a scaling by Z_eff^2 will be performed. */

  long			  phase;								  /**< Phase of the material: condensed or gaseous */
//  const double  m_g_cm2;                                /**< Slope of linear approximation in fluence reduction due to nuclear interactions, not used yet */
} AT_single_material_data_struct;


/**
 * @struct AT_table_of_material_data_struct
 * TODO
 */
typedef struct {
  const long    n;
  const long    material_no[MATERIAL_DATA_N];
  bool          ready[MATERIAL_DATA_N];
  const long    ICRU_ID[MATERIAL_DATA_N];
  double        density_g_cm3[MATERIAL_DATA_N];
  double        I_eV[MATERIAL_DATA_N];
  const double  alpha_g_cm2_MeV[MATERIAL_DATA_N];
  const double  p_MeV[MATERIAL_DATA_N];
  const double  m_g_cm2[MATERIAL_DATA_N];
  double        average_A[MATERIAL_DATA_N];
  double        average_Z[MATERIAL_DATA_N];
  const char*   material_name[MATERIAL_DATA_N];
  long			phase[MATERIAL_DATA_N];
} AT_table_of_material_data_struct;

// This is to suppress "defined but not used" warnings by gcc here
#ifdef __GNUC__
#define VARIABLE_IS_NOT_USED __attribute__ ((unused))
#else
#define VARIABLE_IS_NOT_USED
#endif

static AT_table_of_material_data_struct VARIABLE_IS_NOT_USED AT_Material_Data = {
    MATERIAL_DATA_N,
    // material_no
    {  User_Defined_Material,
       Water_Liquid,      Aluminum_Oxide,   Aluminum,     PMMA,      Alanine,
       LiF,	              Air,              Silicon,      Copper,    Tungsten,
       Gammex_Lung_LN450,		Gammex_AP6_Adipose_RMI453,	Gammex_BR12_Breast_RMI454,		  Gammex_CT_Solid_Water_RMI451,		Gammex_Water,
	   Gammex_Muscle_RMI452,	Gammex_LV1_RMI,				Gammex_SR2_Brain, 		   		  Gammex_IB3_Inner_Bone_RMI456,	   	Gammex_B200_Bone_Mineral,
	   Gammex_CB2_30_CaCO3,		Gammex_CB2_50_CaCO3,	    Gammex_SB3_Cortical_Bone_RMI450,  Lead},

    // ready
    {  false,
       true,              true,             true,         true,      true,
       true,			  true,             true,         true,      true,
	   true,              true,             true,         true,      true,
       true,			  true,             true,         true,      true,
	   true,              true,				true,         true},
    // ICRU_ID
    {  0,
       276,               106,              13,           223,       0,
       185,               104,              14,           29,       74,
       0,                 0,              	0,            0,         0,
       0,              	  0,              	0,            0,         0,
       0,                 0,				0,            82},
    // density_g_cm3
    {  0.0,
       1.00,              3.97,             2.6989,       1.188,     1.42,
       2.64,   			  1.20479E-03,      2.33,         8.96,      19.3,
       0.450,             0.920,            0.980,        1.015,     1.0,
       1.050,          	  1.039,         	1.049,        1.133,   	 1.145,
       1.340,	          1.560,			1.819,        11.35},
    // I_eV
    {  0.0,
       75.0,              145.2,            166.0,        74.0,      71.9,
       10.0,              85.7,             173.0,        322.0,     727.0,
       71.45,             65.38,            66.81,        68.72,     68.86,
       68.55,             68.64,            62.51,        77.01,     77.09,
       77.48,             87.94,			97.40,        823.0},
    // alpha_g_cm2 - TODO No data for LiF, air
    {  0.0,
       0.00231,           0.003058,         0.003266,     0.001988,  0.00216381,
       0.0,	  		      0.0,              0.0,          0.0,       0.0,
       0.0,               0.0,              0.0,          0.0,       0.0,
       0.0,            	  0.0,              0.0,          0.0,       0.0,
       0.0,               0.0,				0.0,          0.0},
	// p_MeV - TODO No data for LiF, air
    {  0.0,
       1.761,             1.748,            1.745,        1.762,     1.79165987,
       0.0,	  		      0.0,              0.0,          0.0,       0.0,
       0.0,               0.0,              0.0,          0.0,       0.0,
       0.0,            	  0.0,              0.0,          0.0,       0.0,
       0.0,               0.0,				0.0,          0.0},
    // m_g_cm2 - TODO No data processed for nuclear interactions in Alanine, hence set to -100
    {  0.0,
       0.01153,           0.01305,          0.01230,      0.01338,   -100.0,
	   0.0,	 			  0.0,              0.0,          0.0,       0.0,
       0.0,               0.0,              0.0,          0.0,       0.0,
       0.0,            	  0.0,              0.0,          0.0,       0.0,
       0.0,               0.0,				0.0,          0.0},
    // average_A
    {  0.0,
       13.0,              21.72,            27.0,         11.556,    12.8088,
       17.7333,           14.78,            28.085,       63.546,    183.84,
       12.33, 		      10.84,		 	11.24,        11.78,     12.98,
       11.76,          	  11.77,          	10.43,        14.45,  	 14.46,
       14.80,             17.62,			19.99,        207.2},
    // average_Z
    {  0.0,
       7.22,              10.637,           13.0,         6.24,      6.44,
       8.0,		          7.375,            14.0,         29.0,      74.0,
       6.68,     	      5.91,      	    6.10,    	  6.36,      7.22,
       6.36,           	  6.36,           	5.78,         7.70,      7.71,
       7.89,              9.23,				10.34,        82.0},
    // material_name
    {  "User defined",
       "Water, Liquid",   "Aluminum Oxide", "Aluminum",   "PMMA",    "Alanine",
       "Lithium Fluoride","Air", "Silicon", "Copper", "Tungsten",
       "Gammex Lung LN450", "Gammex AP6 Adipose RMI453", "BR12 Breast RMI454",	"CT Solid Water RMI451", "Gammex Water",
		"Muscle RMI452" , "LV1 Liver RMI" ,	"SR2 Brain", "IB3 Inner Bone RMI 456", "B200 Bone Mineral",
		"CB2 30% CaCO3", "CB2 50% CaCO3", "SB3 Cortical Bone RMI 450", "Lead" },
	// phase
	{  phase_undefined,
	   phase_condensed, phase_condensed, phase_condensed, phase_condensed, phase_condensed,
	   phase_condensed, phase_gaseous,   phase_condensed, phase_condensed, phase_condensed,
	   phase_condensed, phase_condensed, phase_condensed, phase_condensed, phase_condensed,
	   phase_condensed, phase_condensed, phase_condensed, phase_condensed, phase_condensed,
	   phase_condensed, phase_condensed, phase_condensed, phase_condensed}
};

/* Cucinnotta calculated average A for water as 14.3, but it seems that it is 13.0 (Leszek) *
 * looks like error in his article */


/**
 * Get index of material in AT_Material_Data for given material_no
 * (currently for example material with number 2 has index 1)
 *
 * @param[in] material_number  material number
 * @return                 material index in AT_Material_Data table
 */
long AT_material_index_from_material_number( const long material_number );


/**
 * Get material name
 * @param[in]  material_no
 * @param[out] material_name
 */
void AT_material_name_from_number( const long material_no,
    char* material_name);


/**
 * Get material number
 * @param[in] material_name
 * @return    material number
 */
long AT_material_number_from_name( const char* material_name );


/**
 * Get material density [g/cm3] for single material with number material_no
 * @param[in] material_no
 * @return    material density [g/cm3]
 */
double AT_density_g_cm3_from_material_no( const long   material_no );




/**
 * Get mean ionization potential in eV for single material with number material_no
 * @param[in] material_no
 * @return    mean ionization potential [eV]
 */
double AT_I_eV_from_material_no( const long   material_no );


/**
 * Get fit parameter for power-law representation of stp.power/range/E-dependence for single material with number material_no
 * @param[in] material_no
 * @return    fit parameter for power-law representation of stp.power/range/E-dependence
 */
double AT_alpha_g_cm2_MeV_from_material_no( const long   material_no );


/**
 * Get fit parameter for power-law representation of stp.power/range/E-dependence for single material with number material_no
 * @param[in] material_no
 * @return    fit parameter for power-law representation of stp.power/range/E-dependence
 */
double AT_p_MeV_from_material_no( const long   material_no );


/**
 * Get fit parameter for the linear representation of fluence changes due to nuclear interactions based on data from Janni for single material with number material_no
 * @param[in] material_no
 * @return    fit parameter for the linear representation of fluence changes due to nuclear interactions based on data from Janni
 */
double AT_m_g_cm2_from_material_no( const long   material_no );


/**
 * Get average mass number for single material with number material_no
 * @param[in] material_no
 * @return    average mass number
 */
double AT_average_A_from_material_no( const long   material_no );


/**
 * Get average atomic number for single material with number material_no
 * @param[in] material_no
 * @return    average atomic number
 */
double AT_average_Z_from_material_no( const long   material_no );

/**
 * Get material phase
 * @param[in] material_no
 * @return    material phase index (enumerator)
 */
long AT_phase_from_material_no( const long   material_no );

/**
 * Returns material data for single material
 * @param[in]  material_no           index number of material
 * @param[out] density_g_cm3         physical density [g/cm3]
 *   electron_density = number_of_electron_per_molecule * avogadro_constant / molar_mass * 1m3 * density
 * @param[out] I_eV                  mean ionization potential [eV]
 * @param[out] alpha_g_cm2_MeV       fit parameter for power-law representation of stopping-power (Bortfeld, 1997)
 *   @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
 *   Here, however, we use the mass stopping power. The correct dimension is g/(cm^2 * MeV^p)
 * @param[out] p_MeV                 fit parameter for power-law representation of stopping-power (Bortfeld, 1997)
 *   @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.\n
 *   p is actually dimensionless, it should be nevertheless indicated, that the energy must be given in MeV
 * @param[out] m_g_cm2               fit parameter for the linear representation of fluence changes due to nuclear interactions based on data from Janni 1982  (Bortfeld, 1997)
 *   @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
 * @param[out] average_A             average mass number
 *   let f_i be fraction by weight of the constituent element with atomic number Z_i and atomic weight A_i\n
 *   let us define average_Z/A = sum_i f_i Z_i / A_i \n
 *   then we have: average_A = average_Z / (average_Z/A) \n
 *   for water (H20) we have: average_Z/A = (2/18) * (1/1) + (16/18)*(8/16) = 0.5555 \n
 *   average_A = 7.22 / 0.555 = 13
 * @param[out] average_Z             average atomic number
 *   let f_i be fraction by weight of the constituent element with atomic number Z_i and atomic weight A_i\n
 *   average_Z = sum_i f_i Z_i \n
 *   for water (H20) we have: average_Z = (2/18)*1 + (16/18)*8 = 7.22\n
 *   @see Tabata, T. (1972) Generalized semiempirical equations for the extrapolated range of electrons, Nucl. Instr and Meth. 103, 85-91.
 */
void AT_get_material_data(     const long  material_no,
    double*  density_g_cm3,
    double*  I_eV,
    double*  alpha_g_cm2_MeV,
    double*  p_MeV,
    double*  m_g_cm2,
    double*  average_A,
    double*  average_Z);


/**
 * Returns material data for list of materials
 * @param[in]   number_of_materials  numbers of materials the routine is called for
 * @param[in]   material_no          material indices (array of size number_of_materials)
 * @param[out]  density_g_cm3        material density in g/cm3 (array of size number_of_materials)
 * @param[out]  I_eV                 mean ionization potential in eV (array of size number_of_materials)
 * @param[out]  alpha_g_cm2_MeV      fit parameter for power-law representation of stp.power/range/E-dependence (array of size number_of_materials)
 * @param[out]  p_MeV                fit parameter for power-law representation of stp.power/range/E-dependence (array of size number_of_materials)
 * @param[out]  m_g_cm2              fit parameter for the linear representation of fluence changes due to nuclear interactions based on data from Janni 1982 (array of size number_of_materials)
 * @param[out]  average_A            average mass number (array of size number_of_materials)
 * @param[out]  average_Z            average atomic number (array of size number_of_materials)
 */
void AT_get_materials_data( const long  number_of_materials,
    const long  material_no[],
    double  density_g_cm3[],
    double  I_eV[],
    double  alpha_g_cm2_MeV[],
    double  p_MeV[],
    double  m_g_cm2[],
    double  average_A[],
    double  average_Z[]);

// ROUTINES TO GET DERIVED PARAMETERS FOR A MATERIAL
/**
 * Get electron density [1/m3] for single material with number material_no
 * @param[in] material_no
 * @return    electron density [1/m3]
 */
double AT_electron_density_m3_from_material_no_single( const long   material_no );

/**
 * Get electron density [1/m3] for materials
 * @param[in]     n                    number of materials
 * @param[in]     material_no          material indices (array of size n)
 * @param[out]    electron_density_m3  electron densities per m3 (array of size n)
 */
void AT_electron_density_m3_from_material_no_multi( const long n,
		const long   material_no[],
		double electron_density_m3[]);

/**
 * Returns material's plasma energy needed for Sternheimer
 * computation of density effect in stopping power
 * @param[in]  material_no  material number
 */
double AT_plasma_energy_J_from_material_no( const long material_no );

// ROUTINES FOR COMPUTING DERIVED PARAMETERS
/**
 * Computes the electron density from average A and Z
 * @param[in] density_g_cm3      physical density (in g/cm3) of material
 * @param[in] average_Z          average atomic number of material
 * @param[in] average_A          mass number of material
 * @return                       electron density in 1/m3
 */
double AT_electron_density_m3_single( const double density_g_cm3,
    const double average_Z,
    const double average_A );

/**
 * Returns electron density from average A and Z
 * @param[in]   n                    size of arrays
 * @param[in]   density_g_cm3        material density in g/cm3 (array of size n)
 * @param[in]   average_A            average mass number (array of size n)
 * @param[in]   average_Z            average atomic number (array of size n)
 * @param[out]  electron_density_m3  electron density in 1/m3 (array of size n)
 */
void AT_electron_density_m3_multi( const long n,
    const double density_g_cm3[],
    const double average_Z[],
    const double average_A[],
    double electron_density_m3[]);

/**
 * Returns plasma energy needed for Sternheimer
 * computation of density effect in stopping power
 * @param[in]  electron_density_m3  electron density in 1/m3
 */
double AT_plasma_energy_J_single( const double electron_density_m3 );

/**
 * Computes the electron density for a given material composition
 *
 * @param[in]  n                     number of constituents in material
 * @param[in]  density_g_cm3         physical density (in g per cm3) of material
 * @param[in]  Z                     atomic numbers of constituents (array of size n)
 * @param[in]  A                     mass numbers of constituents (array of size n)
 * @param[in]  weight_fraction       relative fractions of weight of constituents (array of size n)
 * @param[out] electron_density_m3   electron density per m3
 */
void AT_electron_density_m3_from_composition( const long n,
    const double density_g_cm3,
    const long Z[],
    const long A[],
    const double weight_fraction[],
    double* electron_density_m3);


/**
 * Computes the average mass number for a given material composition
 *
 * @param[in]  n                     number of constituents in material
 * @param[in]  A                     mass numbers of constituents (array of size n)
 * @param[in]  weight_fraction       relative fractions of weight of constituents (array of size n)
 * @param[out] average_A             average A
 */
void AT_average_A_from_composition( const long n,
    const long A[],
    const double weight_fraction[],
    double* average_A);


/**
 * Computes the average atomic number for a given material composition
 *
 * @param[in]  n                     number of constituents in material
 * @param[in]  Z                     atomic numbers of constituents (array of size n)
 * @param[in]  weight_fraction       relative fractions of weight of constituents (array of size n)
 * @param[out] average_Z             average Z
 */
void AT_average_Z_from_composition( const long n,
    const long Z[],
    const double weight_fraction[],
    double* average_Z);

/**
 * Computes the effective atomic number for a given material composition
 *
 * @param[in]  n                     	number of constituents in material
 * @param[in]  Z                     	atomic numbers of constituents (array of size n)
 * @param[in]  weight_fraction       	relative fractions of weight of constituents (array of size n)
 * @param[in]  electron_densities_cm3   if not zero, weight fractions will additionally include electron densities per volume (array of size n)
 * @param[in]  exponent              	exponent for additivity rule reflecting the photon energy regime (usually 3.5 at ~ 100 kV)
 * @param[out] effective_Z           	effective Z
 */
void AT_effective_Z_from_composition( const long n,
    const long Z[],
    const double weight_fraction[],
    const double electron_densities_cm3[],
    const double exponent,
    double* effective_Z);

/**
 * Computes the I value for a given material composition
 *
 * @param[in]  n                     number of constituents in material
 * @param[in]  Z                     atomic numbers of constituents (array of size n)
 * @param[in]  A                     mass numbers of constituents (array of size n)
 * @param[in]  weight_fraction       relative fractions of weight of constituents (array of size n)
 * @param[out] I_eV                  I value in eV
 */
void AT_I_eV_from_composition( const long n,
    const long Z[],
    const long A[],
    const double weight_fraction[],
    double* I_eV);

/**
 * Initializes user defined material. The material can then be used with material index number 0. !Be aware! that is
 * definition is only valid during run-time. When the library is reloaded, the default settings are restored and the
 * material is removed.
 * @param[in] density_g_cm3  physical density in g per cm3
 * @param[in] I_eV           I value in eV
 * @param[in] average_A      average mass number
 * @param[in] average_Z      average atomic number
 * @param[out] status        material defined successfully if zero
 */
void AT_set_user_material( const double density_g_cm3,
		const double I_eV,
		const double average_A,
		const double average_Z,
		long* status);

/**
 * Initializes user defined material from composition data. The material can then be used with material index number 0. !Be aware! that is
 * definition is only valid during run-time. When the library is reloaded, the default settings are restored and the
 * material is removed.
 *
 * @param[in]  n                     number of constituents in material
 * @param[in]  density_g_cm3         physical density (in g per cm3) of material
 * @param[in]  Z                     atomic numbers of constituents (array of size n)
 * @param[in]  A                     mass numbers of constituents (array of size n)
 * @param[in]  weight_fraction       relative fractions of weight of constituents (array of size n)
 * @param[out] status                material defined successfully if zero
 */
void AT_set_user_material_from_composition( const long n,
		const double density_g_cm3,
		const long A[],
	    const long Z[],
	    const double weight_fraction[],
		long* status);
#endif

/* AT_DATAMATERIAL_H_ */
