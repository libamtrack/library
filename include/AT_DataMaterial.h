#ifndef AT_DATAMATERIAL_H_
#define AT_DATAMATERIAL_H_

/**
 * @file
 * @brief Material properties
 */

/*
 *    AT_DataMaterial.h
 *    ==================
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

#include <string.h>
#include <stdio.h>

#include "AT_Error.h"
#include "AT_NumericalRoutines.h"

/**
 * Materials code numbers
 */
enum material_no{
  User_Defined_Material= 0, /**< To be defined by the user during runtime >**/
  Water_Liquid         = 1, /**< Liquid water */
  Aluminum_Oxide       = 2, /**< Aluminium oxide */
  Aluminum             = 3, /**< Aluminium */
  PMMA                 = 4, /**< PMMA */
  Alanine              = 5, /**< Alanine */
  LiF                  = 6  /**< Lithium Fluoride */
};


#define MATERIAL_DATA_N    7

// The next two LET-related structures must be declared here rather than in AT_DataLET.h to avoid circular dependencies
/**
 * Stopping power data for a material and a particle type
 */
typedef struct AT_LET_data_single{
  long     n;                                     /**< number of items in the data table */
  long     particle_no;                           /**< particle number this data is for, see AT_DataParticle.h for definition */
  double*  kin_E_MeV;                             /**< Kinetic energy in MeV, pointer to array of size n */
  double*  stp_pow_el_MeV_cm2_g;                  /**< Electronic (Collision) Stopping Power, pointer to array of size n */
  double*  range_cdsa_g_cm2;                      /**< CSDA (continuous-slowing-down approximation) range, pointer to array of size n */
} AT_LET_data_single;

/**
 * Stopping power data for a material
 */
typedef struct AT_LET_data{
  long            n;                              /**< number of data tables in the structure */
  AT_LET_data_single*   LET_data_single;                /**< pointer to LET data tables for a particle type */
} AT_LET_data;

#define MATERIAL_NAME_LENGTH    255

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
  double          electron_density_m3;                    /**< Electron density (in 1/m^3),  will be calculated if not given */

  long            LET_data_source;                        /**< Defines the source for stopping power data, see enum LET_data_source [ESSENTIAL]. The user has to make sure that he provides the necessary data. */

  double          p_MeV;                                  /**< Prefactor for the power-law description of stopping power: S = p*E^alpha. In MeV^(1/alpha) */
  double          alpha_g_cm2_MeV;                        /**< Exponent for the power-law description of stopping */

  AT_LET_data     LET_data;                               /**< Stopping power data for the material, see AT_DataLET.h for definition. Will be calculated in case of power-law or given by the user (or read in from PSTAR data?). */
                                                          /**<  LET_data has to hold at_least the proton stopping powers. For any particle where no explicit stopping power data is given, a scaling by Z_eff^2 will be performed. */

//  const double  m_g_cm2;                                /**< Slope of linear approximation in fluence reduction due to nuclear interactions, not used yet */
} AT_material;


typedef struct {
  const long    n;
  const long    material_no[MATERIAL_DATA_N];
  const bool    ready[MATERIAL_DATA_N];
  const long    ICRU_ID[MATERIAL_DATA_N];
  const double  density_g_cm3[MATERIAL_DATA_N];
  const double  electron_density_m3[MATERIAL_DATA_N];
  const double  I_eV[MATERIAL_DATA_N];
  const double  alpha_g_cm2_MeV[MATERIAL_DATA_N];
  const double  p_MeV[MATERIAL_DATA_N];
  const double  m_g_cm2[MATERIAL_DATA_N];
  const double  average_A[MATERIAL_DATA_N];
  const double  average_Z[MATERIAL_DATA_N];
  const char*   material_name[MATERIAL_DATA_N];
} material_data;


static const material_data AT_Material_Data = {
    MATERIAL_DATA_N,
    {  User_Defined_Material,    Water_Liquid,   Aluminum_Oxide,   Aluminum,     PMMA,      Alanine,     LiF},           // material_no
    {  false,           true,           true,             true,         true,      true,        true},          // ready
    {  0,               276,             106,              13,           223,       0,           185},          // ICRU_ID
    {  0.0,             1.00,            3.97,             2.6989,       1.19,      1.42,        2.64},         // density_g_cm3
    {  0.0,             3.3456e29,       1.1719e30,        7.8314e29,    3.8698e29, 4.60571e29,  7.341e29},     // electron_density_g_cm3
    {  0.0,             75.0,            145.2,            166.0,        74.0,      71.9,        10.0},         // I_eV
    {  0.0,             0.00231,         0.003058,         0.003266,     0.001988,  0.00216381,  0.0},          // alpha_g_cm2 - TODO No data for LiF
    {  0.0,             1.761,           1.748,            1.745,        1.762,     1.79165987,  0.0},          // p_MeV - TODO No data for LiF
    {  0.0,             0.01153,         0.01305,          0.01230,      0.01338,   -100.0,      0.0},          // m_g_cm2 - TODO No data processed for nuclear interactions in Alanine, hence set to -100
    {  0.0,             13.0,            21.72,            27.0,         0.0,       0.0,         17.7333},      // average_A - TODO find average A values for PMMA and alanine
    {  0.0,             7.22,            10.637,           13.0,         0.0,       0.0,         8.0},          // average_Z - TODO find average Z values for PMMA and alanine
    {  "User defined",  "Water, Liquid", "Aluminum Oxide", "Aluminum",   "PMMA",    "Alanine",   "Lithium Fluoride" }
};

/* Cucinnotta calculated average A for water as 14.3, but it seems that it is 13.0 (Leszek) *
 * looks like error in his article */

/**
 * Get index of material in AT_Material_Data for given material_no
 * (currently for example material with number 2 has index 1)
 *
 * @param material_number  material number
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
 * Get electron density [1/m3] for single material with number material_no
 * @param[in] material_no
 * @return    electron density [1/m3]
 */
double AT_electron_density_m3_from_material_no( const long   material_no );


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
 * Returns material data for single material
 * @param[in]  material_no
 * @param[out] density_g_cm3
 * @param[out] electron_density_m3
 * @param[out] I_eV
 * @param[out] alpha_g_cm2_MeV
 * @param[out] p_MeV
 * @param[out] m_g_cm2
 * @param[out] average_A
 * @param[out] average_Z
 */
void AT_get_material_data(     const long  material_no,
    double*  density_g_cm3,
    double*  electron_density_m3,
    double*  I_eV,
    double*  alpha_g_cm2_MeV,
    double*  p_MeV,
    double*  m_g_cm2,
    double*  average_A,
    double*  average_Z);


/**
 * Returns material data for list of materials
 * @param[in]   number_of_materials  numbers of materials the routine is called for (array of length number_of_materials)
 * @param[in]   material_no          material indices (array of length number_of_materials)
 * @param[out]  density_g_cm3        material density in g/cm3 (array of length number_of_materials)
 * @param[out]  electron_density_m3  electron density in 1/m3 (array of length number_of_materials)\n
 *   electron_density = number_of_electron_per_molecule * avogadro_constant / molar_mass * 1m3 * density
 * @param[out]  I_eV                 mean ionization potential in eV (array of length number_of_materials)
 * @param[out]  alpha_g_cm2_MeV      fit parameter for power-law representation of stp.power/range/E-dependence (array of length nnumber_of_materials)
 *   @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
 *   Here, however, we use the mass stopping power. The correct dimension is g/(cm^2 * MeV^p)
 * @param[out]  p_MeV                fit parameter for power-law representation of stp.power/range/E-dependence (array of length number_of_materials)
 *   @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.\n
 *   p is actually dimensionless, it should be nevertheless indicated, that the energy must be given in MeV
 * @param[out]  m_g_cm2              fit parameter for the linear representation of fluence changes due to nuclear interactions based on data from Janni, 1982 (array of length number_of_materials)
 *   @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
 * @param[out]  average_A            average mass number (array of length number_of_materials) \n
 *   let f_i be fraction by weight of the constituent element with atomic number Z_i and atomic weight A_i\n
 *   let us define average_Z/A = \sum_i f_i Z_i / A_i \n
 *   then we have: average_A = average_Z / (average_Z/A) \n
 *   for water (H20) we have: average_Z/A = (2/18) * (1/1) + (16/18)*(8/16) = 0.5555 \n
 *   average_A = 7.22 / 0.555 = 13
 * @param[out]  average_Z            average atomic number (array of length number_of_materials)\n
 *   let f_i be fraction by weight of the constituent element with atomic number Z_i and atomic weight A_i\n
 *   average_Z = \sum_i f_i Z_i \n
 *   for water (H20) we have: average_Z = (2/18)*1 + (16/18)*8 = 7.22\n
 *   @see Tabata, T. (1972) Generalized semiempirical equations for the extrapolated range of electrons, Nucl. Instr and Meth. 103, 85-91.
 */
void AT_get_materials_data( const long  number_of_materials,
    const long  material_no[],
    double  density_g_cm3[],
    double  electron_density_m3[],
    double  I_eV[],
    double  alpha_g_cm2_MeV[],
    double  p_MeV[],
    double  m_g_cm2[],
    double  average_A[],
    double  average_Z[]);

/////////////////////////////////////////////////////////
/* TEST FUNCTIONS FOR NEW MATERIAL / LET DATA HANDLING */
int AT_check_material( AT_material* material);
int AT_establish_material(AT_material* material);
int AT_free_material(AT_material* material);
double AT_electron_density_m3( const long n,
    const double density_g_cm3,
    const long Z[],
    const long A[],
    const double weight_fraction[]);
double AT_average_A( const long n,
    const long A[],
    const double weight_fraction[]);
double AT_average_Z( const long n,
    const long Z[],
    const double weight_fraction[]);
/* END OF TEST FUNCTIONS FOR NEW MATERIAL / LET DATA HANDLING */
////////////////////////////////////////////////////////////////


#endif /* AT_DATAMATERIAL_H_ */
