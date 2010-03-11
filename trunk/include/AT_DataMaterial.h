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

#include "AT_NumericalRoutines.h"

/**
 * Materials code numbers
 */
enum material_no{
  Water_Liquid         = 1, /**< Liquid water */
  Aluminum_Oxide       = 2, /**< Aluminium oxide */
  Aluminum             = 3, /**< Aluminium */
  PMMA                 = 4, /**< PMMA */
  Alanine              = 5  /**< Alanine */
};

#define MATERIAL_DATA_N    5

typedef struct {
  const long    n;
  const long    material_no[MATERIAL_DATA_N];
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

// TODO: replace electron density by conversion routine using A, Z
static const material_data AT_Material_Data = {
    MATERIAL_DATA_N,
    {  Water_Liquid,    Aluminum_Oxide,  Aluminum,     PMMA,      Alanine},
    {  1.00,            3.97,            2.6989,       1.19,      1.42},
    {  3.3456e29,       1.1719e30,       7.8314e+29,   3.8698e29, 4.60571e29},
    {  75.0,            145.2,           166.0,        74.0,      71.9},
    {  0.00231,         0.003058,        0.003266,     0.001988,  0.00216381},
    {  1.761,           1.748,           1.745,        1.762,     1.79165987},
    {  0.01153,         0.01305,         0.01230,      0.01338,   -100.0},        // No data processed for nuclear interactions in Alanine, hence set to -100
    {  13.0,            0.0,             27.0,         0.0,       0.0},           //TODO find average A values
    {  7.22,            0.0,             13.0,         0.0,       0.0},           //TODO find average Z values
    {  "Water, Liquid", "Aluminum Oxide",  "Aluminum",    "PMMA",     "Alanine"     }
};

// Cucinnotta calculated average A for water as 14.3, but it seems that it is 13.0 (Leszek)
// looks like error in his article

/**
 * TODO
 * @param material_no
 * @return
 */
long AT_index_from_material_no( const long material_no );


/**
 * Get material name
 * @param[in] material_no
 * @param[out] material_name
 */
void getMaterialName( const long material_no,
    char* material_name);


/**
 * Get material number
 * @param[in] material_name
 * @return material number
 */
long getMaterialNo( const char* material_name );


/**
 *TODO
 * @param material_no
 * @return
 */
double AT_density_g_cm3_from_material_no( const long   material_no );


/**
 *TODO
 * @param material_no
 * @return
 */
double AT_electron_density_m3_from_material_no( const long   material_no );


/**
 *TODO
 * @param material_no
 * @return
 */
double AT_I_eV_from_material_no( const long   material_no );


/**
 *TODO
 * @param material_no
 * @return
 */
double AT_alpha_g_cm2_MeV_from_material_no( const long   material_no );

/**
 *TODO
 * @param material_no
 * @return
 */
double AT_p_MeV_from_material_no( const long   material_no );


/**
 *TODO
 * @param material_no
 * @return
 */
double AT_m_g_cm2_from_material_no( const long   material_no );

/**
 *TODO
 * @param material_no
 * @return
 */
double AT_average_A_from_material_no( const long   material_no );


/**
 *TODO
 * @param material_no
 * @return
 */
double AT_average_Z_from_material_no( const long   material_no );


/**
 * TODO
 * @param material_no
 * @param density_g_cm3
 * @param electron_density_m3
 * @param I_eV
 * @param alpha_g_cm2_MeV
 * @param p_MeV
 * @param m_g_cm2
 * @param average_A
 * @param average_Z
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
 * TODO
 * @param number_of_materials
 * @param material_no
 * @param density_g_cm3
 */
void AT_densities_g_cm3_from_material_numbers( const long  number_of_materials,
    const long  material_no[],
    double      density_g_cm3[]);


/**
* Returns material data
* @param  number_of_materials[in]   numbers of materials the routine is called for (array of length number_of_materials)
* @param  material_no[in]           material indices (array of length number_of_materials)
* @param  density_g_cm3[out]        material density in g/cm3 (array of length number_of_materials)
* @param  electron_density_m3[out]  electron density in 1/m3 (array of length number_of_materials)
* @param  I_eV[out]                 mean ionization potential in eV (array of length number_of_materials)
* @param  alpha_g_cm2_MeV[out]      fit parameter for power-law representation of stp.power/range/E-dependence (array of length nnumber_of_materials)
* @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
*       Here, however, we use the mass stopping power. The correct dimension is g/(cm^2 * MeV^p)
* @param  p_MeV[out]                fit parameter for power-law representation of stp.power/range/E-dependence (array of length number_of_materials)
* @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
*       p is actually dimensionless, it should be nevertheless indicated, that the energy must be given in MeV
* @param  m_g_cm2[out]              fit parameter for the linear representation of fluence changes due to nuclear interactions based on data from Janni, 1982 (array of length number_of_materials)
* @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
* @param  average_A[out]            average mass number (array of length number_of_materials) \n
* let f_i be fraction by weight of the constituent element with atomic number Z_i and atomic weight A_i\n
* let us define average_Z/A = \sum_i f_i Z_i / A_i \n
* then we have: average_A = average_Z / (average_Z/A) \n
* for water (H20) we have: average_Z/A = (2/18) * (1/1) + (16/18)*(8/16) = 0.5555 \n
* average_A = 7.22 / 0.555 = 13
* @param  average_Z[out]            average atomic number (pointer to float array of length number_of_materials)\n
* let f_i be fraction by weight of the constituent element with atomic number Z_i and atomic weight A_i\n
* average_Z = \sum_i f_i Z_i \n
* for water (H20) we have: average_Z = (2/18)*1 + (16/18)*8 = 7.22
* @see Tabata, T. (1972) Generalized semiempirical equations for the extrapolated range of electrons, Nucl. Instr and Meth. 103, 85-91.
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

#endif /* AT_DATAMATERIAL_H_ */
