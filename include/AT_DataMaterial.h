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
* Returns material data
* @param  n  number of materials the routine is called for (pointer to single variable)
* @param  material_no  material indices (pointer to long array of length n)
* @param  density_g_cm3  material density in g/cm3 (pointer to float array of length n)
* @param  electron_density_m3  electron density in 1/m3 (pointer to float array of length n)
* @param  I_eV  mean ionization potential in eV (pointer to float array of length n)
* @param  alpha_g_cm2_MeV  fit parameter for power-law representation of stp.power/range/E-dependence (pointer to float array of length n)
* @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
*       Here, however, we use the mass stopping power. The correct dimension is g/(cm^2 * MeV^p)
* @param  p_MeV  fit parameter for power-law representation of stp.power/range/E-dependence (pointer to float array of length n)
* @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
*       p is actually dimenionless, it should be nevertheless indicated, that the energy must be given in MeV
* @param  m_g_cm2  fit parameter for the linear representation of fluence changes due to nuclear interactions based on data from Janni, 1982 (pointer to float array of length n)
* @see  Bortfeld, T. (1997), An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24, 2024ff.
* @param  average_A  average mass number (pointer to float array of length n)
* let f_i be fraction by weight of the constituent element with atomic number Z_i and atomic weight A_i<BR>
* let us define average_Z/A = \sum_i f_i Z_i / A_i <BR>
* then we have: average_A = average_Z / (average_Z/A)
* for water (H20) we have: average_Z/A = (2/18) * (1/1) + (16/18)*(8/16) = 0.5555
* average_A = 7.22 / 0.555 = 13
* @param  average_Z  average atomic number (pointer to float array of length n)<BR>
* let f_i be fraction by weight of the constituent element with atomic number Z_i and atomic weight A_i<BR>
* average_Z = \sum_i f_i Z_i <BR>
* for water (H20) we have: average_Z = (2/18)*1 + (16/18)*8 = 7.22
* @see Tabata, T. (1972) Generalized semiempirical equations for the extrapolated range of electrons, Nucl. Instr and Meth. 103, 85-91.
*/
void AT_getMaterialData( const long  n,
    const long*  material_no,
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
 */
void AT_density_g_cm3_from_material_no( const long  n,
    const long*  material_no,
    double*       density_g_cm3);

#endif /* AT_DATAMATERIAL_H_ */
