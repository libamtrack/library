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

enum material_no{
  Water_Liquid         = 1,
  Aluminum_Oxide       = 2,
  Aluminum             = 3,
  PMMA                 = 4,
  Alanine              = 5
};

#define MATERIAL_DATA_N    5

typedef struct {
  const long    n;
  const long    material_no[MATERIAL_DATA_N];
  const float   density_g_cm3[MATERIAL_DATA_N];
  const float   electron_density_m3[MATERIAL_DATA_N];
  const float   I_eV[MATERIAL_DATA_N];
  const float   alpha_g_cm2_MeV[MATERIAL_DATA_N];
  const float   p_MeV[MATERIAL_DATA_N];
  const float   m_g_cm2[MATERIAL_DATA_N];
  const float   average_A[MATERIAL_DATA_N];
  const float   average_Z[MATERIAL_DATA_N];
  const char*   material_name[MATERIAL_DATA_N];
} material_data;

// TODO: replace electron density by conversion routine using A, Z
static const material_data AT_Material_Data = {
    MATERIAL_DATA_N,
    {  Water_Liquid,     Aluminum_Oxide,   Aluminum,      PMMA,       Alanine},
    {  1.00f,            3.97f,            2.6989f,       1.19f,      1.42f},
    {  3.3456e29f,       1.1719e30f,       7.8314e+29f,   3.8698e29f, 4.60571e29f},
    {  75.0f,            145.2f,           166.0f,        74.0f,      71.9f},
    {  0.00231f,         0.003058f,        0.003266f,     0.001988f,  0.00216381f},
    {  1.761f,           1.748f,           1.745f,        1.762f,     1.79165987f},
    {  0.01153f,         0.01305f,         0.01230f,      0.01338f,   -100.0f},        // No data processed for nuclear interactions in Alanine, hence set to -100
    {  13.0f,            0.0f,             27.0f,         0.0f,       0.0f},           //TODO find average A values
    {  7.22f,            0.0f,             13.0f,         0.0f,       0.0f},           //TODO find average Z values
    {  "Water, Liquid", "Aluminum Oxide",  "Aluminum",    "PMMA",     "Alanine"     }
};

//TODO Cucinnotta calculated average A for water as 14.3, but it seems that it is 13.0 (Leszek)

void getMaterialName( const long* material_no,
    char* material_name);

void getMaterialNo( const char* material_name,
    long* material_no);


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
void AT_getMaterialData( const long*  n,
    const long*  material_no,
    float*  density_g_cm3,
    float*  electron_density_m3,
    float*  I_eV,
    float*  alpha_g_cm2_MeV,
    float*  p_MeV,
    float*  m_g_cm2,
    float*  average_A,
    float*  average_Z
    );

void AT_density_g_cm3_from_material_no( const long*  n,
    const long*  material_no,
    float*       density_g_cm3);

#endif /* AT_DATAMATERIAL_H_ */
