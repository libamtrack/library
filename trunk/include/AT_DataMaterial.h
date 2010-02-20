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

static const material_data AT_Material_Data = {
    MATERIAL_DATA_N,
    {  Water_Liquid,     Aluminum_Oxide,   Aluminum,      PMMA,       Alanine},
    {  1.00f,            3.97f,            2.6989f,       1.19f,      1.42f},
    {  3.3456e29f,       1.1719e30f,       7.8314e+29f,   3.8698e29f, 4.60571e29f},
    {  75.0f,            145.2f,           166.0f,        74.0f,      71.9f},
    {  0.00231f,         0.003058f,        0.003266f,     0.001988f,  0.00216381f},
    {  1.761f,           1.748f,           1.745f,        1.762f,     1.79165987f},
    {  0.01153f,         0.01305f,         0.01230f,      0.01338f,   -100.0f},
    {  14.3f,            0.0f,             0.0f,          0.0f,       0.0f},           //TODO find average A values
    {  7.22f,            0.0f,             0.0f,          0.0f,       0.0f},           //TODO find average Z values
    {  "Water, Liquid", "Aluminum Oxide",  "Aluminum",    "PMMA",     "Alanine"     }
};


void getMaterialName( const long* material_no,
    char* material_name);

void getMaterialNo( const char* material_name,
    long* material_no);

void AT_getMaterialData( const long*  n,
    const long*  material_no,
    float*  density_g_cm3,
    float*  electron_density_m3,
    float*  I_eV,                  //TODO describe this !
    float*  alpha_g_cm2_MeV,       //TODO describe this !
    float*  p_MeV,                 //TODO describe this !
    float*  m_g_cm2,               //TODO describe this !
    float*  average_A,
    float*  average_Z
    );

void AT_density_g_cm3_from_material_no( const long*  n,
    const long*  material_no,
    float*       density_g_cm3);

#endif /* AT_DATAMATERIAL_H_ */
