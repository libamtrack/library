#ifndef AT_ELECTRONRANGE_H_
#define AT_ELECTRONRANGE_H_

/**
 * @file
 * @brief Electron range models
 */

/*
*    AT_ElectronRange.h
*    ==============
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
#include "gsl/gsl_pow_int.h"
#include "gsl/gsl_sf_log.h"

#include "AT_NumericalRoutines.h"
#include "AT_DataMaterial.h"
#include "AT_DataParticle.h"
#include "AT_PhysicsRoutines.h"

///////////////////////////////////////////////////////////////////////
// ER DATA

enum ERModels{
  ER_Test                  = 1,
      ER_ButtsKatz         = 2,
      ER_Waligorski        = 3,
      ER_Geiss             = 4,
      ER_Scholz            = 5,
      ER_Edmund            = 6,
      ER_Tabata            = 7
};

#define ER_DATA_N    7

typedef struct {
  long    n;
  long    ER_no[ER_DATA_N];
  char*   ER_name[ER_DATA_N];
} er_data;

static const er_data AT_ER_Data = {
    ER_DATA_N,
    {  ER_Test,                 ER_ButtsKatz,                                  ER_Waligorski,            ER_Geiss,                        ER_Scholz,                           ER_Tabata },
    {  "simple test ER model",  "Butts & Katz' [Katz et al., 1972] ER model",  "Waligorski's ER model",  "Geiss' [Geiss, 1997] ER model", "ER_Scholz' [Scholz, 2001] ER model", "ER_Tabata [Tabata, 1972] ER model"}
};

/**
* Returns name of the electron model from index
*
* @param  ER_no    electron-range-model index
* @return Er_name  string containing the electron-range model name
*/
void  getERName(  const long* ER_no,
    char* ER_name);

/**
* Returns the maximum electron range (track diameter) in m
*
* @param  n  number of particles in the incident field
* @param  E_MeV_u  kinetic energy for particles in the given field (vector of length n)
* @param  particle_no  particle indices for given field (vector of length n)
* @param  material_no  index for detector material
* @param  er_model  index for electron-range model chosen
* @return  max_electron_range_m  electron range (track diameter) in m
*/
void AT_max_electron_range_m( const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const long*  material_no,
    const long*  er_model,
    float*  max_electron_range_m);

#endif /* AT_ELECTRONRANGE_H_ */
