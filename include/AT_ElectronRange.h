#ifndef AT_ELECTRONRANGE_H_
#define AT_ELECTRONRANGE_H_

/**
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

#include "AT_NumericalRoutines.h"
#include "AT_DataMaterial.h"
#include "AT_DataParticle.h"

///////////////////////////////////////////////////////////////////////
// ER DATA

enum ERModels{
  ER_Test                  = 1,
      ER_ButtsKatz         = 2,
      ER_Waligorski        = 3,
      ER_Geiss             = 4,
      ER_Scholz            = 5,
      ER_Edmund                    = 6
};

#define ER_DATA_N    5

typedef struct {
  long    n;
  long    ER_no[ER_DATA_N];
  char*   ER_name[ER_DATA_N];
} er_data;

static const er_data AT_ER_Data = {
    ER_DATA_N,
    {  ER_Test,          ER_ButtsKatz,          ER_Waligorski,        ER_Geiss, ER_Scholz},
    {  "simple test ER model",  "Butts & Katz' [Katz et al., 1972] ER model",  "Waligorski's ER model",    "Geiss' [Geiss, 1997] ER model", "ER_Scholz' [Scholz, 2001] ER model"}
};

void   getERName(  long* ER_no, char* ER_name);


void AT_max_electron_range_m(  long*  n,
                float*  E_MeV_u,
                long*  particle_no,
                long*  material_no,
                long*  er_model,
                float*  max_electron_range_m);

#ifdef _R
void AT_max_electron_range_mS(  int*  n,
                float*  E_MeV_u,
                int*  particle_no,
                int*  material_no,
                int*   er_model,
                float*  max_electron_range_m);
#endif
#ifdef _S
void AT_max_electron_range_mS(  long*  n,
                float*  E_MeV_u,
                long*  particle_no,
                long*  material_no,
                long*  er_model,
                float*  max_electron_range_m);
#endif

#endif
