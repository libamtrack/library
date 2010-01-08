#ifndef AT_FUNCTIONS_H_
#define AT_FUNCTIONS_H_

/**
 *    AT_Function.s
 *    =============
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

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "AT_Constants.h"
#include "AT_Data.h"
#include "AT_Utils.h"
#include "AT_ElectronRange.h"
#include "AT_PhysicsRoutines.h"

extern int indent_counter;
extern char isp[];
extern FILE * debf;

///////////////////////////////////////////////////////////////////////
// DATA ACCESS ROUTINE EXPORT (MAINLY FOR DEBUGGING)
void AT_LET_MeV_cm2_g(  long*  n,
            float*  E_MeV_u,
            long*  particle_no,
            long*  material_no,
            float*  LET_MeV_cm2_g);

void AT_LET_keV_um(  long*  n,
            float*  E_MeV_u,
            long*  particle_no,
            long*  material_no,
            float*  LET_keV_um);

void AT_CSDA_range_g_cm2(  long*  n,
              float*  E_MeV_u,
              long*  particle_no,
              long*  material_no,
              float*  CSDA_range_g_cm2);


void AT_Particle_Properties(  long*  particle_no,
                /* return values*/
                char**  particle_name,
                char**  USRTRACK_name,
                char**  element_name,
                long*  Z,
                long*  A,
                float*  mass);

void AT_getMaterialData(    long*  n,
                long*  material_no,
                float*  density_g_cm3,
                float*  electron_density_m3,
                float*  I_eV,
                float*  alpha_g_cm2_MeV,
                float*  p_MeV,
                float*  m_g_cm2);

void AT_E_MeV_from_LET(  long*  n,
              float*  LET_MeV_cm2_g,
              long*  particle_no,
              long*  material_no,
              float*  E_MeV);

///////////////////////////////////////////////////////////////////////
// MISC CONV. ROUTINES

#endif // AT_FUNCTIONS_H_
