#ifndef AT_TRANSPORT_H_
#define AT_TRANSPORT_H_

/**
 *    AT_Transport.h
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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "AT_Constants.h"
#include "AT_NumericalRoutines.h"
#include "AT_DataMaterial.h"

extern int indent_counter;
extern char isp[];
extern FILE * debf;

void AT_BortfeldTransportProton(  float*  E_initial_MeV,
                  float*  sE_initial_MeV,
                  float*  fluence_initial_cm2,
                  char*  plateau_dose_material_name,
                  long*  n_shielding_slabs,
                  float*  shielding_thickness_m,
                  char**  shielding_material_name,
                  long*  n_detector_slabs,
                  float*  detector_thickness_m,
                  char*  detector_material_name,
                  /* return values */
                  float*  plateau_dose_Gy,
                  float*  plateau_dose_noNuc_Gy,
                  float*  detector_z_cm,
                  float*  detector_E_MeV,
                  float*  detector_sE_MeV,
                  float*  detector_fluence_cm2,
                  float*  detector_dfluencedz_cm,
                  float*  detector_LET_MeV_g_cm2,
                  float*  detector_dose_Gy,
                  float*  detector_dose_noNuc_Gy,
                  float*  geom_factor);


void AT_BortfeldTransportProtonL(  float*  E_initial_MeV,
                  float*  sE_initial_MeV,
                  float*  fluence_initial_cm2,
                  char*  plateau_dose_material_name,
                  long*  n_shielding_slabs,
                  float*  shielding_thickness_m,
                  char*  shielding_material_nameL,
                  long*  n_detector_slabs,
                  float*  detector_thickness_m,
                  char*  detector_material_name,
                  /* return values */
                  float*  plateau_dose_Gy,
                  float*  plateau_dose_noNuc_Gy,
                  float*  detector_z_cm,
                  float*  detector_E_MeV,
                  float*  detector_sE_MeV,
                  float*  detector_fluence_cm2,
                  float*  detector_dfluencedz_cm,
                  float*  detector_LET_MeV_g_cm2,
                  float*  detector_dose_Gy,
                  float*  detector_dose_noNuc_Gy,
                  float*  geom_factor);


void AT_BortfeldTransportProtonS(  float*  E_initial_MeV,
                  float*  sE_initial_MeV,
                  float*  fluence_initial_cm2,
                  char**  plateau_dose_material_name,
                  long*  n_shielding_slabs,
                  float*  shielding_thickness_m,
                  char**  shielding_material_name,
                  long*  n_detector_slabs,
                  float*  detector_thickness_m,
                  char**  detector_material_name,
                  /* return values */
                  float*  plateau_dose_Gy,
                  float*  plateau_dose_noNuc_Gy,
                  float*  detector_z_cm,
                  float*  detector_E_MeV,
                  float*  detector_sE_MeV,
                  float*  detector_fluence_cm2,
                  float*  detector_dfluencedz_cm,
                  float*  detector_LET_MeV_g_cm2,
                  float*  detector_dose_Gy,
                  float*  detector_dose_noNuc_Gy,
                  float*  geom_factor);

void  AT_Dose_Gy(  float*  density_g_cm3,
            float*  fluence_cm2,
            float*  LET_MeV_g_cm2,
            float*  dfluencedz_cm,
            float*  E_MeV,
            float*  dose_Gy,
            float*  dose_noNuc_Gy);

#endif // AT_TRANSPORT_H_
