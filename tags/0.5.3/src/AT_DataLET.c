/**
 * @brief LET tables
 */

/*
 *    AT_DataLET.c
 *    ==============
 *
 *    Created on: 09.01.2010
 *    Creator: kongruencja
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

#include "AT_DataLET.h"
#include "AT_DataParticle.h"
#include <math.h>

void get_table_value(
    const long    n,
    const double  x[],
    const long    subset_no,
    const double  x_table[],
    const double  y_table[],
    double        y[])
{
  // first: find those PSTAR entries that match the material number
  bool*    matches    =  (bool*)calloc(AT_PSTAR_Data.n, sizeof(bool));
  is_element_int(    subset_no,
      AT_PSTAR_Data.material_no,
      AT_PSTAR_Data.n,
      matches);

  long    n_matches  = 0;
  long    i;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      n_matches++;
    }
  }

  // allocate vectors for extracted LET entries
  double*  x_c  =  (double*)calloc(n_matches, sizeof(double));
  double*  y_c  =  (double*)calloc(n_matches, sizeof(double));

  // and get the values
  long     j  = 0;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      x_c[j]  = x_table[i];
      y_c[j]  = y_table[i];
      j++;
    }
  }
  for (i = 0; i < n; i++){
    // Get proton-LET for scaled energy from table E, L using linear interpolation (to be done on logscale TODO)
    y[i] = AT_get_interpolated_y_from_input_table(x_c, y_c, n_matches, x[i]);
  }

  free(x_c);
  free(y_c);
  free(matches);
}

double AT_CSDA_range_g_cm2_single(   const double  E_MeV_u,
    const long    particle_no,
    const long    material_no)
{
  double CSDA_range_g_cm2;
  const long number_of_particles  =  1;
  get_table_value(number_of_particles, &E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.range_cdsa_g_cm2, &CSDA_range_g_cm2);

  // Conversion CSDA_proton => CSDA_ion
  long Z = AT_Z_from_particle_no_single(particle_no);
  long A = AT_A_from_particle_no_single(particle_no);

  return CSDA_range_g_cm2  * (double)(A)/(double)(Z*Z);
}



void AT_CSDA_range_g_cm2(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    double        CSDA_range_g_cm2[])
{
  get_table_value(number_of_particles, E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.range_cdsa_g_cm2, CSDA_range_g_cm2);

  // Conversion CSDA_proton => CSDA_ion
  long*  Z  =  (long*)calloc(number_of_particles, sizeof(long));
  long*  A  =  (long*)calloc(number_of_particles, sizeof(long));
  AT_Z_from_particle_no(        number_of_particles,
                                particle_no,
                                Z);
  AT_A_from_particle_no(        number_of_particles,
                                particle_no,
                                A);
  long i = 0;
  for (i = 0; i < number_of_particles; i++){
    if (particle_no[i] != 1){
      CSDA_range_g_cm2[i]  *=   (double)(A[i])/(double)((Z[i])*(Z[i]));
    }
  }

  free(Z);
  free(A);
}


double AT_CSDA_range_m_single(  const double  E_MeV_u,
		const long    particle_no,
		const long    material_no){

	double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);
	double CSDA_range_g_cm2 = AT_CSDA_range_g_cm2_single(E_MeV_u, particle_no, material_no);

	return CSDA_range_g_cm2 / (material_density_g_cm3 * 100.0);
}


void AT_CSDA_range_m(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    double        CSDA_range_m[])
{
  // Get material density
  double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);

  // Get mass-norm. CSDA range
  AT_CSDA_range_g_cm2(  number_of_particles,
      E_MeV_u,
      particle_no,
      material_no,
      CSDA_range_m);

  long  i;
  for (i = 0; i < number_of_particles; i++){
    CSDA_range_m[i]  /=  material_density_g_cm3 * 100.0;
  }

}


