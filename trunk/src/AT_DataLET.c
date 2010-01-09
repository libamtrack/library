/**
*    AT_DataLET.c
*    ==============
*
*    Created on: 09.01.2010
*    Author: kongruencja
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

#include "AT_DataLET.h"

// get LET-data for given material
void getPSTARvalue(
    const long* n,
    const float* x,
    const long* material_no,
    const float* x_table,
    const float* y_table,
    float* y)
{
  // first: find those PSTAR entries that match the material name
  bool*    matches    =  (bool*)calloc(AT_PSTAR_Data.n, sizeof(bool));
  matchi(    material_no,
      AT_PSTAR_Data.material_no,
      &AT_PSTAR_Data.n,
      matches);

  long    n_matches  = 0;
  long    i;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      n_matches++;
    }
    //    printf(debf,"idx: %i, match: %d\n",i, matches[i]);
  }

  // allocate vectors for extracted LET entries
  float*  x_c  =  (float*)calloc(n_matches, sizeof(float));
  float*  y_c  =  (float*)calloc(n_matches, sizeof(float));

  // and get the values
  long     j  = 0;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      x_c[j]  = x_table[i];
      y_c[j]  = y_table[i];
      j++;
    }
  }
  long  n_pol      = 4 + 1;
  for (i = 0; i < *n; i++){
    // Get proton-LET for scaled energy from table E, L using 4th degree polynomial (n_pol - 1 = 2) interpolation
    float  err_y_tmp  = 0.0f;    // dummy
    interp(    x_c,
        y_c,
        &n_matches,
        &n_pol,
        &x[i],
        &y[i],
        &err_y_tmp);
  }

  free(x_c);
  free(y_c);
  free(matches);
}


void AT_LET_MeV_cm2_g(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const long*  material_no,
    float*  LET_MeV_cm2_g)
{

  // get scaled energies for all given particles and energies
  float*  sE  =  (float*)calloc(*n, sizeof(float));
  AT_scaled_energy(  n,
      E_MeV_u,
      particle_no,
      sE);

  // get effective charge for all given particles and energies
  float*  eC  =  (float*)calloc(*n, sizeof(float));
  AT_effective_charge_from_particle_no(  n,
      E_MeV_u,
      particle_no,
      eC);

  getPSTARvalue(n, sE, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, LET_MeV_cm2_g);
  long   i;
  for (i = 0; i < *n; i++){
    LET_MeV_cm2_g[i] *=   eC[i] * eC[i];
  }

  free(eC);
  free(sE);
}


void AT_LET_keV_um(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const long*  material_no,
    float*  LET_keV_um)
{
  long  match, n_tmp = 1;
  pmatchi(  material_no,
      &n_tmp,
      AT_Material_Data.material_no,
      &AT_Material_Data.n,
      &match);

  // Get mass-norm. LET
  AT_LET_MeV_cm2_g(  n,
      E_MeV_u,
      particle_no,
      material_no,
      LET_keV_um);

  long  i;
  for (i = 0; i < *n; i++){
    LET_keV_um[i]  *=  AT_Material_Data.density_g_cm3[match] * 0.1f;
  }

}

void AT_CSDA_range_g_cm2(  const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const long*  material_no,
    float*  CSDA_range_g_cm2)
{
  getPSTARvalue(n, E_MeV_u, material_no, AT_PSTAR_Data.range_cdsa_g_cm2, AT_PSTAR_Data.kin_E_MeV, CSDA_range_g_cm2);
}

void AT_E_MeV_from_CDSA_range(  const long*  n,
    const float*  CSDA_range_g_cm2,
    const long*  particle_no,
    const long*  material_no,
    float*  E_MeV)
{
  // scaled energies
  float*  sE  =  (float*)calloc(*n, sizeof(float));

  getPSTARvalue(n, CSDA_range_g_cm2, material_no, AT_PSTAR_Data.range_cdsa_g_cm2, AT_PSTAR_Data.kin_E_MeV, sE);

  // scale energy back
  AT_E_MeV_u_from_scaled_energy(n , sE, particle_no, E_MeV);

  free( sE );
}

void AT_E_MeV_from_LET(  const long*  n,
    const float*  LET_MeV_cm2_g,
    const long*  particle_no,
    const long*  material_no,
    float*  E_MeV)
{

  //printf("n = %ld\n", *n);
  //printf("particle_no = %ld\n", *particle_no);
  //printf("material_no = %ld\n", *material_no);
  //printf("AT_E_MeV_from_LET LET = %g\n", *LET_MeV_cm2_g);

  // scaled energies
  float*  sE  =  (float*)calloc(*n, sizeof(float));

  //TODO add effective charge correction !!

  // Conversion LETion => LETproton
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  float*  charge  =  (float*)calloc(*n, sizeof(float));
  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  //printf("AT_E_MeV_from_LET LET ion = %g\n", LET_MeV_cm2_g[0]);

  // loop over n to find charge for all given particles and energies

  float*  LET_MeV_cm2_g_copy = (float*)calloc(*n,sizeof(float));
  memcpy(LET_MeV_cm2_g_copy,LET_MeV_cm2_g,*n);
  long  i;
  for(i = 0; i < *n; i++){
    charge[i]  = AT_Particle_Data.Z[matches[i]];
    LET_MeV_cm2_g_copy[i] /= (charge[i]*charge[i]);
  }

  //printf("AT_E_MeV_from_LET LET proton = %g\n", LET_MeV_cm2_g[0]);

  getPSTARvalue(n, LET_MeV_cm2_g_copy, material_no, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, AT_PSTAR_Data.kin_E_MeV, sE);

  // scale energy back
  AT_E_MeV_u_from_scaled_energy(n , sE, particle_no, E_MeV);

  free( LET_MeV_cm2_g_copy );
  free( sE );
}


