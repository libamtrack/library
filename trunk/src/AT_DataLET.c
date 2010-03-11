/**
 * @file
 * @brief LET tables
 */

/*
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


void getPSTARvalue(
    const long   n,
    const float  x[],
    const long   material_no,
    const float  x_table[],
    const float  y_table[],
    float        y[])
{
  // first: find those PSTAR entries that match the material number
  bool*    matches    =  (bool*)calloc(AT_PSTAR_Data.n, sizeof(bool));
  matchi(    &material_no,
      AT_PSTAR_Data.material_no,
      &AT_PSTAR_Data.n,
      matches);

  long    n_matches  = 0;
  long    i;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      n_matches++;
    }
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
  for (i = 0; i < n; i++){
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


void AT_LET_MeV_cm2_g(  const long  number_of_particles,
    const float  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    float        LET_MeV_cm2_g[])
{

  // TODO think if we shall use scaled energy or just energy per nucleon
  // get scaled energies for all given particles and energies
  float*  sE  =  (float*)calloc(number_of_particles, sizeof(float));
  AT_scaled_energy(  &number_of_particles,
      E_MeV_u,
      particle_no,
      sE);

  // get effective charge for all given particles and energies
  float*  Zeff_proton  =  (float*)calloc(number_of_particles, sizeof(float));
  long*   particle_no_proton  =  (long*)calloc(number_of_particles, sizeof(long));
  long   i;
  for (i = 0; i < number_of_particles; i++){
    particle_no_proton[i] = 1;
  }
  AT_effective_charge_from_particle_no(  &number_of_particles,
      E_MeV_u,
      particle_no_proton,
      Zeff_proton);

  // get effective charge for all given particles and energies
  float*  Zeff_ion  =  (float*)calloc( number_of_particles, sizeof(float));
  AT_effective_charge_from_particle_no(  &number_of_particles,
      E_MeV_u,
      particle_no,
      Zeff_ion);

  getPSTARvalue(number_of_particles, sE, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, LET_MeV_cm2_g);
  for (i = 0; i < number_of_particles; i++){
    if( particle_no[i] != 1){ // for particles other than proton scale LET by (Zeff_ion / Zeff_proton)^2
      LET_MeV_cm2_g[i] *=   gsl_pow_2(Zeff_ion[i] / Zeff_proton[i]);
    }
  }

  free(Zeff_ion);
  free(Zeff_proton);
  free(sE);
}


void AT_LET_keV_um(  const long  number_of_particles,
    const float  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    float        LET_keV_um[])
{
  // Get material density
  double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);

  // Get mass-norm. LET
  AT_LET_MeV_cm2_g(  number_of_particles,
      E_MeV_u,
      particle_no,
      material_no,
      LET_keV_um);

  long  i;
  for (i = 0; i < number_of_particles; i++){
    LET_keV_um[i]  *=  (float)material_density_g_cm3 * 0.1f;
  }

}


void AT_CSDA_range_g_cm2(  const long  number_of_particles,
    const float   E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    float         CSDA_range_g_cm2[])
{
  getPSTARvalue(number_of_particles, E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.range_cdsa_g_cm2, CSDA_range_g_cm2);

  // Conversion CSDA_proton => CSDA_ion
  long*  Z  =  (long*)calloc(number_of_particles, sizeof(long));
  long*  A  =  (long*)calloc(number_of_particles, sizeof(long));
  AT_Z_from_particle_no(        &number_of_particles,
                                particle_no,
                                Z);
  AT_A_from_particle_no(        &number_of_particles,
                                particle_no,
                                A);
  long i = 0;
  for (i = 0; i < number_of_particles; i++){
    if (particle_no[i] != 1){
      CSDA_range_g_cm2[i]  *=   (float)(A[i])/(float)((Z[i])*(Z[i]));
    }
  }

  free(Z);
  free(A);
}


void AT_CSDA_range_m(  const long  number_of_particles,
    const float  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    float        CSDA_range_m[])
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
    CSDA_range_m[i]  /=  (float)material_density_g_cm3 * 100.0f;
  }

}


void AT_E_MeV_from_CDSA_range(  const long  number_of_particles,
    const float  CSDA_range_g_cm2[],
    const long   particle_no[],
    const long   material_no,
    float        E_MeV[])
{
  // scaled energies
  float*  sE  =  (float*)calloc(number_of_particles, sizeof(float));

  getPSTARvalue(number_of_particles, CSDA_range_g_cm2, material_no, AT_PSTAR_Data.range_cdsa_g_cm2, AT_PSTAR_Data.kin_E_MeV, sE);

  // scale energy back
  AT_E_MeV_u_from_scaled_energy(&number_of_particles , sE, particle_no, E_MeV);

  free( sE );
}


void AT_E_MeV_from_LET(  const long  number_of_particles,
    const float  LET_MeV_cm2_g[],
    const long   particle_no[],
    const long   material_no,
    float        E_MeV[])
{
  // scaled energies
  float*  scaled_E_MeV  =  (float*)calloc(number_of_particles, sizeof(float));

  //TODO add effective charge correction !!

  // Conversion LETion => LETproton
  long*  charge  =  (long*)calloc(number_of_particles, sizeof(long));
  AT_Z_from_particle_no(&number_of_particles,particle_no,charge);


  // loop over n to find charge for all given particles and energies
  float*  LET_MeV_cm2_g_copy = (float*)calloc(number_of_particles,sizeof(float));
  memcpy(LET_MeV_cm2_g_copy,LET_MeV_cm2_g,number_of_particles);
  long  i;
  for(i = 0; i < number_of_particles; i++){
    LET_MeV_cm2_g_copy[i] /= gsl_pow_2( (float)charge[i] );;
  }

  getPSTARvalue(number_of_particles, LET_MeV_cm2_g_copy, material_no, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, AT_PSTAR_Data.kin_E_MeV, scaled_E_MeV);

  // scale energy back
  AT_E_MeV_u_from_scaled_energy(&number_of_particles , scaled_E_MeV, particle_no, E_MeV);

  free( charge );
  free( LET_MeV_cm2_g_copy );
  free( scaled_E_MeV );
}
