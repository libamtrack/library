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
    const double  x[],
    const long   material_no,
    const double  x_table[],
    const double  y_table[],
    double        y[])
{
  // first: find those PSTAR entries that match the material number
  bool*    matches    =  (bool*)calloc(AT_PSTAR_Data.n, sizeof(bool));
  is_element_int(    material_no,
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
  long  n_pol      = 4 + 1;
  for (i = 0; i < n; i++){
    // Get proton-LET for scaled energy from table E, L using 4th degree polynomial (n_pol - 1 = 2) interpolation
    double  err_y_tmp  = 0.0;    // dummy
    interp(    x_c,
        y_c,
        n_matches,
        n_pol,
        x[i],
        &y[i],
        &err_y_tmp);
  }

  free(x_c);
  free(y_c);
  free(matches);
}


double AT_LET_MeV_cm2_g_single(  const double  E_MeV_u,
    const long    particle_no,
    const long    material_no){

  // get LET for proton of same energy / nucleon
  const long number_of_particles  =  1;
  double LET_MeV_cm2_g            =  0.0;
  getPSTARvalue(number_of_particles, &E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, &LET_MeV_cm2_g);

  double Zeff_ion    =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u,particle_no);

  double Zeff_proton =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u,1001);

  // scale proton LET by ratio of effective Z
  if( particle_no != 1001){ // for particles other than proton scale LET by (Zeff_ion / Zeff_proton)^2
      LET_MeV_cm2_g *=   gsl_pow_2(Zeff_ion / Zeff_proton);
  }

  return LET_MeV_cm2_g;
}


void AT_LET_MeV_cm2_g(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    double        LET_MeV_cm2_g[])
{
  // get LET for proton of same energy / nucleon
  getPSTARvalue(number_of_particles, E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, LET_MeV_cm2_g);

  // get effective charge for all given particles and energies
  double*  Zeff_ion  =  (double*)calloc( number_of_particles, sizeof(double));
  AT_effective_charge_from_E_MeV_u(  number_of_particles,
      E_MeV_u,
      particle_no,
      Zeff_ion);

  // get effective charge for protons of same energy
  long   i;
  double*  Zeff_proton           =  (double*)calloc(number_of_particles, sizeof(double));
  long*    particle_no_proton    =  (long*)calloc(number_of_particles, sizeof(long));
  for (i = 0; i < number_of_particles; i++){
    particle_no_proton[i] = 1001;
  }

  AT_effective_charge_from_E_MeV_u(  number_of_particles,
      E_MeV_u,
      particle_no_proton,
      Zeff_proton);

  // scale proton LET by ratio of effective Z
  for (i = 0; i < number_of_particles; i++){
    if( particle_no[i] != 1001){ // for particles other than proton scale LET by (Zeff_ion / Zeff_proton)^2
      LET_MeV_cm2_g[i] *=   gsl_pow_2(Zeff_ion[i] / Zeff_proton[i]);
    }
  }

  free(particle_no_proton);
  free(Zeff_proton);

  free(Zeff_ion);
}


void AT_LET_keV_um(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    double        LET_keV_um[])
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
    LET_keV_um[i]  *=  material_density_g_cm3 * 0.1;
  }
}


void AT_CSDA_range_g_cm2(  const long  number_of_particles,
    const double   E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    double         CSDA_range_g_cm2[])
{
  getPSTARvalue(number_of_particles, E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.range_cdsa_g_cm2, CSDA_range_g_cm2);

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


void AT_E_MeV_from_CDSA_range(  const long  number_of_particles,
    const double  CSDA_range_g_cm2[],
    const long   particle_no[],
    const long   material_no,
    double        E_MeV[])
{
  // scaled energies
  double*  sE  =  (double*)calloc(number_of_particles, sizeof(double));

  getPSTARvalue(number_of_particles, CSDA_range_g_cm2, material_no, AT_PSTAR_Data.range_cdsa_g_cm2, AT_PSTAR_Data.kin_E_MeV, sE);

  free( sE );
}


void AT_E_MeV_from_LET(  const long  number_of_particles,
    const double  LET_MeV_cm2_g[],
    const long   particle_no[],
    const long   material_no,
    double        E_MeV[])
{
  // scaled energies
  double*  scaled_E_MeV  =  (double*)calloc(number_of_particles, sizeof(double));

  //TODO add effective charge correction !!

  // Conversion LETion => LETproton
  long*  charge  =  (long*)calloc(number_of_particles, sizeof(long));
  AT_Z_from_particle_no(        number_of_particles,
      particle_no,
      charge);


  // loop over n to find charge for all given particles and energies
  double*  LET_MeV_cm2_g_copy = (double*)calloc(number_of_particles,sizeof(double));
  memcpy(LET_MeV_cm2_g_copy,LET_MeV_cm2_g,number_of_particles);
  long  i;
  for(i = 0; i < number_of_particles; i++){
    LET_MeV_cm2_g_copy[i] /= gsl_pow_2( (double)charge[i] );;
  }

  getPSTARvalue(number_of_particles, LET_MeV_cm2_g_copy, material_no, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, AT_PSTAR_Data.kin_E_MeV, scaled_E_MeV);

  free( charge );
  free( LET_MeV_cm2_g_copy );
  free( scaled_E_MeV );
}
