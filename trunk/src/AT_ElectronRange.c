/**
*    AT_ElectronRange.c
*    ==============
*
*    Created on: 08.01.2010
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

#include "AT_ElectronRange.h"

void getERName(
    const long* ER_no,
    char* ER_name){
  switch( (int)(*ER_no) ){
  case ER_Test:
    strcpy(ER_name,"simple test ER model");
    break;
  case ER_ButtsKatz:
    strcpy(ER_name,"Butts & Katz' [Katz et al., 1972] ER model");
    break;
  case ER_Waligorski:
    strcpy(ER_name,"Waligorski's ER model");
    break;
  case ER_Geiss:
    strcpy(ER_name,"Geiss' [Geiss, 1997] ER model");
    break;
  case ER_Scholz:
    strcpy(ER_name,"ER_Scholz' [Scholz, 2001] ER model");
    break;
  default:
    strcpy(ER_name,"*** invalid choice ***");
    break;
  }
}


#ifdef _R
void AT_max_electron_range_mS(  int*  n,
    float*  E_MeV_u,
    int*  particle_no,
    int*  material_no,
    int*   er_model,
    float*  max_electron_range_m)
{
  long n_long = (long)(*n);
  long material_no_long = (long)(*material_no);
  long er_model_long = (long)(*er_model);

  long i;
  long * particle_no_long = (long*)calloc(*n,sizeof(long));
  for(i = 0 ; i < *n ; i++){
    particle_no_long[i] = (long)particle_no[i];
  }

  AT_max_electron_range_m( &n_long,
      E_MeV_u,
      particle_no_long,
      &material_no_long,
      &er_model_long,
      max_electron_range_m);

  free(particle_no_long);

}
#endif


void AT_max_electron_range_m( const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const long*  material_no,
    const long*   er_model,
    float*  max_electron_range_m)
{
  //TODO change to another routine (avoid pmatchi), like get_density_of_material

  // Get density matching to material_name (only 1 name therefore n_mat = 1)
  long  n_mat  = 1;
  long  match;
  pmatchi(  material_no,
        &n_mat,
        AT_Material_Data.material_no,
        &AT_Material_Data.n,
        &match);

  float* mass    =  (float*)calloc(*n, sizeof(float));

  AT_mass_from_particle_no(n,particle_no,mass);

  long  i;
  for (i = 0; i < *n; i++){
    float tmpE  = E_MeV_u[i];
    if (tmpE < 0) {tmpE *= -1.0f;}  // E can be set neg. if non-PSTAR are given --> use pos. value

    float E_div_E0 = E_MeV_u[i] / (mass[i]*proton_mass_MeV_c2);
    float w_keV;
    if( *er_model == ER_ButtsKatz ){
      w_keV = 2 * electron_mass_MeV_c2 * ( E_div_E0*E_div_E0 + 2*E_div_E0) * 1e3;
      max_electron_range_m[i] = 1e-5 * w_keV;
    }
    if( *er_model == ER_Waligorski ){
      double alpha = 1.667;
      w_keV = 2 * electron_mass_MeV_c2 * ( E_div_E0*E_div_E0 + 2*E_div_E0) * 1e3;
      if( w_keV < 1. ) alpha = 1.079;
      max_electron_range_m[i] =  6* 1e-6 * (float)pow( w_keV, alpha );
    }
    if( *er_model == ER_Edmund ){
      double alpha = 1.67;
      w_keV = 2 * electron_mass_MeV_c2 * ( E_div_E0*E_div_E0 + 2*E_div_E0) * 1e3;
      if( w_keV < 1. ) alpha = 1.079;
      max_electron_range_m[i] =  6.13*1e-6 * (float)pow( w_keV, alpha );
    }
    if( *er_model == ER_Geiss ){
      max_electron_range_m[i] = 4e-5 * (float)pow(tmpE, 1.5);
    }
    if( *er_model == ER_Scholz ){
      max_electron_range_m[i] = 5e-5 * (float)pow(tmpE, 1.7);
    }

    // Scale maximum el. range with material density relative to water (1/rho)
    max_electron_range_m[i]    /= AT_Material_Data.density_g_cm3[match];

    // covert cm to m
    max_electron_range_m[i]    /= 1e2;  // cm to m

  }
  free(mass);
}
