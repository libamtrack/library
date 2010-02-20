/**
 * @file
 * @brief Electron Range models
 */

/*
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

  // find look-up index for ER number in ER data table
  long  match;
  const long n_tmp = 1;

  pmatchi(  ER_no,
      &n_tmp,
      AT_ER_Data.ER_no,
      &AT_ER_Data.n,
      &match);

  if( match != -1){
    strcpy(ER_name, AT_ER_Data.ER_name[match]);
  } else {
    strcpy(ER_name,"*** invalid choice ***");
  }
}



void AT_max_electron_range_m( const long*  n,
    const float*  E_MeV_u,
    const long*  particle_no,
    const long*  material_no,
    const long*   er_model,
    float*  max_electron_range_m)
{
  // Get density matching to material_name (only 1 name therefore n_mat = 1)
  const long  n_mat  = 1;
  float material_density_g_cm3;
  float average_A;
  float average_Z;

  AT_getMaterialData( &n_mat, material_no, &material_density_g_cm3,
      NULL,NULL,NULL,NULL,NULL, &average_A, &average_Z );

  float* mass    =  (float*)calloc(*n, sizeof(float));

  AT_mass_from_particle_no(n,particle_no,mass);

  float* beta    =  (float*)calloc(*n, sizeof(float));

  AT_beta_from_mass(n,E_MeV_u,mass,beta);

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
    if( *er_model == ER_Tabata ){
      // general constants (best fit to experimental data)
      const double b1 = 2.335;
      const double b2 = 1.209;
      const double b3 = 1.78e-4;
      const double b4 = 0.9891;
      const double b5 = 3.01e-4;
      const double b6 = 1.468;
      const double b7 = 1.18e-2;
      const double b8 = 1.232;
      const double b9 = 0.109;
      // average A and Z for given material
      const double A = average_A;
      const double Z = average_Z;
      // constants...
      const double a1_g_cm2 = 0.1*b1*A / pow(Z,b2); // g_cm2
      const double a2 = b3*Z;
      const double a3 = b4 - b5*Z;
      const double a4 = b6 - b7*Z;
      const double a5 = b8 / pow(Z,b9);
      double tau = 2.0 * gsl_pow_2(beta[i]) / (1. - gsl_pow_2(beta[i]));
      max_electron_range_m[i] = (a1_g_cm2)*(((gsl_sf_log(1 + a2 * tau))/a2) - ((a3*tau)/(1 + a4*pow(tau,a5))) );
    }

    // Scale maximum el. range with material density relative to water (1/rho)
    max_electron_range_m[i]    /= material_density_g_cm3;

    // covert cm to m
    max_electron_range_m[i]    /= 1e2;  // cm to m

  }
  free(mass);
  free(beta);
}
