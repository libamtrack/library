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
    const long*  material_no,
    const long*   er_model,
    float*  max_electron_range_m)
{

  /********************************************************
   ********* CALCULATION BEFORE PARTICLE LOOP *************
   *******************************************************/

  // Get density matching to material_name (only 1 name therefore n_mat = 1)
  const long  n_mat  = 1;
  float material_density_g_cm3;
  float average_A;
  float average_Z;
  AT_getMaterialData( &n_mat, material_no, &material_density_g_cm3,
      NULL,NULL,NULL,NULL,NULL, &average_A, &average_Z );

  // Get beta from energy
  float* beta    =  (float*)calloc(*n, sizeof(float));
  AT_beta_from_E(n,E_MeV_u,beta);

  // Get energy of delta-electron from energy of ion
  float* wmax_MeV =  (float*)calloc(*n, sizeof(float));
  AT_max_E_transfer_MeV(n,E_MeV_u,wmax_MeV);

  // a1,..a5 needed in Tabata ER model
  double a1_g_cm2 = 0;
  double a2 = 0;
  double a3 = 0;
  double a4 = 0;
  double a5 = 0;
  if( *er_model == ER_Tabata ){
    // general constants (best fit to experimental data)
    const double b1_g_cm2 = 0.2335;
    const double b2 = 1.209;
    const double b3 = 1.78e-4;
    const double b4 = 0.9891;
    const double b5 = 3.01e-4;
    const double b6 = 1.468;
    const double b7 = 1.18e-2;
    const double b8 = 1.232;
    const double b9 = 0.109;
    // constants...
    a1_g_cm2 = b1_g_cm2*average_A / pow(average_Z,b2); // g_cm2
    a2 = b3*average_A;
    a3 = b4 - b5*average_Z;
    a4 = b6 - b7*average_Z;
    a5 = b8 / pow(average_Z,b9);
  }

  float max_electron_range_g_cm2 = 0.0f;

  /********************************************************
   *********************  PARTICLE LOOP *******************
   *******************************************************/
  long  i;
  for (i = 0; i < *n; i++){
    float tmpE_MeV_u  = E_MeV_u[i];
    if (tmpE_MeV_u < 0) {tmpE_MeV_u *= -1.0f;}  // E can be set neg. if non-PSTAR are given --> use pos. value //TODO what does it mean ?

    float wmax_keV = wmax_MeV[i] * 1000.0f;

    if( *er_model == ER_ButtsKatz ){
      max_electron_range_g_cm2 = 1e-5 * wmax_keV;
    }
    if( *er_model == ER_Waligorski ){
      double alpha = 1.667;
      if( wmax_keV < 1. ) alpha = 1.079;
      max_electron_range_g_cm2 =  6* 1e-6 * (float)pow( wmax_keV, alpha );
    }
    if( *er_model == ER_Edmund ){
      double alpha = 1.67;
      if( wmax_keV < 1. ) alpha = 1.079;
      max_electron_range_g_cm2 =  6.13*1e-6 * (float)pow( wmax_keV, alpha );
    }
    if( *er_model == ER_Geiss ){
      max_electron_range_g_cm2 = 4e-5 * (float)pow(tmpE_MeV_u, 1.5);
    }
    if( *er_model == ER_Scholz ){
      max_electron_range_g_cm2 = 5e-5 * (float)pow(tmpE_MeV_u, 1.7);
    }
    if( *er_model == ER_Tabata ){
      double tau = 2.0f * gsl_pow_2(beta[i]) / (1. - gsl_pow_2(beta[i]));
      max_electron_range_g_cm2 = (a1_g_cm2)*(((gsl_sf_log(1 + a2 * tau))/a2) - ((a3*tau)/(1 + a4*pow(tau,a5))) );
    }

    // Scale maximum el. range with material density relative to water (1/rho) and convert cm to m
    max_electron_range_m[i]    = 1e-2 * max_electron_range_g_cm2 / material_density_g_cm3;

  }
  free(beta);
  free(wmax_MeV);
}
