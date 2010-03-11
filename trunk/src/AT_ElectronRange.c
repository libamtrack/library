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
    const int ER_no,
    char* ER_name){

  // find look-up index for ER number in ER data table
  int  match;

  pmatchi_int(  &ER_no,
      1,
      AT_ER_Data.ER_no,
      AT_ER_Data.n,
      &match);

  if( match != -1){
    strcpy(ER_name, AT_ER_Data.ER_name[match]);
  } else {
    strcpy(ER_name,"*** invalid choice ***");
  }
}


inline double AT_ER_ButtsKatz_range_g_cm2(double wmax_keV){
  return 1e-5 * wmax_keV;
}


inline double AT_ER_Waligorski_range_g_cm2(double wmax_keV){
  double alpha = 1.667;
  if( wmax_keV < 1. ) alpha = 1.079;
  return  6e-6 * pow( wmax_keV, alpha );
}


inline double AT_ER_Edmund_range_g_cm2(double wmax_keV){
  double alpha = 1.667;
  if( wmax_keV < 1. ) alpha = 1.079;
  return 6.13*1e-6  * pow( wmax_keV, alpha );
}


inline double AT_ER_Geiss_range_g_cm2(double E_MeV_u){
  return 4e-5 * pow(E_MeV_u, 1.5);
}


inline double AT_ER_Scholz_range_g_cm2(double E_MeV_u){
  return 5e-5 * pow(E_MeV_u, 1.7);
}


inline double AT_ER_Tabata_range_g_cm2(double beta, double a1_g_cm2, double a2, double a3, double a4, double a5){
  double tau = 2.0 * gsl_pow_2(beta) / (1.0 - gsl_pow_2(beta));
  return (a1_g_cm2)*(((gsl_sf_log(1.0 + a2 * tau))/a2) - ((a3*tau)/(1.0 + a4*pow(tau,a5))) );
}


inline void AT_ER_Tabata_constants(const double average_A, const double average_Z, double * a1_g_cm2, double * a2, double * a3, double * a4, double * a5){
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
  *a1_g_cm2 = b1_g_cm2*average_A / pow(average_Z,b2);
  *a2 = b3*average_A;
  *a3 = b4 - b5*average_Z;
  *a4 = b6 - b7*average_Z;
  *a5 = b8 / pow(average_Z,b9);
}


//TODO replace float by double
void AT_max_electron_ranges_m( const long  number_of_particles,
    const float  E_MeV_u[],
    const int    material_no,
    const int    er_model,
    float  max_electron_range_m[])
{

  /********************************************************
   ********* CALCULATION BEFORE PARTICLE LOOP *************
   *******************************************************/

  // Get density matching to material_name (only 1 name therefore n_mat = 1)
  const long  number_of_materials_tmp  =  1;
  double material_density_g_cm3        =  0.0;
  double average_A                     =  0.0;
  double average_Z                     =  0.0;
  long  material_number_long = (long)material_no;
  AT_get_materials_data( number_of_materials_tmp, &material_number_long, &material_density_g_cm3,
      NULL,NULL,NULL,NULL,NULL, &average_A, &average_Z );

  // Get beta from energy
  float* beta     =  (float*)calloc(number_of_particles, sizeof(float));

  // Get energy of delta-electron from energy of ion
  float* wmax_MeV  =  (float*)calloc(number_of_particles, sizeof(float));
  AT_max_E_transfer_MeV(&number_of_particles,E_MeV_u,wmax_MeV);

  // a1,..a5 needed in Tabata ER model
  double a1_g_cm2 = 0.0;
  double a2 = 0.0;
  double a3 = 0.0;
  double a4 = 0.0;
  double a5 = 0.0;
  if( er_model == ER_Tabata ){
    AT_beta_from_E(&number_of_particles,E_MeV_u,beta);

    // general constants (best fit to experimental data)
    AT_ER_Tabata_constants(average_A , average_Z, &a1_g_cm2, &a2, &a3, &a4, &a5);
  }

  float max_electron_range_g_cm2 = 0.0f;

  /********************************************************
   *********************  PARTICLE LOOP *******************
   *******************************************************/
  long  i;
  for (i = 0; i < number_of_particles; i++){
    float wmax_keV = wmax_MeV[i] * 1000.0f;

    switch( er_model ){
      case ER_ButtsKatz :
        max_electron_range_g_cm2  =  AT_ER_ButtsKatz_range_g_cm2(wmax_keV);
        break;
      case ER_Waligorski :
        max_electron_range_g_cm2  =  AT_ER_Waligorski_range_g_cm2(wmax_keV);
        break;
      case ER_Edmund :
        max_electron_range_g_cm2  =  AT_ER_Edmund_range_g_cm2(wmax_keV);
        break;
      case ER_Geiss :
        max_electron_range_g_cm2  =  AT_ER_Geiss_range_g_cm2(E_MeV_u[i]);
        break;
      case ER_Scholz :
        max_electron_range_g_cm2  =  AT_ER_Scholz_range_g_cm2(E_MeV_u[i]);
        break;
      case ER_Tabata :
        max_electron_range_g_cm2  =  AT_ER_Tabata_range_g_cm2(beta[i], a1_g_cm2, a2, a3, a4, a5);
        break;
      default:
        max_electron_range_g_cm2  = 0.0;
        break;
    }

    // Scale maximum el. ranges with material density relative to water (1/rho) and convert cm to m
    max_electron_range_m[i]      =  1e-2 * max_electron_range_g_cm2 / (float)material_density_g_cm3;

  }
  free(beta);
  free(wmax_MeV);
}

double AT_max_electron_range_m(  const double E_MeV_u,
    const int    material_no,
    const int    er_model){

  // Get density matching to material_name (only 1 name therefore n_mat = 1)
  const long  number_of_materials_tmp  = 1;
  double material_density_g_cm3        = 0.0;
  double average_A                     = 0.0;
  double average_Z                     = 0.0;
  long  material_number_long = (long)material_no;
  AT_get_materials_data( number_of_materials_tmp, &material_number_long, &material_density_g_cm3,
      NULL,NULL,NULL,NULL,NULL, &average_A, &average_Z );


  // Get energy of delta-electron from energy of ion
  float  wmax_MeV = 0.0;
  const long  number_of_particles_tmp  = 1;
  const float E_MeV_u_float = (float)E_MeV_u;
  AT_max_E_transfer_MeV(&number_of_particles_tmp,&E_MeV_u_float,&wmax_MeV);
  double wmax_keV = (double)wmax_MeV * 1000.0;

  // a1,..a5 needed in Tabata ER model
  double a1_g_cm2 = 0.0;
  double a2 = 0.0;
  double a3 = 0.0;
  double a4 = 0.0;
  double a5 = 0.0;
  // Get beta from energy - needed only in Tabata model
  float  beta  =  0.0;
  if( er_model == ER_Tabata ){
    AT_beta_from_E(&number_of_particles_tmp,&E_MeV_u_float,&beta);

    // general constants (best fit to experimental data)
    AT_ER_Tabata_constants(average_A , average_Z, &a1_g_cm2, &a2, &a3, &a4, &a5);
  }

  double max_electron_range_g_cm2 = 0.0;

  switch( er_model ){
    case ER_ButtsKatz :
      max_electron_range_g_cm2  =  AT_ER_ButtsKatz_range_g_cm2(wmax_keV);
      break;
    case ER_Waligorski :
      max_electron_range_g_cm2  =  AT_ER_Waligorski_range_g_cm2(wmax_keV);
      break;
    case ER_Edmund :
      max_electron_range_g_cm2  =  AT_ER_Edmund_range_g_cm2(wmax_keV);
      break;
    case ER_Geiss :
      max_electron_range_g_cm2  =  AT_ER_Geiss_range_g_cm2(E_MeV_u);
      break;
    case ER_Scholz :
      max_electron_range_g_cm2  =  AT_ER_Scholz_range_g_cm2(E_MeV_u);
      break;
    case ER_Tabata :
      max_electron_range_g_cm2  =  AT_ER_Tabata_range_g_cm2(beta, a1_g_cm2, a2, a3, a4, a5);
      break;
    default:
      max_electron_range_g_cm2  =  0.0;
      break;
  }

  // Scale maximum el. range with material density relative to water (1/rho) and convert cm to m
  double max_electron_range_m   = 1e-2 * max_electron_range_g_cm2 / material_density_g_cm3;

  return max_electron_range_m;
}
