/**
 * @file
 * @brief ...
 */

/*
*    AT_DataMaterial.c
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

#include "AT_DataMaterial.h"

void getMaterialName(
    const long material_no,
    char* material_name){

  // find look-up index for material number in material data table
  long  match;
  const long n_tmp = 1;

  pmatchi(  &material_no,
      &n_tmp,
      AT_Material_Data.material_no,
      &AT_Material_Data.n,
      &match);

  if( match != -1){
    strcpy(material_name, AT_Material_Data.material_name[match]);
  } else {
    strcpy(material_name,"*** invalid choice ***");
  }
}

long getMaterialNo(
    const char* material_name){

  // find look-up index for material name in material data table
  long  match;
  const long n_tmp = 1;

  pmatchc(  &material_name,
      &n_tmp,
      AT_Material_Data.material_name,
      &AT_Material_Data.n,
      &match);

  if( match != -1){
    return AT_Material_Data.material_no[match];
  } else {
    return -1;
  }
}

//TODO function to get properties of only one material can be useful
void AT_getMaterialData( const long  n,
    const long*  material_no,
    double*  density_g_cm3,
    double*  electron_density_m3,
    double*  I_eV,
    double*  alpha_g_cm2_MeV,
    double*  p_MeV,
    double*  m_g_cm2,
    double*  average_A,
    double*  average_Z)
{
  long*  match  =  (long*)calloc(n, sizeof(long));
  pmatchi(  material_no,
      &n,
      AT_Material_Data.material_no,
      &AT_Material_Data.n,
      match);

  long i;
  if( density_g_cm3 != NULL ){
    for(i = 0; i < n; i++){
      density_g_cm3[i]      = AT_Material_Data.density_g_cm3[match[i]];
    }
  }
  if( electron_density_m3 != NULL ){
    for(i = 0; i < n; i++){
      electron_density_m3[i]      = AT_Material_Data.electron_density_m3[match[i]];
    }
  }
  if( I_eV != NULL ){
    for(i = 0; i < n; i++){
      I_eV[i]      = AT_Material_Data.I_eV[match[i]];
    }
  }
  if( alpha_g_cm2_MeV != NULL ){
    for(i = 0; i < n; i++){
      alpha_g_cm2_MeV[i]      = AT_Material_Data.alpha_g_cm2_MeV[match[i]];
    }
  }
  if( p_MeV != NULL ){
    for(i = 0; i < n; i++){
      p_MeV[i]      = AT_Material_Data.p_MeV[match[i]];
    }
  }
  if( m_g_cm2 != NULL ){
    for(i = 0; i < n; i++){
      m_g_cm2[i]      = AT_Material_Data.m_g_cm2[match[i]];
    }
  }
  if( average_A != NULL ){
    for(i = 0; i < n; i++){
      average_A[i]      = AT_Material_Data.average_A[match[i]];
    }
  }
  if( average_Z != NULL ){
    for(i = 0; i < n; i++){
      average_Z[i]      = AT_Material_Data.average_Z[match[i]];
    }
  }
  free(match);
}

void AT_density_g_cm3_from_material_no( const long  n,
    const long*  material_no,
    double*       density_g_cm3)
{
  AT_getMaterialData(n,material_no,density_g_cm3,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}

