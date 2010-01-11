/**
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
    const long* material_no,
    char* material_name){
  switch( (int)(*material_no) ){
  case Water_Liquid:
    strcpy(material_name,"Water, Liquid");
    break;
  case Aluminum_Oxide:
    strcpy(material_name,"Aluminum Oxide");
    break;
  case Aluminum:
    strcpy(material_name,"Aluminum");
    break;
  case PMMA:
    strcpy(material_name,"PMMA");
    break;
  default:
    strcpy(material_name,"*** invalid choice ***");
    break;
  }
}

void getMaterialNo(
    const char* material_name,
    long* material_no){
  *material_no  = -1;
  if( strcmp(material_name,"Water, Liquid") == 0)
    *material_no = Water_Liquid;
  if( strcmp(material_name,"Aluminum Oxide") == 0)
    *material_no = Aluminum_Oxide;
  if( strcmp(material_name,"Aluminum") == 0)
    *material_no = Aluminum;
  if( strcmp(material_name,"PMMA") == 0)
    *material_no = PMMA;
}

//TODO function to get properties of only one material can be useful
void AT_getMaterialData(    const long*  n,
    const long*   material_no,
    float*  density_g_cm3,
    float*  electron_density_m3,
    float*  I_eV,
    float*  alpha_g_cm2_MeV,
    float*  p_MeV,
    float*  m_g_cm2)
{
  long*  match  =  (long*)calloc(*n, sizeof(long));
  pmatchi(  material_no,
      n,
      AT_Material_Data.material_no,
      &AT_Material_Data.n,
      match);

  long i;
  if( density_g_cm3 != NULL ){
    for(i = 0; i < *n; i++){
      density_g_cm3[i]      = AT_Material_Data.density_g_cm3[match[i]];
    }
  }
  if( electron_density_m3 != NULL ){
    for(i = 0; i < *n; i++){
      electron_density_m3[i]      = AT_Material_Data.electron_density_m3[match[i]];
    }
  }
  if( I_eV != NULL ){
    for(i = 0; i < *n; i++){
      I_eV[i]      = AT_Material_Data.I_eV[match[i]];
    }
  }
  if( alpha_g_cm2_MeV != NULL ){
    for(i = 0; i < *n; i++){
      alpha_g_cm2_MeV[i]      = AT_Material_Data.alpha_g_cm2_MeV[match[i]];
    }
  }
  if( p_MeV != NULL ){
    for(i = 0; i < *n; i++){
      p_MeV[i]      = AT_Material_Data.p_MeV[match[i]];
    }
  }
  if( m_g_cm2 != NULL ){
    for(i = 0; i < *n; i++){
      m_g_cm2[i]      = AT_Material_Data.m_g_cm2[match[i]];
    }
  }
  free(match);
}

void AT_density_g_cm3_from_material_no( const long*  n,
    const long*  material_no,
    float*       density_g_cm3)
{
  AT_getMaterialData(n,material_no,density_g_cm3,NULL,NULL,NULL,NULL,NULL);
}

//TODO should be removed and also all references
void AT_density_g_cm3_from_material_name( const long*  n,
    const char*  material_name,
    float*       density_g_cm3)
{
  long  i;
  long  match;
  long  n_mat  = 1;
  pmatchc(  &material_name,
      &n_mat,
      AT_Material_Data.material_name,
      &AT_Material_Data.n,
      &match);
  for(i = 0; i < *n; i++){
    density_g_cm3[i]    = AT_Material_Data.density_g_cm3[match];
  }
}

void AT_density_g_cm3_from_material_nameS(  char**  material_name,
    float*  density_g_cm3){
  long  n;
  n    = 1;
  AT_density_g_cm3_from_material_name(  &n,
      *material_name,
      density_g_cm3);
}

//TODO should be removed and also all references
void AT_electron_density_m3_from_material_name(  const long*  n,
    const char*  material_name,
    float*  electron_density_m3)
{
  long  i;
  long  match;
  long  n_mat  = 1;
  pmatchc(  &material_name,
      &n_mat,
      AT_Material_Data.material_name,
      &AT_Material_Data.n,
      &match);
  for(i = 0; i < *n; i++){
    electron_density_m3[i]  = AT_Material_Data.electron_density_m3[match];
  }
}

void AT_electron_density_m3_from_material_nameS(  char**  material_name,
    float*  electron_density_m3){
  long  n;
  n    = 1;
  AT_electron_density_m3_from_material_name(  &n,
      *material_name,
      electron_density_m3);
}
