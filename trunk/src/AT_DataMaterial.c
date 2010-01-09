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

void getMaterialName(long* material_no, char* material_name){
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

void getMaterialNo(char* material_name, long* material_no){
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


///////////////////////////////////////////////////////////////////////
// Routines to access MATERIAL data
///////////////////////////////////////////////////////////////////////
void AT_getMaterialData(    long*  n,
    long*  material_no,
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
  for(i = 0; i < *n; i++){
    density_g_cm3[i]      = AT_Material_Data.density_g_cm3[match[i]];
    electron_density_m3[i]= AT_Material_Data.electron_density_m3[match[i]];
    I_eV[i]               = AT_Material_Data.I_eV[match[i]];
    alpha_g_cm2_MeV[i]    = AT_Material_Data.alpha_g_cm2_MeV[match[i]];
    p_MeV[i]              = AT_Material_Data.p_MeV[match[i]];
    m_g_cm2[i]            = AT_Material_Data.m_g_cm2[match[i]];
  }

  free(match);
}


#define matchIt      long  n_mat  = 1;                  \
                pmatchc(  &material_name,                \
                                &n_mat,                    \
                                AT_Material_Data.material_name,      \
                                &AT_Material_Data.n,            \
                                &match);

void AT_density_g_cm3( long*  n,
    char*  material_name,
    float*  density_g_cm3)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    density_g_cm3[i]    = AT_Material_Data.density_g_cm3[match];
  }
}

void AT_density_g_cm3S(    char**  material_name,
    float*  density_g_cm3){
  long  n;
  n    = 1;
  AT_density_g_cm3(  &n,
      *material_name,
      density_g_cm3);
}

void AT_electron_density_m3(  long*  n,
    char*  material_name,
    float*  electron_density_m3)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    electron_density_m3[i]  = AT_Material_Data.electron_density_m3[match];
  }
}

void AT_electron_density_m3S(  char**  material_name,
    float*  electron_density_m3){
  long  n;
  n    = 1;
  AT_electron_density_m3(  &n,
      *material_name,
      electron_density_m3);
}

void AT_alpha_g_cm2_MeV(    long*  n,
    char*  material_name,
    float*  alpha_g_cm2_MeV)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    alpha_g_cm2_MeV[i]    = AT_Material_Data.alpha_g_cm2_MeV[match];
  }
}

void AT_p_MeV(          long*  n,
    char*  material_name,
    float*  p_MeV)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    p_MeV[i]        = AT_Material_Data.p_MeV[match];
  }
}

void AT_m_g_cm2(        long*  n,
    char*  material_name,
    float*  m_g_cm2)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    m_g_cm2[i]        = AT_Material_Data.m_g_cm2[match];
  }
}


