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

///////////////////////////////////////////// SINGLE MATERIAL IMPLEMENTATION ///////////////////////////////////////////////////////


long AT_index_from_material_no( const long material_no ){
  long  index                =  -1;
  long  number_of_materials  =  1;
  find_elements_int(  &material_no,
      number_of_materials,
      AT_Material_Data.material_no,
      AT_Material_Data.n,
      &index);   // TODO replace call to pmatchi by call to simpler function which will find the index
  return index;
}


void getMaterialName(
    const long material_no,
    char* material_name){

  long  index = AT_index_from_material_no( material_no );

  if( index != -1){
    strcpy(material_name, AT_Material_Data.material_name[index]);
  } else {
    strcpy(material_name,"*** invalid choice ***");
  }
}


long getMaterialNo(
    const char* material_name){

  // find look-up index for material name in material data table
  long  match;
  const long n_tmp = 1;

  find_elements_char(  &material_name,
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



double AT_density_g_cm3_from_material_no( const long   material_no )
{
  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    return 0.0;
  }
  return AT_Material_Data.density_g_cm3[index];
}


double AT_electron_density_m3_from_material_no( const long   material_no )
{
  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    return 0.0;
  }
  return AT_Material_Data.electron_density_m3[index];
}


double AT_I_eV_from_material_no( const long   material_no )
{
  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    return 0.0;
  }
  return AT_Material_Data.I_eV[index];
}


double AT_alpha_g_cm2_MeV_from_material_no( const long   material_no )
{
  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    return 0.0;
  }
  return AT_Material_Data.alpha_g_cm2_MeV[index];
}

double AT_p_MeV_from_material_no( const long   material_no )
{
  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    return 0.0;
  }
  return AT_Material_Data.p_MeV[index];
}


double AT_m_g_cm2_from_material_no( const long   material_no )
{
  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    return 0.0;
  }
  return AT_Material_Data.m_g_cm2[index];
}

double AT_average_A_from_material_no( const long   material_no )
{
  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    return 0.0;
  }
  return AT_Material_Data.average_A[index];
}


double AT_average_Z_from_material_no( const long   material_no )
{
  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    return 0.0;
  }
  return AT_Material_Data.average_Z[index];
}


void AT_get_material_data(     const long  material_no,
    double*  density_g_cm3,
    double*  electron_density_m3,
    double*  I_eV,
    double*  alpha_g_cm2_MeV,
    double*  p_MeV,
    double*  m_g_cm2,
    double*  average_A,
    double*  average_Z){

  long  index = AT_index_from_material_no( material_no );
  if( index == -1){
    printf("Material no %ld not found\n", material_no);
    *density_g_cm3        =  0.0;
    *electron_density_m3  =  0.0;
    *I_eV                 =  0.0;
    *alpha_g_cm2_MeV      =  0.0;
    *p_MeV                =  0.0;
    *m_g_cm2              =  0.0;
    *average_A            =  0.0;
    *average_Z            =  0.0;
    return ;
  }
  if( density_g_cm3 != NULL ){
      *density_g_cm3            =  AT_Material_Data.density_g_cm3[index];
  }
  if( electron_density_m3 != NULL ){
      *electron_density_m3      =  AT_Material_Data.electron_density_m3[index];
  }
  if( I_eV != NULL ){
      *I_eV                     =  AT_Material_Data.I_eV[index];
  }
  if( alpha_g_cm2_MeV != NULL ){
      *alpha_g_cm2_MeV          =  AT_Material_Data.alpha_g_cm2_MeV[index];
  }
  if( p_MeV != NULL ){
      *p_MeV                    =  AT_Material_Data.p_MeV[index];
  }
  if( m_g_cm2 != NULL ){
      *m_g_cm2                  =  AT_Material_Data.m_g_cm2[index];
  }
  if( average_A != NULL ){
      *average_A                =  AT_Material_Data.average_A[index];
  }
  if( average_Z != NULL ){
      *average_Z                =  AT_Material_Data.average_Z[index];
  }
}


///////////////////////////////////////////// MULTIPLE MATERIAL IMPLEMENTATION ///////////////////////////////////////////////////////


void AT_get_materials_data( const long  number_of_materials,
    const long  material_no[],
    double  density_g_cm3[],
    double  electron_density_m3[],
    double  I_eV[],
    double  alpha_g_cm2_MeV[],
    double  p_MeV[],
    double  m_g_cm2[],
    double  average_A[],
    double  average_Z[])
{
  long*  match  =  (long*)calloc(number_of_materials, sizeof(long));
  find_elements_int(  material_no,
      number_of_materials,
      AT_Material_Data.material_no,
      AT_Material_Data.n,
      match);

  long i;
  if( density_g_cm3 != NULL ){
    for(i = 0; i < number_of_materials; i++){
      density_g_cm3[i]      = AT_Material_Data.density_g_cm3[match[i]];
    }
  }
  if( electron_density_m3 != NULL ){
    for(i = 0; i < number_of_materials; i++){
      electron_density_m3[i]      = AT_Material_Data.electron_density_m3[match[i]];
    }
  }
  if( I_eV != NULL ){
    for(i = 0; i < number_of_materials; i++){
      I_eV[i]      = AT_Material_Data.I_eV[match[i]];
    }
  }
  if( alpha_g_cm2_MeV != NULL ){
    for(i = 0; i < number_of_materials; i++){
      alpha_g_cm2_MeV[i]      = AT_Material_Data.alpha_g_cm2_MeV[match[i]];
    }
  }
  if( p_MeV != NULL ){
    for(i = 0; i < number_of_materials; i++){
      p_MeV[i]      = AT_Material_Data.p_MeV[match[i]];
    }
  }
  if( m_g_cm2 != NULL ){
    for(i = 0; i < number_of_materials; i++){
      m_g_cm2[i]      = AT_Material_Data.m_g_cm2[match[i]];
    }
  }
  if( average_A != NULL ){
    for(i = 0; i < number_of_materials; i++){
      average_A[i]      = AT_Material_Data.average_A[match[i]];
    }
  }
  if( average_Z != NULL ){
    for(i = 0; i < number_of_materials; i++){
      average_Z[i]      = AT_Material_Data.average_Z[match[i]];
    }
  }
  free(match);
}


