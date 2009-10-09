/**
*    AT_Constants.c
*    ==============
*
*    Created on: 28.07.2009
*    Author: greilich
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

#include "AT_Constants.h"

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

void getRDDName(long* RDD_no, char* RDD_name){
  strcpy(RDD_name,"*** invalid choice ***");
  long i;
  for (i = 0; i < RDD_DATA_N; i++){
    if (AT_RDD_Data.RDD_no[i] == *RDD_no){
      strcpy(RDD_name, AT_RDD_Data.RDD_name[i]);
    }
  }
}

void getRDDNo(char* RDD_name, long* RDD_no){
  *RDD_no = 0;
  long i;
  for (i = 0; i < RDD_DATA_N; i++){
    if (strcmp(RDD_name, AT_RDD_Data.RDD_name[i]) == 0){
      *RDD_no = AT_RDD_Data.RDD_no[i];
      break;
    }
  }
}

void getERName(long* ER_no, char* ER_name){
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

void getGammaName(long* Gamma_no, char* Gamma_name){
  switch( (int)(*Gamma_no) ){
  case GR_Test:
    strcpy(Gamma_name,"simple test gamma response");
    break;
  case GR_GeneralTarget:
    strcpy(Gamma_name,"generalized multi-target/multi-hit gamma response");
    break;
  case GR_Radioluminescence:
    strcpy(Gamma_name,"radioluminescence gamma response");
    break;
  case GR_ExpSaturation:
    strcpy(Gamma_name,"exp.-sat. gamma response (obsolete, use gen. target/hit instead)");
    break;
  case GR_LinQuad:
    strcpy(Gamma_name,"linear-quadratic gamma response");
    break;
  case GR_LinQuad_Log:
    strcpy(Gamma_name,"lethal events number response");
    break;
  default:
    strcpy(Gamma_name,"*** invalid choice ***");
    break;
  }
}


