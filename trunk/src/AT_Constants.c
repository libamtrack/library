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



