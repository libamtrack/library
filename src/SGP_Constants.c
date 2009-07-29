/*
 * SGP_Constants.c
 *
 *  Created on: 28.07.2009
 *      Author: greilich
 */
#include "SGP_Constants.h"

#include <string.h>

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
			strcpy(material_name,"");
			break;
		}
}

void getMaterialNo(char* material_name, long* material_no){
	*material_no	= -1;
	if( strcmp(material_name,"Water, Liquid") == 0)
		*material_no = Water_Liquid;
	if( strcmp(material_name,"Aluminum Oxide") == 0)
		*material_no = Aluminum_Oxide;
	if( strcmp(material_name,"Aluminum") == 0)
		*material_no = Aluminum;
	if( strcmp(material_name,"PMMA") == 0)
		*material_no = PMMA;
}
