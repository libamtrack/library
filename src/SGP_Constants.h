#ifndef SGP_CONSTANTS_H_
#define SGP_CONSTANTS_H_

#include <string.h>

double	c_m_s						=	299792458.0; 						// speed of light [m/s]
double	c_cm_s						=	29979245800; 						// speed of light [cm/s]
double	proton_mass_MeV_c2			=	938.272029;							// proton mass
double	proton_mass_kg				=	1.67262171e-27;
double	electron_mass_MeV_c2		=	0.510998918;						// electron mass
double	electron_mass_kg			=	9.10938215e-31;
double	electron_mass_g				=	9.10938215e-28;
double	e_C							=	1.60217653e-19;						// elementary charge
double	e_esu						=	1.60217653e-19 / 3.33564e-10;						// elementary charge
double	e0_F_m						=	8.8541878176e-12;					// electrical permittivity of the vacuum

double	g_cm3_to_kg_m3				=	1000.0;								// replace later by system of units
double	MeV_to_J					=	1.60217646e-13;
double	MeV_to_keV					=	1000.0;
double	MeV_g_to_J_kg				=	1.60217646e-10;
double	GeV_g_to_J_kg				=	1.60217646e-7;
double	keV_um_to_MeV_m				=	1000.0;
double	m_to_cm						=	100.0;
double	m_to_um						=	1e6;

double	pi							=	3.14159265;

enum GammaResponseModels{
    GR_Test					= 0,
    GR_GeneralTarget		= 1,
    GR_Radioluminescence	= 2,
    GR_ExpSaturation		= 3
};

enum RDDModels{
	RDD_Test				= 0,
	RDD_KatzPoint			= 1,			/* parameters: 0 - beta,	1 - Z,			 2 - material-no */
	RDD_Geiss				= 2,			/* parameters: 0 - E_MeV_u, 1 - particle_no, 2 - material_no, 3 - a0 [m] */
	RDD_Site				= 3				/* parameters: 0 - E_MeV_u, 1 - particle_no, 2 - material_no, 3 - a0 [m] */ //as defined in Edmund et al. , 2007
};

enum material_no{
	Water_Liquid			= 0,
	Aluminum_Oxide			= 1,
 	Aluminum				= 2,
 	PMMA					= 3
};

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

enum ERModels{
	ER_Test					= 0,
	ER_PowerLawZBeta		= 1,
	ER_PowerLawLET			= 2
};

#endif // SGP_CONSTANTS_H_
