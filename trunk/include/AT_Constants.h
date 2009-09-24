#ifndef AT_CONSTANTS_H_
#define AT_CONSTANTS_H_

/**
 *    AT_Constants.h
 *    ==============
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

#include <string.h>
#include <stdio.h>

#define _LINUX // [_LINUX or _WINDOWS] : in Linux we have isnan function while in Windows we have _isnan
#define _R // [_S or _R] in S we can pass long type to the function via as.single, but in R we pass int type
//#define _DEBUG // debugging printouts

extern int indent_counter;
extern char isp[];
extern FILE * debf;

static  double  c_m_s                 =  299792458.0;                          // speed of light [m/s]
static  double  c_cm_s                =  29979245800.0;                        // speed of light [cm/s]
static  double  proton_mass_MeV_c2    =  938.272029;                           // proton mass
static  double  proton_mass_kg        =  1.67262171e-27;
static  double  electron_mass_MeV_c2  =  0.510998918;                          // electron mass
static  double  electron_mass_kg      =  9.10938215e-31;
static  double  electron_mass_g       =  9.10938215e-28;
static  double  e_C                   =  1.60217653e-19;                       // elementary charge
static  double  e_esu                 =  1.60217653e-19 / 3.33564e-10;         // elementary charge
static  double  e0_F_m                =  8.8541878176e-12;                     // electrical permittivity of the vacuum

static  double  g_cm3_to_kg_m3        =  1000.0;                               // replace later by system of units
static  double  MeV_to_J              =  1.60217646e-13;
static  double  MeV_to_keV            =  1000.0;
static  double  MeV_g_to_J_kg         =  1.60217646e-10;
static  double  GeV_g_to_J_kg         =  1.60217646e-7;
static  double  keV_um_to_MeV_m       =  1000.0;
static  double  m_to_cm               =  100.0;
static  double  m_to_um               =  1e6;

static  double  pi                    =  3.14159265;

enum Method{
	ME_Grid           = 1,
			ME_SPIFF      = 2,
			ME_Katz       = 3
};

enum GammaResponseModels{
  GR_Test                  = 1,      /* no parameters */
      GR_GeneralTarget     = 2,      /* */
      GR_Radioluminescence = 3,      /* 0 - Smax, 1 - D0, 2 - dyn */
      GR_ExpSaturation     = 4,      /* 0 - Smax, 1 - D0 */
      GR_LinQuad           = 5        /* 0 - alpha, 1 - beta */
};

enum RDDModels{
  RDD_Test                 = 1,      /* no parameters */
      RDD_KatzPoint        = 2,      /* parameters: 0 - r_min [m] (lower integration limit), 1 - d_min_Gy (lower dose cut-off) */
      RDD_Geiss            = 3,      /* parameters: 0 - a0 [m] (core diameter) */
      RDD_Site             = 4,      /* parameters: 0 - a0 [m] (core diameter), 1 - d_min_Gy (lower dose cut-off)  */ // after Edmund et al., 2007, but modified with dose-cut off
      RDD_ExtTarget        = 5       /* parameters: 0 - r_min [m] (core diameter), 1 - a0 [m] (target diameter), 2 - D_min [Gy] (cut-off dose) */ //as defined in Edmund et al. , 2007
};

enum material_no{
  Water_Liquid             = 1,
      Aluminum_Oxide       = 2,
      Aluminum             = 3,
      PMMA                 = 4
};

enum ERModels{
  ER_Test                  = 0,
      ER_ButtsKatz         = 1,
      ER_Waligorski        = 2,
      ER_Geiss             = 3,
      ER_Scholz            = 4
};

void   getMaterialName(  long* material_no, char* material_name);
void   getMaterialNo(    char* material_name, long* material_no);

void   getRDDName(  long* RDD_no, char* RDD_name);
void   getERName(  long* ER_no, char* ER_name);
void   getGammaName(  long* Gamma_no, char* Gamma_name);
void   getMethodName(  long* Method_no, char* Method_name);


#endif // AT_CONSTANTS_H_
