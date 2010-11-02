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

///////////////////////////////////////////////////////////////////////
// GammaResponse DATA

enum GammaResponseModels{
  GR_Test                  = 0,      /* no parameters */
      GR_GeneralTarget     = 1,      /* */
      GR_Radioluminescence = 2,      /* 0 - Smax, 1 - D0, 2 - dyn */
      GR_ExpSaturation     = 3,      /* 0 - Smax, 1 - D0 */
      GR_LinQuad           = 4        /* 0 - alpha, 1 - beta */
};

#define GR_DATA_N    5

typedef struct {
  long    n;
  long    GR_no[GR_DATA_N];
  long	  n_parameters[GR_DATA_N];
  char**  parameter_name[GR_DATA_N][4];
  float   parameter_default[GR_DATA_N][4];
  char*   GR_name[GR_DATA_N];
} gr_data;

static const gr_data AT_GR_Data = {
		GR_DATA_N,
    {  GR_Test,          GR_GeneralTarget,          GR_Radioluminescence,        GR_ExpSaturation, GR_LinQuad},
    {  0, 4, 3, 2, 2},
    {  {"","","",""},{"S_max", "D0_Gy", "c", "m"},{"S_max","D0_Gy","dyn",""},{"S_max","D0_Gy","",""},{"alpha","beta","",""}},
    {  {0,0,0,0}, {1, 10, 1, 1}, {1,10,5,0}, {1,10,0,0}, {1, 1, 0, 0}},
    {  "simple test gamma response",  "generalized multi-target/multi-hit gamma response",  "radioluminescence gamma response",    "exp.-sat. gamma response (obsolete, use gen. target/hit instead)", "linear-quadratic gamma response"}
};


///////////////////////////////////////////////////////////////////////
// RDD DATA

enum RDDModels{
  RDD_Test                 = 0,      /* no parameters */
      RDD_KatzPoint        = 1,      /* parameters: 0 - r_min [m] (lower integration limit), 1 - d_min_Gy (lower dose cut-off) */
      RDD_Geiss            = 2,      /* parameters: 0 - a0 [m] (core diameter) */
      RDD_Site             = 3,      /* parameters: 0 - a0 [m] (core diameter), 1 - d_min_Gy (lower dose cut-off)  */ // after Edmund et al., 2007, but modified with dose-cut off
      RDD_ExtTarget        = 4       /* parameters: 0 - r_min [m] (core diameter), 1 - a0 [m] (target diameter), 2 - D_min [Gy] (cut-off dose) */ //as defined in Edmund et al. , 2007
};

#define RDD_DATA_N    5

typedef struct {
  long    n;
  long    RDD_no[RDD_DATA_N];
  long	  n_parameters[RDD_DATA_N];
  char**  parameter_name[RDD_DATA_N][3];
  float   parameter_default[RDD_DATA_N][3];
  char*   RDD_name[RDD_DATA_N];
} rdd_data;

static const rdd_data AT_RDD_Data = {
		RDD_DATA_N,
    {  RDD_Test,          RDD_KatzPoint,          RDD_Geiss,        RDD_Site, RDD_ExtTarget},
    {  0, 2, 1, 2, 3},
    {  {"","",""},{"r_min_m", "d_min_Gy",""},{"a0_m","",""},{"a0_m","d_min_Gy",""},{"r_min_m","a0_m","D_min_Gy"}},
    {  {0,0,0}, {1e-10, 1e-10,0}, {5e-8,0,0}, {5e-8,1e-10,0}, {1e-10, 5e-8, 1e-10}},
    {  "Simple step test function",  "Katz' point target RDD [Katz et al., 1972]",  "Geiss' RDD [Geiss et al., 1998]",    "Site RDD, as defined in [Edmund et al., 2007]", "Katz' extended target, as defined in [Edmund et al., 2007]"}
};


enum material_no{
  Water_Liquid             = 1,
      Aluminum_Oxide       = 2,
      Aluminum             = 3,
      PMMA                 = 4
};

///////////////////////////////////////////////////////////////////////
// ER DATA

enum ERModels{
  ER_Test                  = 0,
      ER_ButtsKatz         = 1,
      ER_Waligorski        = 2,
      ER_Geiss             = 3,
      ER_Scholz            = 4
};

#define ER_DATA_N    5

typedef struct {
  long    n;
  long    ER_no[ER_DATA_N];
  char*   ER_name[ER_DATA_N];
} er_data;

static const er_data AT_ER_Data = {
		ER_DATA_N,
    {  ER_Test,          ER_ButtsKatz,          ER_Waligorski,        ER_Geiss, ER_Scholz},
    {  "simple test ER model",  "Butts & Katz' [Katz et al., 1972] ER model",  "Waligorski's ER model",    "Geiss' [Geiss, 1997] ER model", "ER_Scholz' [Scholz, 2001] ER model"}
};



void   getMaterialName(  long* material_no, char* material_name);
void   getMaterialNo(    char* material_name, long* material_no);

void   getRDDName(  long* RDD_no, char* RDD_name);
void   getRDDNo(char* RDD_name, long* RDD_no);
void   getERName(  long* ER_no, char* ER_name);
void   getGammaName(  long* Gamma_no, char* Gamma_name);
void   getMethodName(  long* Method_no, char* Method_name);


#endif // AT_CONSTANTS_H_
