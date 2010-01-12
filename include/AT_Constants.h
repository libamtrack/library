#ifndef AT_CONSTANTS_H_
#define AT_CONSTANTS_H_

/**
 * @file
 * @brief ...
 */


/*
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
  ME_Grid      = 1,
      ME_SPIFF = 2,
      ME_Katz  = 3
};


#endif /* AT_CONSTANTS_H_ */
