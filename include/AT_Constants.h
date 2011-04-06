#ifndef AT_CONSTANTS_H_
#define AT_CONSTANTS_H_

/**
 * @brief This files contains physical constants and unit conversions factors
 *        used by libamtrack
 */


/*
 *    AT_Constants.h
 *    ==============
 *
 *    Copyright 2006, 2010 The libamtrack team
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

static const double  proton_mass_MeV_c2      =  938.272013;                           /* proton mass CODATA http://physics.nist.gov/constants */
static const double  atomic_mass_unit_MeV_c2 = 931.494028;     /*P.J. Mohr and D.B. Newell, “Resource Letter FC-1: The Physics of Fundamental Constants,” Am. J. Phys, 78 (2010) 338 */
static const double  neutron_mass_MeV_c2      =  939.565346;                           /* proton mass */
static const double  electron_mass_MeV_c2    =  0.510998918;                          /* electron mass */
static const double  e_C                     =  1.60217653e-19;                       /* elementary charge */
static const double  e0_F_m                  =  8.8541878176e-12;                     /* electrical permitivity of the vacuum */

static const double  MeV_to_J                =  1.60217646e-13;
static const double  MeV_g_to_J_kg           =  1.60217646e-10;
static const double  m_to_cm                 =  100.0;
static const double  cm_to_mm                =  10.0;

static const double	 c_m_s				     =  299792458;

static const double  Avogadro_constant_1_mol =  6.02214179e23;

static const double  classical_electron_radius_m = 2.8179402894e-15;
static const double  Planck_constant_J_s	 =  6.62606896e-34;
static const double  Dirac_constant_J_s		 =  1.054571628e-34;
static const double	 fine_structure_constant =  7.297353e-3;


#endif /* AT_CONSTANTS_H_ */
