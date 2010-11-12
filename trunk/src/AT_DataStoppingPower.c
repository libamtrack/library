/**
 * @brief Stopping Power
 */

/*
 *    AT_DataStoppingPower.c
 *    ==============
 *
 *    Created on: 12.11.2010
 *    Creator: kongruencja
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

#include "AT_DataStoppingPower.h"


void AT_stopping_power_source_model_name_from_number( const long source_no, char* source_name){
	// TODO
}


long AT_stopping_power_source_model_number_from_name( const char* material_name ){
	// TODO
	return -1;
}


double AT_Bethe_Stopping_Number_single(	const double 	E_MeV_u,
										const long 		particle_no,
										const long 		material_no,
										const double	E_restricted_keV)
{
	  const double beta2 	= gsl_pow_2(AT_beta_from_E_single(E_MeV_u));
	  const double I_MeV	= AT_I_eV_from_material_no(material_no) / 1e6;
	  double Wm_MeV			= AT_max_relativistic_E_transfer_MeV_single(E_MeV_u);

	  /* Restricted stopping number requested? */
	  bool restricted = false;
	  if(	(E_restricted_keV > 0) &
			(E_restricted_keV / 1000.0 < Wm_MeV)){
		  restricted = true;}

	  /* First part of stopping number */
	  double SN11			=  2.0 * electron_mass_MeV_c2 * beta2 / (1.0 - beta2);
	  SN11					/= I_MeV;

	  if(	restricted){
		  Wm_MeV				= E_restricted_keV / 1000.0;
	  }
	  double SN12			= Wm_MeV / I_MeV;

	  /* Second part of stopping number */
	  double SN2			= beta2;
	  if(	restricted){
		  SN2					/= 2;
		  SN2					+= (1.0 - beta2) * Wm_MeV / (4.0 * electron_mass_MeV_c2);
	  }

	  /* Third part of stopping number (shell correction) */
	  // TODO: implement

	  /* Forth part of stopping number (density correction */
	  // TODO: implement

	  return (0.5 * log(SN11 * SN12) - SN2);
}

double AT_Bethe_Mass_Stopping_Power_MeV_cm2_g_single(	const double E_MeV_u,
														const long particle_no,
														const long material_no,
														const double E_restricted_keV)
{
	  const double beta2 		= gsl_pow_2(AT_beta_from_E_single(E_MeV_u));
	  const double Z			= AT_average_Z_from_material_no(material_no);
	  const double A			= AT_average_A_from_material_no(material_no);
	  const double z_eff		= AT_effective_charge_from_E_MeV_u_single(E_MeV_u,
																		particle_no);
	  const double k_MeV_cm2_g	= 0.307075;												// ICRU49, p.6, after Cohen and Taylor (1986)
	  const double SN			= AT_Bethe_Stopping_Number_single(	E_MeV_u,
																	particle_no,
																	material_no,
																	E_restricted_keV);

	  return( k_MeV_cm2_g * Z / A * z_eff * z_eff * SN / (beta2));
}

void AT_Bethe_Mass_Stopping_Power_MeV_cm2_g(	const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const double E_restricted_keV,
		double Mass_Stopping_Power_MeV_cm2_g[])
{
	long i;
	for (i = 0; i < n; i++){
		Mass_Stopping_Power_MeV_cm2_g[i]	=	AT_Bethe_Mass_Stopping_Power_MeV_cm2_g_single(	E_MeV_u[i],
																								particle_no[i],
																								material_no,
																								E_restricted_keV);
	}
}
