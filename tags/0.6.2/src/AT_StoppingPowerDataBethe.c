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

#include "AT_StoppingPowerDataBethe.h"

int AT_Bethe_wrapper( const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const char info[],
		double mass_stopping_power_MeV_cm2_g[]){
	AT_Bethe_energy_loss_MeV_cm2_g(n,
			E_MeV_u,
			particle_no,
			material_no,
			-1.0,
			true,
			mass_stopping_power_MeV_cm2_g);
	return AT_Success;
}

double AT_el_energy_loss_leading_term_MeV_cm2_g(	const double 	E_MeV_u,
						const long 		particle_no,
						const long 		material_no,
						const bool		use_effective_charge){

	const double Z			= AT_average_Z_from_material_no(material_no);
	const double A			= AT_average_A_from_material_no(material_no);
	assert( A > 0 );

	const double beta2 		= gsl_pow_2(AT_beta_from_E_single(E_MeV_u));
	assert( beta2 > 0);

	double z;
	if(use_effective_charge){
		z = AT_effective_charge_from_E_MeV_u_single(E_MeV_u, particle_no);
	}else{
		z = AT_Z_from_particle_no_single(particle_no);
	}

	// ICRU49, p.6, after Cohen and Taylor (1986), k_MeV_cm2_g	= 0.307075
	return(0.307075 * (Z / A) * (z*z) / (beta2));

}

double AT_Bethe_Stopping_Number(	const double 	E_MeV_u,
									const long 		particle_no,
									const long 		material_no,
									const double	E_restricted_keV)
{
	  const double beta2 	= gsl_pow_2(AT_beta_from_E_single(E_MeV_u));
	  const double I_eV		= AT_I_eV_from_material_no(material_no);
	  const double I_MeV	= I_eV * 1e-6;
	  double Wm_MeV			= AT_max_relativistic_E_transfer_MeV_single(E_MeV_u);

	  /* Restricted stopping number requested? */
	  bool restricted = false;
	  if(	(E_restricted_keV > 0.0) && (E_restricted_keV / 1000.0 < Wm_MeV))
		  restricted = true;

	  /* First part of stopping number */
	  double SN11			=  2.0 * electron_mass_MeV_c2 * beta2 / (1.0 - beta2);
	  assert( I_MeV > 0. );
	  SN11					/= I_MeV;

	  if(	restricted){
		  Wm_MeV				= E_restricted_keV * 1e-3;
	  }
	  double SN12			= Wm_MeV / I_MeV;

	  /* Second part of stopping number */
	  double SN2			= beta2;
	  if(	restricted){
		  SN2					/= 2;
		  SN2					+= (1.0 - beta2) * Wm_MeV / (4.0 * electron_mass_MeV_c2);
	  }

	  /* Third part of stopping number (density correction following Sternheimer, 1971) */
	  double delta				= 0.0;
	  long   phase              = AT_phase_from_material_no(material_no);
	  if( (phase =! phase_undefined) ){
		  double kinetic_variable	= AT_kinetic_variable_single(E_MeV_u);
		  double plasma_energy_J 	= AT_plasma_energy_J_from_material_no(material_no);
		  double I_J				= I_MeV * MeV_to_J;

		  double C					= 1.0 + 2.0 * log(I_J / plasma_energy_J);

		  // Find x_0 and x_1 dependent on phase, I-value and C
		  double x_0 = 0.0;
		  double x_1 = 0.0;
		  if( phase == phase_condensed){
			  if(I_eV < 100){
				  x_1	= 2.0;
				  if(C <= 3.681){
					  x_0 	= 0.2;
				  }else{
					  x_0	= 0.326 * C - 1.0;
				  }
			  }else{ // I_eV >= 100
				  x_1	= 3.0;
				  if(C <= 5.215){
					  x_0	= 0.2;
				  }else{
					  x_0	= 0.326 * C - 1.5;
				  }
			  }
		  }else{ // gaseous material
			  x_0	= 0.326 * C - 2.5;
			  x_1	= 5.0;
			  if(C < 10.0){
				  x_0	= 1.6;
				  x_1	= 4.0;
			  }
			  if(C >= 10.0 && C < 10.5){
				  x_0	= 1.7;
				  x_1	= 4.0;
			  }
			  if(C >= 10.5 && C < 11.0){
				  x_0	= 1.8;
				  x_1	= 4.0;
			  }
			  if(C >= 11.0 && C < 11.5){
				  x_0	= 1.9;
				  x_1	= 4.0;
			  }
			  if(C >= 11.5 && C < 12.25){
				  x_0	= 2.0;
				  x_1	= 4.0;
			  }
			  if(C >= 12.25 && C < 13.804){
				  x_0	= 2.0;
				  x_1	= 5.0;
			  }
		  }

		  double x_a				= C / 4.606;
		  double m					= 3.0;
		  double a					= 4.606 * (x_a - x_0) / pow(x_1 - x_0, m);

		  if( kinetic_variable >= x_0 && kinetic_variable <= x_1){
			  delta						= 4.606 * kinetic_variable - C + a * pow(x_1 - kinetic_variable, m);
		  }
		  if( kinetic_variable > x_1){
			  delta						= 4.606 * kinetic_variable - C;
		  }
	  }
	  double SN3		= delta;

	  /* Forth part of stopping number (shell correction) TODO: implement */


	  assert( SN11 > 0. );
	  assert( SN12 > 0. );

	  return (0.5 * log(SN11 * SN12) - SN2 - SN3);
}


double AT_Bethe_energy_loss_MeV_cm2_g_single(	const double E_MeV_u,
												const long particle_no,
												const long material_no,
												const double E_restricted_keV,
												const bool use_effective_charge)
{
	double result = 0.0;
	/* Compute only above 1.0 MeV, otherwise theory is too wrong
	 * below return zero */
	// TODO: Find smarter criterion because this may cause problems in the code (as it did
	// TODO: with the inappropriatly set lower limit for CSDA range integration (was 0, now 1.0 MeV)
	if(E_MeV_u >= BETHE_LOWER_LIMIT_E_MEV_U){
		result =  AT_el_energy_loss_leading_term_MeV_cm2_g(	E_MeV_u, particle_no, material_no, use_effective_charge) \
				* AT_Bethe_Stopping_Number(	E_MeV_u, particle_no, material_no, E_restricted_keV);
	}
    return result;

}

void AT_Bethe_energy_loss_MeV_cm2_g(	const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const double E_restricted_keV,
		const bool use_effective_charge,
		double Bethe_energy_loss_MeV_cm2_g[])
{
	long i;
	for (i = 0; i < n; i++){
		Bethe_energy_loss_MeV_cm2_g[i]	=	AT_Bethe_energy_loss_MeV_cm2_g_single(	E_MeV_u[i],
																					particle_no[i],
																					material_no,
																					E_restricted_keV,
																					use_effective_charge);
	}
}


