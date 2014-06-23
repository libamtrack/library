/*
 *    AT_EnergyLoss.c
 *    ===============
 *
 *    Copyright 2006, 2014 The libamtrack team
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

#include "AT_EnergyLoss.h"


#ifdef HAVE_CERNLIB

#include "cfortran.h"
#include "gen.h"

void AT_Vavilov_PDF(const long n, const double lambda_V[], const double kappa, const double beta,
		double density[]){
	VAVSET(kappa, beta, 1);

	for(int i = 0; i < n; i++){
		density[i] = VAVDEN(lambda_V[i]);
	}
}

void AT_Landau_PDF(const long n, const double lambda[], double density[]){
	for(int i = 0; i < n; i++){
		density[i] = DENLAN(lambda[i]);
	}
}

#endif /* HAVE_CERNLIB */


int AT_Rutherford_SDCS(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const long n,
		const double T_MeV[],
		double dsdT_m2_MeV[]){

	// Get particle data
	long I 			= AT_nuclear_spin_from_particle_no_single(particle_no);
	double Z_eff  	= AT_effective_charge_from_E_MeV_u_single(E_MeV_u, particle_no);
	long Z			= AT_Z_from_particle_no_single(particle_no);
	long A			= AT_A_from_particle_no_single(particle_no);

	// Get material Z
	double Z_material = AT_average_Z_from_material_no(material_no);

	// We follow the formulation of Sawakuchi, 2007
	double particle_mass_MeV_c2	= Z * proton_mass_MeV_c2 + (A-Z) * neutron_mass_MeV_c2;
	double total_E_MeV_u		= particle_mass_MeV_c2 / A + E_MeV_u;
	double Q_c_MeV_c2			= particle_mass_MeV_c2 * particle_mass_MeV_c2 / electron_mass_MeV_c2;
	double K_MeV_m2  			= 2 * M_PI * pow(classical_electron_radius_m, 2) * electron_mass_MeV_c2;

	double beta				= AT_beta_from_E_single(E_MeV_u);
	double beta2			= beta * beta;
	double T_max_MeV		= AT_max_E_transfer_MeV_single(E_MeV_u);

	double term_0, term_1, term_2;

	long i;
	for (i = 0; i < n; i++){
		if(T_MeV[i] > T_max_MeV){
			dsdT_m2_MeV[i]		= 0.0;
		}else{
			term_0				= K_MeV_m2 * Z_material * Z_eff * Z_eff / (beta2 * T_MeV[i] * T_MeV[i]);
			term_1				= 1.0 - beta2 * T_MeV[i] / T_max_MeV;
			if(I == 0.0){
				dsdT_m2_MeV[i]			= term_0 * term_1;
			}
			if(I == 0.5){
				term_2				= T_MeV[i] * T_MeV[i] / (2.0 * total_E_MeV_u * total_E_MeV_u);
				dsdT_m2_MeV[i]			= term_0 * (term_1 + term_2);
			}
			if(I == 1.0){
				term_2				=  (1.0 + T_MeV[i] / (2.0 * Q_c_MeV_c2));
				term_2				*= T_MeV[i] * T_MeV[i] / (3.0 * total_E_MeV_u * total_E_MeV_u);
				term_2				+= term_1 * (1.0 + T_MeV[i] / (3.0 * Q_c_MeV_c2));
				dsdT_m2_MeV[i]			= term_0 * term_2;
			}
		}
	}

	return EXIT_SUCCESS;
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


double AT_Bethe_mean_energy_loss_MeV(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){
	  return( AT_Bethe_energy_loss_MeV_cm2_g_single(E_MeV_u,
				particle_no,
				material_no,
				-1.0,
				false) *
				AT_density_g_cm3_from_material_no(material_no) *
				(slab_thickness_um / 1e4));
}

double AT_kappa( const double 	E_MeV_u,
		const long      particle_no,
		const long 		material_no,
		const double    slab_thickness_um){
	// Factor 1/2 results from the difference between ICRU49, where
	// dE/dx = K * ... * (0.5 ln() - beta² - delta/2 - ...) and
	// Seltzer & Berger where
	// dE/dx = 1/2 * K * ... * (ln() - 2*beta² - delta - ...)
	return( 0.5 * AT_el_energy_loss_leading_term_MeV_cm2_g(E_MeV_u, particle_no, material_no, false) *
			AT_density_g_cm3_from_material_no(material_no) *
			(slab_thickness_um / 1e4) /
			AT_max_E_transfer_MeV_single(E_MeV_u));
}

#ifdef HAVE_CERNLIB
void AT_lambda_from_energy_loss( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double lambda_V[]){

	double 	kappa	      = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double 	beta          = AT_beta_from_E_single(E_MeV_u);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	double  mean_loss_MeV = AT_Bethe_mean_energy_loss_MeV( E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);

	for(long i = 0; i < n; i++){
		lambda_V[i] =  (energy_loss_keV[i] / 1000.0 - mean_loss_MeV) / xi - 0.42278433509 - beta*beta - log(kappa);
	}

	return;
}

double AT_energy_loss_from_lambda_single( const double lambda,
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double 	kappa	      = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double 	beta          = AT_beta_from_E_single(E_MeV_u);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	double  mean_loss_MeV = AT_Bethe_mean_energy_loss_MeV( E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);

	return(1000.0 * (xi * (lambda + 0.42278433509 + beta*beta + log(kappa)) + mean_loss_MeV));
}

void AT_energy_loss_from_lambda( const long n,
		const double lambda[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double energy_loss_keV[]){

	double 	kappa	      = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double 	beta          = AT_beta_from_E_single(E_MeV_u);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	double  mean_loss_MeV = AT_Bethe_mean_energy_loss_MeV( E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);

	for(long i = 0; i < n; i++){
		energy_loss_keV[i] =  1000.0 * (xi * (lambda[i] + 0.42278433509 + beta*beta + log(kappa)) + mean_loss_MeV);
	}

	return;
}

void AT_Landau_energy_loss_distribution( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double fDdD[]){

	double 	kappa	      = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	double* lambda        = (double*)calloc(n, sizeof(double));

	AT_lambda_from_energy_loss( n,
			energy_loss_keV,
			E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um,
			lambda);

	AT_Landau_PDF( n,
			lambda,
			fDdD);

	for(int i = 0; i < n; i++){
		fDdD[i] /= xi;
	}

	return;
}

void AT_Vavilov_energy_loss_distribution( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double fDdD[]){

	double 	kappa	      = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double 	beta          = AT_beta_from_E_single(E_MeV_u);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	double* lambda        = (double*)calloc(n, sizeof(double));

	AT_lambda_from_energy_loss( n,
			energy_loss_keV,
			E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um,
			lambda);

	AT_Vavilov_PDF( n,
			lambda,
		    kappa,
		    beta,
			fDdD);

	for(int i = 0; i < n; i++){
		fDdD[i] /= xi;
	}

	return;
}

void AT_energy_loss_distribution( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double fDdD[]){

	double 	kappa	      = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	if (kappa <= 0.01){
		AT_Landau_energy_loss_distribution( n,
				energy_loss_keV,
				E_MeV_u,
				particle_no,
				material_no,
				slab_thickness_um,
				fDdD);
		return;
	}
	if (kappa < 10){
		AT_Vavilov_energy_loss_distribution( n,
				energy_loss_keV,
				E_MeV_u,
				particle_no,
				material_no,
				slab_thickness_um,
				fDdD);
		return;
	}
}
/**
 * This is code from ROOT (CERN) 5.3.1
 */
double AT_lambda_Vavilov_Mode(const double kappa, const double beta){

   VAVSET(kappa, beta, 1);
   double x = -4.22784335098467134e-01 - log(kappa) - beta*beta;
   if (x>-0.223172) x = -0.223172;
   double eps = 0.01;
   double dx;

   do {
      double p0 = VAVDEN(x - eps);
      double p1 = VAVDEN(x);
      double p2 = VAVDEN(x + eps);
      double y1 = 0.5*(p2-p0)/eps;
      double y2 = (p2-2*p1+p0)/(eps*eps);
      if(y2!=0){
    	  dx = - y1/y2;
      }else{
    	  dx = 0.0;
      }
      x += dx;
      if (fabs(dx) < eps) eps = 0.1*fabs(dx);
   } while (fabs(dx) > 1E-5);
   return x;
}

double AT_lambda_Vavilov_FWHM_left(const double kappa, const double beta){

   double x = AT_lambda_Vavilov_Mode(kappa, beta);
   VAVSET(kappa, beta, 1);
   double p = VAVDEN(x) * 0.5;

   x -= 1.3637;
   double eps = 0.01;
   double dx;

   do {
      double p0 = VAVDEN(x);
      double p1 = VAVDEN(x - eps);
      double p2 = VAVDEN(x + eps);
      double y1 = p0 - p;
      double y2 = 0.5*(p2-p1)/eps;
      if(y2!=0){
    	  dx = - y1/y2;
      }else{
    	  dx = 0.0;
      }
      x += dx;
      if (fabs(dx) < eps) eps = 0.1*fabs(dx);
   } while (fabs(dx) > 1E-5);
   return x;
}

double AT_lambda_Vavilov_FWHM_right(const double kappa, const double beta){

   double x = AT_lambda_Vavilov_Mode(kappa, beta);
   VAVSET(kappa, beta, 1);
   double p = VAVDEN(x) * 0.5;

   x += 2.655;
   double eps = 0.01;
   double dx;

   do {
      double p0 = VAVDEN(x);
      double p1 = VAVDEN(x - eps);
      double p2 = VAVDEN(x + eps);
      double y1 = p0 - p;
      double y2 = 0.5*(p2-p1)/eps;
      if(y2!=0){
    	  dx = - y1/y2;
      }else{
    	  dx = 0.0;
      }
      x += dx;
      if (fabs(dx) < eps) eps = 0.1*fabs(dx);
   } while (fabs(dx) > 1E-5);
   return x;
}

double AT_lambda_Vavilov_FWHM(const double kappa, const double beta){
	return(AT_lambda_Vavilov_FWHM_right(kappa, beta) - AT_lambda_Landau_FWHM_left(kappa, beta));
}

double AT_energy_loss_keV_Vavilov_FWHM(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double kappa = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double beta  = AT_beta_from_E_single(E_MeV_u);

	return AT_energy_loss_from_lambda_single(AT_lambda_Vavilov_FWHM_right(kappa, beta),
			E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um) -
			AT_energy_loss_from_lambda_single(AT_lambda_Vavilov_FWHM_left(kappa, beta),
						E_MeV_u,
						particle_no,
						material_no,
						slab_thickness_um);
}

double AT_energy_loss_keV_Vavilov_Mode(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double kappa = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double beta  = AT_beta_from_E_single(E_MeV_u);

	return AT_energy_loss_from_lambda_single(AT_lambda_Vavilov_Mode(kappa, beta),
			E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);
}

double AT_lambda_Vavilov_Mean(const double kappa, const double beta) {
   return -0.42278433509 -1.0 - log(kappa) - beta*beta;
}

double AT_lambda_Vavilov_Variance(const double kappa, const double beta) {
   return (1.0 - 0.5*beta*beta) / kappa;
}

double AT_Vavilov_FWHM(){
	return 0.0;
}

double AT_lambda_Vavilov_Skewness(const double kappa, const double beta) {
   return (0.5-(beta*beta)/3)/(kappa*kappa) * pow ((1-0.5*beta*beta)/kappa, -1.5);
}

double AT_lambda_Landau_Mode(){
	return -0.2258;
}

double AT_energy_loss_keV_Landau_Mode(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){
	return AT_energy_loss_from_lambda_single(AT_lambda_Landau_Mode(),
			E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);
}

double AT_lambda_Landau_FWHM_left(){
	return -0.2258-1.3637;
}

double AT_lambda_Landau_FWHM_right(){
	return -0.2258+2.655;
}

double AT_lambda_Landau_FWHM(){
	return(AT_lambda_Landau_FWHM_right() - AT_lambda_Landau_FWHM_left());
}

double AT_energy_loss_keV_Landau_FWHM(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){
	double kappa = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	return(4.02 * xi * 1000);
}

double AT_lambda_Landau_Mean(const double kappa, const double beta) {
   return -0.42278433509 - log(kappa) - beta*beta;
}



double AT_Gauss_Mode(){
	return 0;
}

double AT_Gauss_Mean(){
	return 0;
}

double AT_Gauss_FWHM(){
	return 0;
}

double AT_energy_loss_mode(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double 	kappa	      = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	if (kappa <= 0.01){
		return AT_energy_loss_keV_Landau_Mode(E_MeV_u,
				particle_no,
				material_no,
				slab_thickness_um);
	}
	if (kappa < 10){
		return AT_energy_loss_keV_Vavilov_Mode(E_MeV_u,
				particle_no,
				material_no,
				slab_thickness_um);
	}
	// kappa >= 10
	return 	AT_Bethe_mean_energy_loss_MeV( E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);
}

double AT_energy_loss_FWHM(const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double 	kappa	      = AT_kappa(E_MeV_u, particle_no, material_no, slab_thickness_um);
	if (kappa <= 0.01){
		return AT_energy_loss_keV_Landau_FWHM(E_MeV_u,
				particle_no,
				material_no,
				slab_thickness_um);
	}
	if (kappa < 10){
		return AT_energy_loss_keV_Vavilov_FWHM(E_MeV_u,
				particle_no,
				material_no,
				slab_thickness_um);
	}
//	// kappa >= 10
//	return 	AT_Bethe_mean_energy_loss_MeV( E_MeV_u,
//			particle_no,
//			material_no,
//			slab_thickness_um);
	return 0.0;
}

#endif /* HAVE_CERNLIB */
