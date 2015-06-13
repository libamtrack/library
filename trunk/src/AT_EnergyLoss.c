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


void AT_kappa_multi( const long n,
		const double 	E_MeV_u[],
		const long      particle_no[],
		const long 		material_no,
		const double    slab_thickness_um[],
		double          kappa[])
{
	long i;
	for(i = 0; i < n; i++){
		kappa[i] = AT_kappa_single(E_MeV_u[i], particle_no[i], material_no, slab_thickness_um[i]);
	}
}


void AT_lambda_mean_multi( const long n,
		const double	E_MeV_u[],
		const long      particle_no[],
		const long 		material_no,
		const double    slab_thickness_um[],
		double 			lambda_mean[])
{
	long i;
	for(i = 0; i < n; i++){
		lambda_mean[i] = AT_lambda_mean_single(E_MeV_u[i], particle_no[i], material_no, slab_thickness_um[i]);
	}

}

double AT_lambda_mean_single( const double	E_MeV_u,
		const long      particle_no,
		const long 		material_no,
		const double    slab_thickness_um)
{
	double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double 	beta          = AT_beta_from_E_single(E_MeV_u);
	return(-0.42278433509 - beta*beta - log(kappa));
}

void AT_lambda_max_multi( const long n,
		const double	E_MeV_u[],
		const long      particle_no[],
		const long 		material_no,
		const double    slab_thickness_um[],
		double 			lambda_max[])
{
	long i;
	for(i = 0; i < n; i++){
		double lambda_mean = AT_lambda_mean_single(E_MeV_u[i], particle_no[i], material_no, slab_thickness_um[i]);
		lambda_max[i] = AT_lambda_max_single( lambda_mean);
	}
}

double AT_lambda_max_single( double lambda_mean)
{
	return(0.60715 + 1.1934 * lambda_mean + (0.67794 + 0.052382 * lambda_mean)*exp(0.94753 + 0.74442*lambda_mean));
}


double AT_kappa_single( const double 	E_MeV_u,
		const long      particle_no,
		const long 		material_no,
		const double    slab_thickness_um)
{
	// Factor 1/2 results from the difference between ICRU49, where
	// dE/dx = K * ... * (0.5 ln() - beta² - delta/2 - ...) and
	// Seltzer & Berger where
	// dE/dx = 1/2 * K * ... * (ln() - 2*beta² - delta - ...)
	return( 0.5 * AT_el_energy_loss_leading_term_MeV_cm2_g(E_MeV_u, particle_no, material_no, false) *
			AT_density_g_cm3_from_material_no(material_no) *
			(slab_thickness_um / 1e4) /
			AT_max_E_transfer_MeV_single(E_MeV_u));
}


double AT_lambda_landau_from_energy_loss_single( const double energy_loss_keV,
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double 	beta          = AT_beta_from_E_single(E_MeV_u);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	double  mean_loss_MeV = AT_Bethe_mean_energy_loss_MeV( E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);

	return( (energy_loss_keV / 1000.0 - mean_loss_MeV) / xi - 0.42278433509 - beta*beta - log(kappa));
}

void AT_lambda_landau_from_energy_loss_multi( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double lambda_landau[]){

	for(long i = 0; i < n; i++){
		lambda_landau[i] = AT_lambda_landau_from_energy_loss_single( energy_loss_keV[i], E_MeV_u, particle_no, material_no, slab_thickness_um);
	}

	return;
}


double AT_lambda_vavilov_from_energy_loss_single( const double energy_loss_keV,
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double 	beta          = AT_beta_from_E_single(E_MeV_u);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	double  mean_loss_MeV = AT_Bethe_mean_energy_loss_MeV( E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);

	return( kappa*((energy_loss_keV / 1000.0 - mean_loss_MeV) / xi - 0.42278433509 - beta*beta) );
}

void AT_lambda_vavilov_from_energy_loss_multi( const long n,
		const double energy_loss_keV[],
		const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um,
		double lambda_vavilov[]){

	for(long i = 0; i < n; i++){
		lambda_vavilov[i] =  AT_lambda_vavilov_from_energy_loss_single( energy_loss_keV[i],
				E_MeV_u,
				particle_no,
				material_no,
				slab_thickness_um);
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

	double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
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

	double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
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

	double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
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
double AT_lambda_Vavilov_Mode( const double kappa, const double beta){

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

double AT_lambda_Vavilov_FWHM_left( const double kappa, const double beta){

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

double AT_lambda_Vavilov_FWHM_right( const double kappa, const double beta){

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

double AT_lambda_Vavilov_FWHM( const double kappa, const double beta){
	return(AT_lambda_Vavilov_FWHM_right(kappa, beta) - AT_lambda_Landau_FWHM_left(kappa, beta));
}

double AT_energy_loss_keV_Vavilov_FWHM( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
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

double AT_energy_loss_keV_Vavilov_Mode( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double beta  = AT_beta_from_E_single(E_MeV_u);

	return AT_energy_loss_from_lambda_single(AT_lambda_Vavilov_Mode(kappa, beta),
			E_MeV_u,
			particle_no,
			material_no,
			slab_thickness_um);
}

double AT_lambda_Vavilov_Mean( const double kappa, const double beta) {
   return -0.42278433509 -1.0 - log(kappa) - beta*beta;
}

double AT_lambda_Vavilov_Variance( const double kappa, const double beta) {
   return (1.0 - 0.5*beta*beta) / kappa;
}

double AT_Vavilov_FWHM(){
	return 0.0;
}

double AT_lambda_Vavilov_Skewness( const double kappa, const double beta) {
   return (0.5-(beta*beta)/3)/(kappa*kappa) * pow ((1-0.5*beta*beta)/kappa, -1.5);
}

double AT_lambda_Landau_Mode(){
	return -0.2258;
}

double AT_energy_loss_keV_Landau_Mode( const double E_MeV_u,
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

double AT_energy_loss_keV_Landau_FWHM( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){
	double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
	double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
	return(4.02 * xi * 1000);
}

double AT_lambda_Landau_Mean( const double kappa, const double beta) {
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

double AT_energy_loss_mode( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
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

double AT_energy_loss_FWHM( const double E_MeV_u,
		const long particle_no,
		const long material_no,
		const double slab_thickness_um){

	double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
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
