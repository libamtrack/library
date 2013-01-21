/**
 * @brief Moliere theory of multiple Coulomb scattering
 */

/*
 *    AT_MultipleCoulombScattering.c
 *    ==============
 *
 *    Created on: 10.02.2012
 *    Creator: klimpki
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

#include "AT_MultipleCoulombScattering.h"


double AT_characteristic_single_scattering_angle_single( const double E_MeV_u,
														 const int    particle_charge_e,
														 const double target_thickness_cm,
														 const char   element_acronym[PARTICLE_NAME_NCHAR]){

	double	beta				=	AT_beta_from_E_single( E_MeV_u );
	double	gamma				=	AT_gamma_from_E_single( E_MeV_u );
	double	momentum_kg_m_s		= 	beta*gamma*atomic_mass_unit_MeV_c2*MeV_to_J/c_m_s;
	int		Z					=	AT_Z_from_element_acronym_single( element_acronym );
	double	electron_density	=	AT_electron_density_cm3_from_element_acronym_single( element_acronym );
	double	chi_c_2				=	4*M_PI*gsl_pow_2(particle_charge_e)*(Z+1)*electron_density*target_thickness_cm*gsl_pow_2((100*fine_structure_constant*Planck_constant_J_s)/(2*M_PI*beta*momentum_kg_m_s));

	return 	sqrt(chi_c_2);
}


int AT_characteristic_single_scattering_angle( const long  		n,
											   const double 	E_MeV_u[],
											   const int		particle_charge_e[],
											   const double		target_thickness_cm[],
											   char*			element_acronym[],
											   double        	chi_c[]){

	// loop over n to find chi_c for all particles
	long  i;
	for(i = 0; i < n; i++){
	    chi_c[i] = AT_characteristic_single_scattering_angle_single( E_MeV_u[i],
	    									particle_charge_e[i],
	    									target_thickness_cm[i],
	    									element_acronym[i] );
	  	}

	return 0;
}


double AT_screening_angle_single( const double E_MeV_u,
								  const int    particle_charge_e,
								  const char   element_acronym[PARTICLE_NAME_NCHAR]){

	double	beta				=	AT_beta_from_E_single( E_MeV_u );
	double	gamma				=	AT_gamma_from_E_single( E_MeV_u );
	double	momentum_MeV_c		= 	beta*gamma*atomic_mass_unit_MeV_c2;
	int		Z					=	AT_Z_from_element_acronym_single( element_acronym );
	double	chi_0				=	(fine_structure_constant*electron_mass_MeV_c2*pow(Z, 2.0/3.0))/(0.889*momentum_MeV_c);

	return chi_0*sqrt(1.13+3.76*gsl_pow_2(fine_structure_constant*particle_charge_e*Z/beta));
}


int AT_screening_angle( const long 		n,
						const double 	E_MeV_u[],
						const int		particle_charge_e[],
    					char*			element_acronym[],
    					double        	chi_a[]){

	// loop over n to find chi_a for all particles
	long  i;
	for(i = 0; i < n; i++){
	    chi_a[i] = AT_screening_angle_single( E_MeV_u[i],
	    									particle_charge_e[i],
	    									element_acronym[i] );
	  	}

	return 0;
}

double AT_effective_collision_number_single( const double E_MeV_u,
										     const int    particle_charge_e,
										     const double target_thickness_cm,
										     const char   element_acronym[PARTICLE_NAME_NCHAR]){

	double	chi_c	=	AT_characteristic_single_scattering_angle_single( E_MeV_u,
												particle_charge_e,
												target_thickness_cm,
												element_acronym );
	double	chi_a	=	AT_screening_angle_single( E_MeV_u,
												particle_charge_e,
												element_acronym );
	double	exp_b	=	gsl_pow_2(chi_c)/(1.167*gsl_pow_2(chi_a));

	if( exp_b < 1.14 ){
#ifndef NDEBUG
		printf("Moliere theory cannot be applied because the number of collisions in the target material is too small.");
#endif
		return	0;
	}
	else{
		return	exp_b;
	}
}


int AT_effective_collision_number( const long  		n,
								   const double 	E_MeV_u[],
								   const int		particle_charge_e[],
								   const double		target_thickness_cm[],
								   char*			element_acronym[],
								   double        	exp_b[]){

	// loop over n to find exp_b for all particles
	long  i;
	for(i = 0; i < n; i++){
	    exp_b[i] = AT_effective_collision_number_single( E_MeV_u[i],
	    									particle_charge_e[i],
	    									target_thickness_cm[i],
	    									element_acronym[i] );
	  	}

	return 0;
}


double AT_reduced_target_thickness_single( const double E_MeV_u,
										   const int    particle_charge_e,
										   const double target_thickness_cm,
										   const char   element_acronym[PARTICLE_NAME_NCHAR]){

	double	exp_b	=	AT_effective_collision_number_single( E_MeV_u,
													particle_charge_e,
													target_thickness_cm,
													element_acronym );

	return	2.6 + 2.3863*log10(1.167*exp_b) - 3.234/(log10(1.167*exp_b)+0.994);
}


int AT_reduced_target_thickness( const long  	n,
								 const double 	E_MeV_u[],
								 const int		particle_charge_e[],
								 const double	target_thickness_cm[],
								 char*			element_acronym[],
								 double        	B[]){

	// loop over n to find B for all particles
	long  i;
	for(i = 0; i < n; i++){
	    B[i] = AT_reduced_target_thickness_single( E_MeV_u[i],
	    									particle_charge_e[i],
	    									target_thickness_cm[i],
	    									element_acronym[i] );
	  	}

	return 0;
}


double AT_characteristic_multiple_scattering_angle_single( const double E_MeV_u,
														   const int    particle_charge_e,
														   const double target_thickness_cm,
														   const char   element_acronym[PARTICLE_NAME_NCHAR]){

	double	B		=	AT_reduced_target_thickness_single( E_MeV_u,
													particle_charge_e,
													target_thickness_cm,
													element_acronym );
	double	chi_c	=	AT_characteristic_single_scattering_angle_single( E_MeV_u,
													particle_charge_e,
													target_thickness_cm,
													element_acronym );

	return	chi_c*sqrt(B/2);
}


int AT_characteristic_multiple_scattering_angle( const long  	n,
												 const double 	E_MeV_u[],
												 const int		particle_charge_e[],
												 const double	target_thickness_cm[],
												 char*			element_acronym[],
												 double        	Theta_M[]){

	// loop over n to find Theta_M for all particles
	long  i;
	for(i = 0; i < n; i++){
	    Theta_M[i] = AT_characteristic_multiple_scattering_angle_single( E_MeV_u[i],
	    									particle_charge_e[i],
	    									target_thickness_cm[i],
	    									element_acronym[i] );
	  	}

	return 0;
}


double AT_Moliere_function_f0( double red_Theta){

	const double	offset1		=	0.00000;
	const double	mean1		=	0.00000;
	const double	width1		=	1.41334;
	const double	amplitude1	=	3.54605;
	double			peak1		=	offset1 + (amplitude1/(width1*sqrt(M_PI/2)))*exp(-2*gsl_pow_2((red_Theta-mean1)/width1));

	return	peak1;
}


double AT_Moliere_function_f1( double red_Theta){

	const double	offset1		=	0.00363;
	const double	mean1		=	0.05913;
	const double	width1		=	0.79315;
	const double	amplitude1	=	0.79611;
	double			peak1		=	offset1 + (amplitude1/(width1*sqrt(M_PI/2)))*exp(-2*gsl_pow_2((red_Theta-mean1)/width1));

	const double	offset2		=	0.00363;
	const double	mean2		=	0.96017;
	const double	width2		=	0.98618;
	const double	amplitude2	=  -1.04677;
	double			peak2		=	offset2 + (amplitude2/(width2*sqrt(M_PI/2)))*exp(-2*gsl_pow_2((red_Theta-mean2)/width2));

	const double	offset3		=	0.00363;
	const double	mean3		=	0.96017;
	const double	width3		=	2.15012;
	const double	amplitude3	=   0.71399;
	double			peak3		=	offset3 + (amplitude3/(width3*sqrt(M_PI/2)))*exp(-2*gsl_pow_2((red_Theta-mean3)/width3));

	return	peak1 + peak2 + peak3;
}


double AT_Moliere_function_f2( double red_Theta){

	const double	offset1		=	0.00064;
	const double	mean1		=	0.08484;
	const double	width1		=	0.85501;
	const double	amplitude1	=	4.22102;
	double			peak1		=	offset1 + (amplitude1/(width1*sqrt(M_PI/2)))*exp(-2*gsl_pow_2((red_Theta-mean1)/width1));

	const double	offset2		=	0.00064;
	const double	mean2		=	0.92230;
	const double	width2		=	0.01138;
	const double	amplitude2	=  -0.00286;
	double			peak2		=	offset2 + (amplitude2/(width2*sqrt(M_PI/2)))*exp(-2*gsl_pow_2((red_Theta-mean2)/width2));

	const double	offset3		=	0.00064;
	const double	mean3		=	0.92230;
	const double	width3		=	1.45644;
	const double	amplitude3	=  -6.11463;
	double			peak3		=	offset3 + (amplitude3/(width3*sqrt(M_PI/2)))*exp(-2*gsl_pow_2((red_Theta-mean3)/width3));

	const double	offset4		=	0.000064;
	const double	mean4		=	1.32168;
	const double	width4		=	1.08239;
	const double	amplitude4	=   3.71941;
	double			peak4		=	offset4 + (amplitude4/(width4*sqrt(M_PI/2)))*exp(-2*gsl_pow_2((red_Theta-mean4)/width4));

	return	peak1 + peak2 + peak3 + peak4;
}


double AT_scattering_angle_distribution_single( const double E_MeV_u,
												const int    particle_charge_e,
												const double target_thickness_cm,
												const char   element_acronym[PARTICLE_NAME_NCHAR],
												const double Theta){

	double	Theta_M			=	AT_characteristic_multiple_scattering_angle_single( E_MeV_u,
													particle_charge_e,
													target_thickness_cm,
													element_acronym );
	double	chi_c			=	AT_characteristic_single_scattering_angle_single( E_MeV_u,
													particle_charge_e,
													target_thickness_cm,
													element_acronym );
	double	B				=	AT_reduced_target_thickness_single( E_MeV_u,
													particle_charge_e,
													target_thickness_cm,
													element_acronym );

	double	red_Theta_pos	=	Theta/(chi_c*sqrt(B));
	double	correction0_pos	=	AT_Moliere_function_f0(red_Theta_pos);
	double	correction1_pos	=	AT_Moliere_function_f1(red_Theta_pos);
	double	correction2_pos	=	AT_Moliere_function_f2(red_Theta_pos);

	double	red_Theta_neg	=  -red_Theta_pos;
	double	correction0_neg	=	AT_Moliere_function_f0(red_Theta_neg);
	double	correction1_neg	=	AT_Moliere_function_f1(red_Theta_neg);
	double	correction2_neg	=	AT_Moliere_function_f2(red_Theta_neg);

	if(Theta>0){
		return	(1/(4*M_PI*gsl_pow_2(Theta_M)))*(correction0_pos + correction1_pos/B + correction2_pos/(B*B));
	}
	else{
		return	(1/(4*M_PI*gsl_pow_2(Theta_M)))*(correction0_neg + correction1_neg/B + correction2_neg/(B*B));
	}
}


int AT_scattering_angle_distribution( const long  	n,
									  const double 	E_MeV_u,
									  const int		particle_charge_e,
									  const double	target_thickness_cm,
									  const char	element_acronym[PARTICLE_NAME_NCHAR],
									  const double	Theta[],
									  double        distribution[]){

	// loop over n to find scattering angle distributions for all particles
	long  i;
	for(i = 0; i < n; i++){
	    distribution[i] = AT_scattering_angle_distribution_single( E_MeV_u,
	    												particle_charge_e,
	    												target_thickness_cm,
	    												element_acronym,
	    												Theta[i]);
	  	}

	return 0;
}


double AT_Highland_angle_single( const double E_MeV_u,
								 const int    particle_charge_e,
								 const double l_over_lR){

	double	beta				=	AT_beta_from_E_single( E_MeV_u );
	double	gamma				=	AT_gamma_from_E_single( E_MeV_u );
	double	momentum_MeV_c		=	beta*gamma*atomic_mass_unit_MeV_c2;

	return	(14.1/(beta*momentum_MeV_c))*particle_charge_e*sqrt(l_over_lR)*(1 + log10(l_over_lR)/9);
}


int AT_Highland_angle( const long  		n,
					   const double 	E_MeV_u[],
					   const int		particle_charge_e[],
					   const double		l_over_lR[],
    				   double        	Theta0[]){

	// loop over n to find Highland angles for all particles
	long  i;
	for(i = 0; i < n; i++){
	    Theta0[i] = AT_Highland_angle_single( E_MeV_u[i],
	    									  particle_charge_e[i],
	    									  l_over_lR[i] );
	  	}

	return 0;
}
