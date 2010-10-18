/**
 * @brief Physics related routines
 */

/*
 *    AT_PhysicsRoutines.c
 *    ==============
 *
 *    Created on: 8.01.2010
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

#include "AT_PhysicsRoutines.h"


inline double AT_beta_from_E_single( const double E_MeV_u ){ //TODO is energy per nucleon really defined like that ?
  assert( E_MeV_u > 0.);
  return sqrt(1.0 - 1.0/gsl_pow_2(1.0 + E_MeV_u/(1.0079*proton_mass_MeV_c2)));
}


int AT_beta_from_E( const long  n,
    const double  E_MeV_u[],
    double        beta[])
{
  // loop over n to find beta for all energies
  long  i;
  for(i = 0; i < n; i++){
    beta[i]        =  AT_beta_from_E_single(E_MeV_u[i]);
  }
  return 0;
}

inline double AT_gamma_from_E_single( const double E_MeV_u ){
  double beta 	= AT_beta_from_E_single( E_MeV_u );
  assert( beta < 1.0);
  return (1.0 / sqrt(1.0 - beta * beta));
}

int AT_gamma_from_E( const long  n,
    const double  E_MeV_u[],
    double        gamma[])
{
  // loop over n to find gamma for all energies
  long  i;
  for(i = 0; i < n; i++){
    gamma[i]        =  AT_gamma_from_E_single(E_MeV_u[i]);
  }
  return 0;
}


inline double AT_E_from_beta_single(  const double beta ){
  assert( beta < 1.0);
  return 1.0079*proton_mass_MeV_c2 * (sqrt(1.0 / (1.0 - gsl_pow_2(beta))) - 1.0);
}


int AT_E_from_beta(  const long  n,
    const double  beta[],
    double        E_MeV_u[])
{
  // loop over n to find E for all betas
  long  i;
  for(i = 0; i < n; i++){
    E_MeV_u[i]      =  AT_E_from_beta_single(beta[i]);
  }
  return 0;
}

inline double AT_E_MeV_u_from_momentum_single( 	const double momentum_MeV_c_u){
	double total_E_MeV_u = sqrt(momentum_MeV_c_u * momentum_MeV_c_u + 1.0079 * 1.0079 * proton_mass_MeV_c2 * proton_mass_MeV_c2);
	return (total_E_MeV_u - 1.0079 * proton_mass_MeV_c2);
}

int AT_E_MeV_u_from_momentum(  const long  n,
    const double  momentum_MeV_c_u[],
    double        E_MeV_u[])
{
  // loop over n to find E for all betas
  long  i;
  for(i = 0; i < n; i++){
    E_MeV_u[i]      =  AT_E_MeV_u_from_momentum_single(momentum_MeV_c_u[i]);
  }
  return 0;
}

inline double AT_effective_charge_from_beta_single(  const double beta,
    const long Z){
  // Return effective charge according to Barkas-Bethe-approximation
  if (Z!=1){
    return (double)(Z) * (1.0 - exp(-125.0 * beta / (pow(Z, 2.0/3.0))));
  }else{
    return 1.0 - exp(-125.0 * beta);
  }
}


int AT_effective_charge_from_beta( const long  n,
    const double  beta[],
    const long   Z[],
    double        effective_charge[])
{
  // loop over n particles
  long  i;
  for (i = 0; i < n; i++){
    effective_charge[i]    =  AT_effective_charge_from_beta_single(beta[i],Z[i]);
  }
  return 0;
}


double AT_effective_charge_from_E_MeV_u_single(  const double E_MeV_u,
    const long  particle_no){
  double beta  =  AT_beta_from_E_single(E_MeV_u);
  long Z       =  AT_Z_from_particle_no_single(particle_no);
  return AT_effective_charge_from_beta_single(beta,Z);
}


int AT_effective_charge_from_E_MeV_u( const  long  n,
    const double  E_MeV_u[],
    const long   particle_no[],
    double        effective_charge[])
{
  // loop over n particles
  long  i;
  for (i = 0; i < n; i++){
    effective_charge[i]    =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u[i],particle_no[i]);
  }
  return 0;
}


inline double AT_max_relativistic_E_transfer_MeV_single( const double E_MeV_u ){
  const double beta = AT_beta_from_E_single(E_MeV_u);
  // TODO what does it mean MeV_c2, are units correct ?
  // TODO add m/M terms (ICRU49, p. 6, eq. 2.4
  return 2.0 * electron_mass_MeV_c2 * gsl_pow_2(beta) / (1.0 - gsl_pow_2(beta));
}


inline double AT_max_classic_E_transfer_MeV_single( const double E_MeV_u ){
  assert( E_MeV_u > 0);
  return 4.0 * electron_mass_MeV_c2 / proton_mass_MeV_c2 * E_MeV_u;
}


inline double AT_max_E_transfer_MeV_single( const double E_MeV_u){
  /**
   * if E_MeV_u < 0:    use non-relativistic formula
   * if E_MeV_u > 0:    use relativistic formula
   */
  // TODO instead of using negative values of the energy switch parameter "relativistic" should be added to argument list
  if(E_MeV_u >= 0){
    return AT_max_relativistic_E_transfer_MeV_single(E_MeV_u);
  }else{
    return AT_max_classic_E_transfer_MeV_single( -1.0 * E_MeV_u);
  }
}


int AT_max_E_transfer_MeV(  const long  n,
    const double  E_MeV_u[],
    double        max_E_transfer_MeV[])
{
  // TODO instead of using negative values of the energy switch parameter "relativistic" should be added to argument list
  long  i;
  for (i = 0; i < n; i++){
    max_E_transfer_MeV[i]  =  AT_max_E_transfer_MeV_single(E_MeV_u[i]);
  }
  return 0;
}

inline double AT_momentum_from_E_MeV_c_u_single( const double E_MeV_u){
	double	beta		=	AT_beta_from_E_single(E_MeV_u);
	double	gamma		=	AT_gamma_from_E_single(E_MeV_u);
	double	m_MeV_c2	=	1.0079 * proton_mass_MeV_c2;
	return 	gamma * m_MeV_c2 * beta * 1.0;			// Here: c = 1
}

int AT_momentum_from_E_MeV_c_u( const long  n,
    const double  E_MeV_u[],
    double        momentum_MeV_c_u[])
{
  // loop over n
  long  i;
  for(i = 0; i < n; i++){
	  momentum_MeV_c_u[i]       =  AT_momentum_from_E_MeV_c_u_single(	E_MeV_u[i]);
  }
  return 0;
}

void AT_Bohr_Energy_Straggling_g_cm2(  const long*  n,
    const long*  material_no,
    double*  dsE2dz)
{ // TODO shall this function be defined for one or many materials ?
  assert( n != NULL);
  assert( *n > 0);
  long  i;
  double tmp;
  double  electron_density_m3;
  for (i = 0; i < *n; i++){

    AT_get_materials_data(  *n,
        material_no,
        NULL,
        &electron_density_m3,
        NULL, NULL, NULL, NULL, NULL, NULL);

    tmp                  =  gsl_pow_4(e_C) * electron_density_m3;
    tmp                 /=  4.0 * M_PI * gsl_pow_2(e0_F_m);
    tmp                 /=  gsl_pow_2(MeV_to_J) * m_to_cm;
    dsE2dz[i]            =  tmp;
  }
}


double AT_dose_Gy_from_fluence_cm2_single(  const double  E_MeV_u,
    const long    particle_no,
    const double  fluence_cm2,
    const long    material_no){

	double LET_MeV_cm2_g = AT_LET_MeV_cm2_g_single( E_MeV_u, particle_no, material_no );

	// Multiply by fluence, convert from MeV/g to Gy
	return LET_MeV_cm2_g * fluence_cm2 * MeV_g_to_J_kg;
}


void AT_dose_Gy_from_fluence_cm2(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    double        dose_Gy[])
{
  // Multiply by fluence, convert from MeV/g to Gy
  long  i;
  for (i = 0; i < n; i++){
    dose_Gy[i] =  AT_dose_Gy_from_fluence_cm2_single(  E_MeV_u[i], particle_no[i], fluence_cm2[i], material_no );
  }

}


double AT_fluence_cm2_from_dose_Gy_single( const double  E_MeV_u,
    const long    particle_no,
    const double  D_Gy,
    const long    material_no )
{
  double LET_MeV_cm2_g = AT_LET_MeV_cm2_g_single( E_MeV_u, particle_no, material_no );

  return (D_Gy / MeV_g_to_J_kg) / LET_MeV_cm2_g;
}


void AT_fluence_cm2(  const long  n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  D_Gy[],
    const long    material_no,
    double        fluence_cm2[])
{
  long i;
  for (i = 0; i < n; i++){
    fluence_cm2[i] =  AT_fluence_cm2_from_dose_Gy_single(  E_MeV_u[i], particle_no[i], D_Gy[i], material_no );
  }
}


void AT_interparticleDistance_m(       const long   n,
    const double  LET_MeV_cm2_g[],
    const double  fluence_cm2[],
    double        results_m[])
{
  long i;
  double fluence;
  for( i = 0 ; i < n ; i++ ){
    if( fluence_cm2[i] > 0 ){
      results_m[i] = 2.0 / sqrt(M_PI*1e4*fluence_cm2[i]);
    } else {
      fluence = (-fluence_cm2[i]) / (LET_MeV_cm2_g[i] * MeV_g_to_J_kg);
      results_m[i] = 2.0 / sqrt(M_PI*1e4*fluence);
    }
  }
}

// TODO shall it be split into separate functions, like FWHM_to_sigma_cm and so on... ?
void AT_convert_beam_parameters(  const long  n,
    double fluence_cm2[],
    double sigma_cm[],
    double N[],
    double FWHM_mm[])
{
  long  i;
  if(sigma_cm[0] == 0.0){
    for (i = 0; i < n; i++){
      sigma_cm[i]    = FWHM_mm[i] / (2.354820046 * cm_to_mm);                                // 2 * sqrt(2*ln(2))
    }
  }else{
    for (i = 0; i < n; i++){
      FWHM_mm[i]     = sigma_cm[i] * (2.354820046 * cm_to_mm);
    }
  }

  if(fluence_cm2[0] == 0.0){
    for (i = 0; i < n; i++){
      if(sigma_cm[i] != 0.0){
        fluence_cm2[i] = N[i] / (gsl_pow_2(sigma_cm[i]) * 2.0 * M_PI);
      }
    }
  }else{
    for (i = 0; i < n; i++){
      N[i]           = fluence_cm2[i] * gsl_pow_2(sigma_cm[i]) * 2.0 * M_PI;
    }
  }
}

void AT_inv_interparticleDistance_Gy(  const long   n,
    const double   LET_MeV_cm2_g[],
    const double   distance_m[],
    double         results_Gy[])
{
  long i;
  double fluence;
  for( i = 0 ; i < n ; i++ ){
    fluence = gsl_pow_2(2.0/distance_m[i]) * M_1_PI * 1e-4;
    results_Gy[i] = fluence * (LET_MeV_cm2_g[i] * MeV_g_to_J_kg);
  }
}


inline double AT_single_impact_fluence_cm2_single( const double E_MeV_u,
    const long material_no,
    const long er_model){

  double max_electron_range_m = AT_max_electron_range_m(E_MeV_u,(int)material_no,(int)er_model);
  return M_1_PI / gsl_pow_2( max_electron_range_m * m_to_cm ) ; // pi * r_max_m^2 = Track area -> single_impact_fluence [1/cm2]
}

void AT_single_impact_fluence_cm2( const long n,
    const double  E_MeV_u[],
    const long    material_no,
    const long    er_model,
    double        single_impact_fluence_cm2[])
{
  long i;
  for( i = 0 ; i < n ; i++ ){
    single_impact_fluence_cm2[i] = AT_single_impact_fluence_cm2_single(E_MeV_u[i],material_no,er_model);
  }
}


inline double AT_single_impact_dose_Gy_single( const double LET_MeV_cm2_g,
    const double single_impact_fluence_cm2){
  return LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;        // LET * fluence
}

void AT_single_impact_dose_Gy( const long n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    const long    er_model,
    double        single_impact_dose_Gy[])
{
  long i;
  for( i = 0 ; i < n ; i++ ){
    single_impact_dose_Gy[i] = AT_single_impact_dose_Gy_single(       AT_LET_MeV_cm2_g_single(  E_MeV_u[i],
                                                                                                particle_no[i],
                                                                                                material_no),
                                                                      AT_single_impact_fluence_cm2_single(      E_MeV_u[i],
                                                                                                                material_no,
                                                                                                                er_model));
  }
}

double  AT_total_D_Gy( const long  n,
    const double  E_MeV_u[],
    const long   particle_no[],
    const double  fluence_cm2[],
    const long   material_no)
{
  double   total_dose_Gy    =  0.0;
  double*  single_doses_Gy  =  (double*)calloc(n, sizeof(double));

  AT_dose_Gy_from_fluence_cm2(      n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      single_doses_Gy);

  long i;
  for (i = 0; i < n; i++){
    total_dose_Gy       += single_doses_Gy[i];
  }
  free(single_doses_Gy);

  return total_dose_Gy;
}


double AT_total_fluence_cm2( const long n,
    const double   E_MeV_u[],
    const long    particle_no[],
    const double   D_Gy[],
    const long    material_no)
{

  double*  single_fluences_cm2        =  (double*)calloc(n, sizeof(double));

  AT_fluence_cm2(      n,
      E_MeV_u,
      particle_no,
      D_Gy,
      material_no,
      single_fluences_cm2);

  double  total_fluence_cm2 = 0.0;
  long i;
  for (i = 0; i < n; i++){
    total_fluence_cm2       += single_fluences_cm2[i];
  }
  free(single_fluences_cm2);

  return total_fluence_cm2;
}


double AT_fluence_weighted_E_MeV_u( const long    n,
    const double E_MeV_u[],
    const double fluence_cm2[])
 {
  long i;

  double total_fluence_cm2 = 0.0;
  for (i = 0; i < n; i++){
    total_fluence_cm2 += fluence_cm2[i];
  }

  double average_E_MeV_u      = 0.0;
  for (i = 0; i < n; i++){
     average_E_MeV_u += fluence_cm2[i] * E_MeV_u[i];
   }

   average_E_MeV_u /= total_fluence_cm2;

   return average_E_MeV_u;
 }


double AT_dose_weighted_E_MeV_u( const long   n,
    const double  E_MeV_u[],
    const long   particle_no[],
    const double  fluence_cm2[],
    const long   material_no)
 {
  long i;

  double*  single_doses_Gy        =  (double*)calloc(n, sizeof(double));

  AT_dose_Gy_from_fluence_cm2(      n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      single_doses_Gy);

  double total_dose_Gy = 0.0;

  for (i = 0; i < n; i++){
    total_dose_Gy       += single_doses_Gy[i];
  }

  double  doseweighted_E_MeV_u      = 0.0;
  for (i = 0; i < n; i++){
     doseweighted_E_MeV_u += single_doses_Gy[i] * E_MeV_u[i];
   }

   doseweighted_E_MeV_u /= total_dose_Gy;

   free(single_doses_Gy);

   return doseweighted_E_MeV_u;
}

double AT_fluence_weighted_LET_MeV_cm2_g( const long     n,
    const double  E_MeV_u[],
    const long   particle_no[],
    const double  fluence_cm2[],
    const long   material_no)
 {
  long i;

  double*  single_LETs_MeV_cm2_g        =  (double*)calloc(n, sizeof(double));

  double total_fluence_cm2 = 0.0;
  for (i = 0; i < n; i++){
    total_fluence_cm2 += fluence_cm2[i];
  }

  AT_LET_MeV_cm2_g(  n,
      E_MeV_u,
      particle_no,
      material_no,
      single_LETs_MeV_cm2_g);

  double average_LET_MeV_cm2_g      = 0.0;
  for (i = 0; i < n; i++){
     average_LET_MeV_cm2_g += fluence_cm2[i] * single_LETs_MeV_cm2_g[i];
   }

   average_LET_MeV_cm2_g /= total_fluence_cm2;

   free(single_LETs_MeV_cm2_g);

   return average_LET_MeV_cm2_g;
}


double AT_dose_weighted_LET_MeV_cm2_g( const long  n,
    const double  E_MeV_u[],
    const long   particle_no[],
    const double  fluence_cm2[],
    const long   material_no)
 {
  long i;

  double*  single_LETs_MeV_cm2_g  =  (double*)calloc(n, sizeof(double));
  double*  single_doses_Gy        =  (double*)calloc(n, sizeof(double));

  AT_LET_MeV_cm2_g(  n,
      E_MeV_u,
      particle_no,
      material_no,
      single_LETs_MeV_cm2_g);

  AT_dose_Gy_from_fluence_cm2(      n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      single_doses_Gy);

  double total_dose_Gy = 0.0;

  for (i = 0; i < n; i++){
    total_dose_Gy       += single_doses_Gy[i];
  }

  double doseweighted_LET_MeV_cm2_g      = 0.0;
  for (i = 0; i < n; i++){
     doseweighted_LET_MeV_cm2_g += single_doses_Gy[i] * single_LETs_MeV_cm2_g[i];
   }

   doseweighted_LET_MeV_cm2_g /= total_dose_Gy;

   free(single_LETs_MeV_cm2_g);
   free(single_doses_Gy);

   return doseweighted_LET_MeV_cm2_g;
 }

double AT_fluence_weighted_stopping_power_ratio( const long     n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2[],
    const long    material_no,
    const long	  reference_material_no){

	double* LET_MeV_cm2_g			= (double*)calloc(n, sizeof(double));
	double* reference_LET_MeV_cm2_g	= (double*)calloc(n, sizeof(double));

	AT_LET_MeV_cm2_g(  n,
	    E_MeV_u,
	    particle_no,
	    material_no,
	    LET_MeV_cm2_g);

	AT_LET_MeV_cm2_g(  n,
	    E_MeV_u,
	    particle_no,
	    reference_material_no,
	    reference_LET_MeV_cm2_g);

	long i;
	double	stopping_power_ratio	= 0.0;
	double	total_fluence_cm2		= 0.0;

	for (i = 0; i < n; i++){
		if(reference_LET_MeV_cm2_g[i] != 0){
			stopping_power_ratio	+=	fluence_cm2[i] * LET_MeV_cm2_g[i] / reference_LET_MeV_cm2_g[i];
		}
		total_fluence_cm2		+= fluence_cm2[i];
	}

	free(LET_MeV_cm2_g);
	free(reference_LET_MeV_cm2_g);

	return(stopping_power_ratio / total_fluence_cm2);
}


double AT_total_u(    const long n,
                const double E_MeV_u[],
                const long particle_no[],
                const double fluence_cm2[],
                const long material_no,
                const long er_model)
{
  double* norm_fluence    =  (double*)calloc(n, sizeof(double));
  double total_D_Gy       =  AT_total_D_Gy( n, E_MeV_u, particle_no, fluence_cm2, material_no);

  AT_normalize( n, fluence_cm2, norm_fluence);

  double u                =  0.0;
  long i;
  for (i = 0; i < n; i++){
    double single_impact_fluence_cm2 =  AT_single_impact_fluence_cm2_single(  E_MeV_u[i], material_no, er_model);
    double LET_MeV_cm2_g             = AT_LET_MeV_cm2_g_single(  E_MeV_u[i], particle_no[i], material_no);
    u += norm_fluence[i] * AT_single_impact_dose_Gy_single(  LET_MeV_cm2_g, single_impact_fluence_cm2 );
  }

  free(norm_fluence);
  return(total_D_Gy / u);
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
