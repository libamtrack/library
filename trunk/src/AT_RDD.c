/**
 * @file
 * @brief Radial Dose Distribution models
 */

/*
 *    AT_RDD.c
 *    ========
 *
 *    Created on: 28.07.2009
 *    Author: greilich
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

#include "AT_RDD.h"


void getRDDName(const long* RDD_no, char* RDD_name){
  strcpy(RDD_name,"*** invalid choice ***");
  long i;
  for (i = 0; i < RDD_DATA_N; i++){
    if (AT_RDD_Data.RDD_no[i] == *RDD_no){
      strcpy(RDD_name, AT_RDD_Data.RDD_name[i]);
    }
  }
}


void getRDDNo(const char* RDD_name, long* RDD_no){
  *RDD_no = 0;
  long i;
  for (i = 0; i < RDD_DATA_N; i++){
    if (strcmp(RDD_name, AT_RDD_Data.RDD_name[i]) == 0){
      *RDD_no = AT_RDD_Data.RDD_no[i];
      break;
    }
  }
}

//////////////////////////////////////////// Katz Model Coefficients /////////////////////////////////////////////

inline double AT_RDD_Katz_C_J_m( const double electron_density_m3){
  // C = 2 pi N e^4 / ( m c^2 (4 pi eps_0)^2)
  // C = N e^4 / ( 8 m pi (c * eps_0)^2)
  return electron_density_m3 * gsl_pow_4(GSL_CONST_MKSA_ELECTRON_CHARGE) * M_1_PI * 0.125 / (GSL_CONST_MKSA_MASS_ELECTRON * gsl_pow_2(GSL_CONST_MKSA_SPEED_OF_LIGHT * GSL_CONST_MKSA_VACUUM_PERMITTIVITY));
}


inline double AT_RDD_Katz_coeff_Gy(const double C_J_m,
    const double Z_eff,
    const double beta,
    const double density_kg_m3,
    const double r_max_m){

  // (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
  return C_J_m * 0.5 * M_1_PI * gsl_pow_2(Z_eff/beta)/ (density_kg_m3 * gsl_pow_2(r_max_m));
}

//////////////////////////////////////// Linear ER model calculations /////////////////////////////////////////////


inline double   AT_RDD_Katz_LinearER_Dpoint_Gy(        const double r_m,
    const double r_max_m,
    const double Katz_point_coeff_Gy){

  // D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/rmax^2 * rmax/r * (rmax/r - 1.)
  return Katz_point_coeff_Gy * r_max_m/r_m * (r_max_m/r_m - 1.0);
}


inline double   AT_RDD_Katz_LinearER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double r_max_m,
    const double Katz_point_coeff_Gy){

  // Dav(r1,r2) = 2 * coeff * ( log(r2/r1) - (r2 - r1)/rmax ) / ((r2/rmax)^2 - (r1/rmax)^2)
  return 2.0 * Katz_point_coeff_Gy * ( log(r2_m/r1_m) - (r2_m - r1_m)/r_max_m ) / (gsl_pow_2(r2_m/r_max_m) - gsl_pow_2(r1_m/r_max_m));
}


double   AT_RDD_Katz_LinearER_dEdx_J_m(  const double a0_m,
    const double r_max_m,
    const double material_density_kg_m3,
    const double Katz_point_coeff_Gy){

// dEdx = rho \int_a0^rmax r D(r) dr =
//      = rho * (pi rmax^2 - pi a0^2) * D_av(a0,rmax)

  return material_density_kg_m3 * M_PI * \
          (gsl_pow_2(r_max_m) - gsl_pow_2(a0_m)) * \
          AT_RDD_Katz_LinearER_Daverage_Gy(a0_m, r_max_m, r_max_m, Katz_point_coeff_Gy);
}


double   AT_RDD_Katz_LinearER_DSite_Gy( const double r_m,
    const double a0_m,
    const double r_max_m,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy){

  //Dsite(r) = 1 / (rho pi a0^2) * (LET - dEdx)   for r < a0
  //Dsite(r) = D(r)                               for r >= a0
  if( r_m < a0_m ){
    return M_1_PI * (LET_J_m - dEdx_J_m)/ (material_density_kg_m3 * gsl_pow_2(a0_m));
  } else {
    return AT_RDD_Katz_LinearER_Dpoint_Gy(r_m, r_max_m, Katz_point_coeff_Gy);
  }

}


//////////////////////////////////////// Power law ER model calculations /////////////////////////////////////////////


inline double AT_RDD_Katz_PowerLawER_Dpoint_Gy(const double r_m,
    const double alpha,
    const double r_max_m,
    const double Katz_point_coeff_Gy){

  // D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/r^2 * 1/alpha * (1 - r/rmax)^(1/alpha)
  return Katz_point_coeff_Gy * (1.0/alpha) * gsl_pow_2(r_max_m/r_m) * pow(1.0 - r_m/r_max_m, 1.0 / (alpha));
}


double   AT_RDD_Katz_PowerLawER_DaverageKernel(  const double x1,
    const double x2,
    const double alpha){

  // C1 = (1-x1)^(1/alpha) ((x1-1)/x1)^(-1/alpha)
  // HGF1 =  _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x1)
  // C2 = (1-x2)^(1/alpha) ((x2-1)/x2)^(-1/alpha)
  // HGF2 =  _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x2)
  // kernel =  2/(x2^2 - x1^2) * (C2* HGF2 - C1 * HGF1)
  const double alpha_inv = 1.0/alpha;
  const double C1 = pow( 1.0-x1, alpha_inv ) * pow( (x1-1.0)/x1, -alpha_inv );
  const double HGF1 = gsl_sf_hyperg_2F1( -alpha_inv, -alpha_inv, 1.0-alpha_inv, 1.0/x1 );
  const double C2 = pow( 1.0-x2, alpha_inv ) * pow( (x2-1.0)/x2, -alpha_inv );
  const double HGF2 = gsl_sf_hyperg_2F1( -alpha_inv, -alpha_inv, 1.0-alpha_inv, 1.0/x2 );
  return 2.0 * (C2 * HGF2 - C1 * HGF1) / (gsl_pow_2(x2) - gsl_pow_2(x1));
}


double   AT_RDD_Katz_PowerLawER_DaverageKernel_approx(  const double x1,
    const double x2,
    const double alpha){

   // F1 = x1 / alpha^2 ( (x1 / 4alpha) * (1/alpha - 1) - 1 ) + log(x1)
   // F2 = x2 / alpha^2 ( (x2 / 4alpha) * (1/alpha - 1) - 1 ) + log(x2)
   // kernel =  2/(x2^2 - x1^2) * (F2 - F1)
   const double F1 = (x1 / gsl_pow_2(alpha)) * ( (x1 / (4.0 * alpha)) * (1.0/alpha - 1.0) - 1.0 ) + log(x1);
   const double F2 = (x2 / gsl_pow_2(alpha)) * ( (x2 / (4.0 * alpha)) * (1.0/alpha - 1.0) - 1.0 ) + log(x2);
   return 2.0 * (F2 - F1) / (gsl_pow_2(x2) - gsl_pow_2(x1));
 }


double   AT_RDD_Katz_PowerLawER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double r_max_m,
    const double alpha,
    const double Katz_point_coeff_Gy){

  //Dav(r1,r2) = coeff * kernel_av( x1, x2 )
  const double x1 = r1_m / r_max_m;
  const double x2 = r2_m / r_max_m;
  return Katz_point_coeff_Gy * AT_RDD_Katz_PowerLawER_DaverageKernel_approx(x1,x2,alpha);
}


double   AT_RDD_Katz_PowerLawER_dEdx_J_m(  const double a0_m,
    const double r_max_m,
    const double material_density_kg_m3,
    const double alpha,
    const double Katz_point_coeff_Gy){

  // dEdx = rho \int_a0^rmax r D(r) dr =
  //      = rho * (pi rmax^2 - pi a0^2) D_av(a0,rmax)
  return material_density_kg_m3 * M_PI * \
  (gsl_pow_2(r_max_m) - gsl_pow_2(a0_m)) * \
  AT_RDD_Katz_PowerLawER_Daverage_Gy(a0_m, r_max_m, r_max_m, alpha, Katz_point_coeff_Gy);
}


double   AT_RDD_Katz_PowerLawER_DSite_Gy( const double r_m,
    const double a0_m,
    const double r_max_m,
    const double material_density_kg_m3,
    const double alpha,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy){
  //Dsite(r) = 1 / (rho pi a0^2) * (LET - dEdx)   for r < a0
  //Dsite(r) = D(r)                               for r >= a0
  if( r_m < a0_m ){
    return M_1_PI * (LET_J_m - dEdx_J_m)/ (material_density_kg_m3 * gsl_pow_2(a0_m));
  } else {
    return AT_RDD_Katz_PowerLawER_Dpoint_Gy(r_m, alpha, r_max_m, Katz_point_coeff_Gy);
  }
}

//////////////////////////////////////// REST /////////////////////////////////////////////

inline float AT_RDD_Katz_dEdx_coeff_J_m(const float r_max_m,
    const float density_kg_m3,
    const float Katz_point_coeff_Gy){

  return 2 * M_PI * density_kg_m3 * r_max_m * r_max_m * Katz_point_coeff_Gy;
}


float AT_RDD_Katz_PowerLawER_dEdx_directVersion_J_m(  const float alpha,
    const float r_min_m,
    const float r_max_m,
    const float Katz_dEdx_coeff_J_m){

  double dEdx_integral = 0.0;
  if( r_min_m < r_max_m){
    dEdx_integral = (alpha/(1.+alpha)) * pow( 1.-r_min_m/r_max_m , 1. +
        1./alpha ) *
        gsl_sf_hyperg_2F1(1.,1.+1./alpha,2.+1./alpha,1.-r_min_m/r_max_m);
  }
  return Katz_dEdx_coeff_J_m*(float)dEdx_integral;
}

//////////////////////////////////////// Cucinotta model calculations /////////////////////////////////////////////

inline double   AT_RDD_Cucinotta_f_shortRange( const double r_m,
    const double beta){

  // fS(r) = 1.0/( r0/r + 0.6 + 1.7 beta + 1.1 beta^2)           [here r0 = 10^(-9) [m]]
  return 1.0/(1e-9/r_m + 0.6 + 1.7 * beta + 1.1 * gsl_pow_2(beta));
}


inline double   AT_RDD_Cucinotta_f_longRange( const double r_m,
    const double r_max_m){

  //fL(r) = exp( -(r/(0.37rmax))^2 )
  return exp( - gsl_pow_2( r_m / (0.37 * r_max_m ) ) );
}


inline double AT_RDD_Cucinotta_Ddelta_Gy( const double r_m,
    const double r_max_m,
    const double beta,
    const double Katz_point_coeff_Gy){

  // Ddelta(r) = coeff * fS(r) * fL(r) * rmax^2/r^2
  return Katz_point_coeff_Gy * AT_RDD_Cucinotta_f_longRange(r_m,r_max_m) * AT_RDD_Cucinotta_f_shortRange(r_m,beta) * gsl_pow_2(r_max_m/r_m);
}


double AT_RDD_Cucinotta_Ddelta_average_integrand_m(  double r_m,
    void * params){

  double r_max_m  =  ((double*)params)[0];
  double beta     =  ((double*)params)[1];
  double res      =  1.0/r_m;
  res            *=  AT_RDD_Cucinotta_f_shortRange( r_m, beta);
  res            *=  AT_RDD_Cucinotta_f_longRange( r_m, r_max_m);
  return res;
}


double   AT_RDD_Cucinotta_Ddelta_average_Gy(  const double r1_m,
    const double r2_m,
    const double r_max_m,
    const double beta,
    const double Katz_point_coeff_Gy){

  // fL(r) = exp( -(r/(0.37rmax))^2 )
  // fS(r) = 1.0/( r0/r + 0.6 + 1.7 beta + 1.1 beta^2)           [here r0 = 10^(-9) [m]]
  // Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 Ddelta(r) 2 pi r dr
  //  thus: ...
  // Dav(r1,r2) = 2 * coeff/ ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 fS(r) * fL(r) * 1/r dr


  if( (r2_m > r_max_m) || (r1_m > r_max_m) || (r1_m > r2_m)){
    printf("wrong parameters given to AT_RDD_Cucinotta_Ddelta_average_Gy\n");
    return 0.0;
  }

  gsl_set_error_handler_off(); // TODO improve error handling

  double delta_average_integral;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (10000);
  gsl_function F;
  F.function = &AT_RDD_Cucinotta_Ddelta_average_integrand_m;
  double params[] = {r_max_m, beta};
  F.params = params;
  int status = gsl_integration_qags (&F, r1_m, r2_m, 1e-11, 1e-7, 10000, w1, &delta_average_integral, &error);
  //printf("integral = %g , error = %g, status = %d\n", delta_average_integral, error, status);
  if (status > 0){
    printf("integration error %d in AT_RDD_Cucinotta_Ddelta_average_Gy\n", status);
    delta_average_integral = -1.0f;
  }
  gsl_integration_workspace_free (w1);

  return (2.0 * Katz_point_coeff_Gy / (gsl_pow_2(r2_m/r_max_m) - gsl_pow_2(r1_m/r_max_m))) * delta_average_integral ;
}


double   AT_RDD_Cucinotta_Dexc_average_Gy(  const double r1_m,
    const double r2_m,
    const double r_max_m,
    const double beta,
    const double Katz_point_coeff_Gy){

  if( (r2_m > r_max_m) || (r1_m > r_max_m) || (r1_m > r2_m) || (r1_m <= 0.0f) ){
    printf("wrong parameters given to AT_RDD_Cucinotta_Dexc_average_Gy\n");
    return 0.0f;
  }

  //  Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 Dexc(r) 2 pi r dr
  //  Dav(r1,r2) = 2 coeff/ ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 exp( - r / 2d ) * 1/r dr
  //
  // here:
  // integral = \int_r1^r2 exp( - r / 2d ) * 1/r * dr = Ei( -r2 / 2d ) - Ei( -r1 / 2d )
  // Ei is the exponential integral function
  const double wr = 13.0 * GSL_CONST_MKSA_ELECTRON_VOLT;
  const double d_m = 0.5 * beta * GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR * GSL_CONST_MKSA_SPEED_OF_LIGHT / wr; // TODO maybe d_m could be passed as argument
  double FA = gsl_sf_expint_Ei( -0.5 * r1_m / d_m );
  double FB = 0;
  if( r2_m / d_m < 100.0f ){ // if x < -50 then Ei(x) < 1e-24, so we can assume than FB = 0
    FB = gsl_sf_expint_Ei( -0.5 * r2_m / d_m );
  }
  return (2.0 * Katz_point_coeff_Gy / (gsl_pow_2(r2_m/r_max_m) - gsl_pow_2(r1_m/r_max_m))) * (FB - FA) ;
}


double   AT_RDD_Cucinotta_Cnorm( const double r_min_m,
    const double r_max_m,
    const double beta,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double Katz_point_coeff_Gy){

  const double LETfactor_Gy = LET_J_m * M_1_PI / (material_density_kg_m3 * (gsl_pow_2(r_max_m) - gsl_pow_2(r_min_m)));
  const double Ddelta_average_Gy = AT_RDD_Cucinotta_Ddelta_average_Gy(r_min_m,r_max_m,r_max_m,beta,Katz_point_coeff_Gy);
  const double Dexc_average_Gy = AT_RDD_Cucinotta_Dexc_average_Gy(r_min_m,r_max_m,r_max_m,beta,Katz_point_coeff_Gy);
  if( (LETfactor_Gy > 0.0) && (Ddelta_average_Gy > 0.0) && (Dexc_average_Gy > 0.0)){
    return (LETfactor_Gy - Ddelta_average_Gy) / Dexc_average_Gy;
  } else {
    printf("problem in AT_RDD_Cucinotta_Cnorm\n");
    return 0.0;
  }
}


double   AT_RDD_Cucinotta_Dexc_Gy( const double r_m,
    const double r_max_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy){

  // Dexc(r) = C exp( - r / 2d ) / r^2            [where d = (beta/2) * (hbar * c / wr) and wr = 13eV ]
  // Dexc(r) = Cnorm * coeff * exp( - r / 2d ) * (rmax/r)^2
  const double wr = 13.0 * GSL_CONST_MKSA_ELECTRON_VOLT;
  const double d_m = 0.5 * beta * GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR * GSL_CONST_MKSA_SPEED_OF_LIGHT / wr;
  return C_norm * Katz_point_coeff_Gy * exp( - 0.5 * r_m / d_m ) * gsl_pow_2(r_max_m/r_m);
}


inline double   AT_RDD_Cucinotta_Dpoint_Gy( const double r_m,
    const double r_max_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy){

  // D(r)    = Dexc(r) + Ddelta(r)
  return AT_RDD_Cucinotta_Dexc_Gy(r_m, r_max_m, beta, C_norm, Katz_point_coeff_Gy) + AT_RDD_Cucinotta_Ddelta_Gy(r_m, r_max_m, beta, Katz_point_coeff_Gy);
}


void AT_RDD_f1_parameters(  /* radiation field parameters */
    const float*  E_MeV_u,
    const long*  particle_no,
    /* detector parameters */
    const long*  material_no,
    /* radial dose distribution model */
    const long*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const long*  er_model,
    const float*  er_parameter,
    /* calculated parameters */
    float * f1_parameters)
{
  float LET_MeV_cm2_g = 0.0f;
  float max_electron_range_m = 0.0f;
  float r_min_m = 0.0f;
  float single_impact_fluence_cm2 = 0.0f;
  float single_impact_dose_Gy = 0.0f;
  float norm_constant_Gy = 0.0f;
  float d_min_Gy = 0.0f;
  float d_max_Gy = 0.0f;
  float dEdx_MeV_cm2_g = 0.0f;


  /*******************************************************************************
   *********************** GET MODEL INDEPENDENT VARIABLES ***********************
   *******************************************************************************/

  double beta   =  AT_beta_from_E_single( (double)(*E_MeV_u) );


  // Get density
  double   density_g_cm3, electron_density_m3;
  AT_get_material_data(
              *material_no,
              &density_g_cm3,
              &electron_density_m3,
              NULL, NULL, NULL, NULL, NULL, NULL);
  double density_kg_m3      =  density_g_cm3 * 1000.0;

  ///////////////////////////////////////////////////////////////////////////////
  // PARAMETER 0: Get the LET (same for all models)
  long  n_tmp    = 1;
  AT_LET_MeV_cm2_g(  n_tmp,
            E_MeV_u,
            particle_no,
            *material_no,
            &LET_MeV_cm2_g);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 2: Get the maximum electron range (same for all RDD models)
  max_electron_range_m = (float)AT_max_electron_range_m( (double)(*E_MeV_u), (int)(*material_no), (int)(*er_model));

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 6: Get the single impact fluence (same for all RDD models)
  AT_single_impact_fluence_cm2( &n_tmp,
      E_MeV_u,
      material_no,
      er_model,
      &single_impact_fluence_cm2);


  // get alpha for some ER models
  double alpha               =  0.0;
  if( (*er_model == ER_Waligorski) || (*er_model == ER_Edmund) ){
    double wmax_MeV     =  AT_max_E_transfer_MeV_single((double)(*E_MeV_u));
    alpha                   =  1.667;
    if(wmax_MeV <= 1e-3){  // if wmax < 1keV
      alpha                 =  1.079;
    } // end if
  }

  double C_J_m              =  0.0;
  double Katz_point_coeff_Gy  =  0.0;
  if( (*rdd_model == RDD_KatzPoint) || (*rdd_model == RDD_Site) || (*rdd_model == RDD_Edmund) || (*rdd_model == RDD_Cucinotta)){
    // TODO shall we move from J_m and m to more reasonable units ?
    // we have for water C_J_m = 1.22e-12 and r_m usually ~ 1e-8
    // calculations in C_J_um and r_um would be more precise
    const long Z            =  AT_Z_from_particle_no_single(*particle_no);
    const double Z_eff      =  AT_effective_charge_from_beta_single(beta, Z);
    C_J_m                   =  AT_RDD_Katz_C_J_m(electron_density_m3);
    Katz_point_coeff_Gy     =  AT_RDD_Katz_coeff_Gy(C_J_m, Z_eff, beta, density_kg_m3, (double)max_electron_range_m);
  }


  /*******************************************************************************
   *********************** GET MODEL SPECIFIC VARIABLES **************************
   *******************************************************************************
   ********* For every model following items should be calculated:     ***********
   *********  1. r_min_m or whatever that goes to (f1_parameters[1])    **********
   ********** 2. f1_parameters[3] <> d_min_Gy                           **********
   ********** 3. f1_parameters[4] <> d_max_Gy                           **********
   ********** 4. f1_parameters[5] <> norm_constant_Gy (or whatever)     **********
   ********** 5. f1_parameters[7] <> single_impact_dose_Gy;             **********
   ********** 6. f1_parameters[8] <> dEdx_MeV_cm2_g;                    **********
   *******************************************************************************
   *******************************************************************************/

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Test
  if( *rdd_model == RDD_Test){
    // 1. set minimum r_min_m (f1_parameters[1])
    r_min_m = 0.0f;                                                                             // r_min_m

    // 2. calculate minimum dose d_min_Gy (f1_parameters[3])
    // 3. calculate maximum dose d_max_Gy (f1_parameters[4])
    // 4. set norm_constant_Gy (f1_parameters[5])
    norm_constant_Gy = LET_MeV_cm2_g * single_impact_fluence_cm2 * MeV_to_J * 1000.0f;          // LET  / track area = Norm.constant k //TODO check units
    d_min_Gy = norm_constant_Gy;                                                                // d_min_Gy = k
    d_max_Gy = norm_constant_Gy;                                                                // d_max_Gy = k

    // 5. calculate single_impact_dose_Gy (f1_parameters[7])
    single_impact_dose_Gy = LET_MeV_cm2_g * single_impact_fluence_cm2 * MeV_g_to_J_kg;          // single_impact_dose = LET * fluence; //TODO check units
    // 6. calculate dEdx_MeV_cm2_g (f1_parameters[8])
    dEdx_MeV_cm2_g = LET_MeV_cm2_g;                                                             // dEdx = LET

    // Save parameters to f1_parameters table
    f1_parameters[1]  = r_min_m;
  } // end RDD_Test


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzPoint
  if( *rdd_model == RDD_KatzPoint){

    // 1. set minimum r_min_m (f1_parameters[1])
    r_min_m               =  rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m             =  max_electron_range_m;
    }                  // If r.max < r.min, r.min = r.max

    // 2. calculate minimum dose d_min_Gy (f1_parameters[3])
    d_min_Gy              =  rdd_parameter[1];                                                 // d_min_Gy given as model parameter

    // 3. calculate maximum dose d_max_Gy (f1_parameters[4])
    if( (*er_model == ER_Waligorski) || (*er_model == ER_Edmund) ){ // "new" Katz RDD
      d_max_Gy            =  (float)AT_RDD_Katz_PowerLawER_Dpoint_Gy((double)r_min_m, alpha, (double)max_electron_range_m, Katz_point_coeff_Gy); // d_max_Gy given as dose for r_min_m
    } else if (*er_model == ER_ButtsKatz){ // "old" Katz RDD
      d_max_Gy            =  (float)AT_RDD_Katz_LinearER_Dpoint_Gy((double)r_min_m, (double)max_electron_range_m, Katz_point_coeff_Gy); // d_max_Gy given as dose for r_min_m
    } else { // not supported ER model
      d_max_Gy            =  0.0f;
    }
    // 4. set norm_constant_Gy (f1_parameters[5])
    norm_constant_Gy      =  0.0f;                                                              // not used here as this RDD model is not normalized

    // 5. calculate single_impact_dose_Gy (f1_parameters[7])
    single_impact_dose_Gy =  LET_MeV_cm2_g * single_impact_fluence_cm2 * MeV_g_to_J_kg;        // single_impact_dose = LET * fluence; TODO check how it should be defined

    // 6. calculate dEdx_MeV_cm2_g (f1_parameters[8])
    dEdx_MeV_cm2_g        =  LET_MeV_cm2_g;                                                    // dEdx = LET ?? // TODO here we can calculate what is \int D(r) r dr

    // Save parameters to f1_parameters table
    f1_parameters[1]      =  r_min_m;
  } // end RDD_KatzPoint


  //TODO check if there is any difference between RDD_Site and RDD_Edmund
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Site
  if( *rdd_model == RDD_Site){

    // 1. set minimum r_min_m (f1_parameters[1])
    r_min_m               =  rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m             =  max_electron_range_m;
    }                  // If r.max < r.min, r.min = r.max

    // 2. calculate minimum dose d_min_Gy (f1_parameters[3])
    d_min_Gy              =  rdd_parameter[1];                                                  // d_min_Gy given as model parameter

    // 3. calculate maximum dose d_max_Gy (f1_parameters[4])
    double  LET_J_m       =  (double)LET_MeV_cm2_g * density_g_cm3; // [MeV / cm]
    LET_J_m              *=  100.0;       // [MeV / m]
    LET_J_m              *=  MeV_to_J;     // [J/m]
    float  dEdx_J_m       =  0.0f;

    if( (*er_model == ER_Waligorski) || (*er_model == ER_Edmund) ){ // calculate dEdx_MeV_cm2_g from "new" Katz RDD
      dEdx_J_m            =  (float)AT_RDD_Katz_PowerLawER_dEdx_J_m((double)r_min_m, (double)max_electron_range_m, density_kg_m3, alpha, Katz_point_coeff_Gy);
      d_max_Gy            =  (float)AT_RDD_Katz_PowerLawER_DSite_Gy(0.0, (double)r_min_m, (double)max_electron_range_m, density_kg_m3, alpha, LET_J_m, (double)dEdx_J_m, Katz_point_coeff_Gy);

    } else if (*er_model == ER_ButtsKatz){ // calculate dEdx_MeV_cm2_g from "old" Katz RDD
      dEdx_J_m            =  (float)AT_RDD_Katz_LinearER_dEdx_J_m((double)r_min_m, (double)max_electron_range_m, density_kg_m3, Katz_point_coeff_Gy);
      d_max_Gy            =  (float)AT_RDD_Katz_LinearER_DSite_Gy(0.0, (double)r_min_m, (double)max_electron_range_m, density_kg_m3, LET_J_m, (double)dEdx_J_m, Katz_point_coeff_Gy);
    } else {
      d_max_Gy            =  0.0f;
    }

    // 4. set norm_constant_Gy (f1_parameters[5])
    norm_constant_Gy      =  0.0f;

    // 5. calculate single_impact_dose_Gy (f1_parameters[7])
    single_impact_dose_Gy =  LET_MeV_cm2_g * single_impact_fluence_cm2 * MeV_g_to_J_kg;   // single_impact_dose = LET * fluence; TODO check how it should be defined

    // 6. calculate dEdx_MeV_cm2_g (f1_parameters[8])
    dEdx_MeV_cm2_g        =  dEdx_J_m / 100.0f / density_g_cm3 / MeV_to_J;

    // Save parameters to f1_parameters table
    f1_parameters[1]      =  r_min_m;
  } // end RDD_Site

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Edmund
  //TODO rewrite in better way
  if( *rdd_model == RDD_Edmund){

    // 1. set minimum r_min_m (f1_parameters[1])
    r_min_m               =  rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m = max_electron_range_m;
    }                  // If r.max < r.min, r.min = r.max

    // 2. calculate minimum dose d_min_Gy (f1_parameters[3])
    d_min_Gy              =  rdd_parameter[1];

    // 3. calculate maximum dose d_max_Gy (f1_parameters[4])
    if( (*er_model == ER_Waligorski) || (*er_model == ER_Edmund) ){ // This model will work only with ER_Waligorski and ER_Edmund
      // Get dEdx by simple integration from r_min_m to r_max_m
      float  dEdx_J_m     =  0.0f;

      float Katz_dEdx_coeff_J_m = AT_RDD_Katz_dEdx_coeff_J_m(max_electron_range_m, density_kg_m3, (float)Katz_point_coeff_Gy);
      dEdx_J_m            =  AT_RDD_Katz_PowerLawER_dEdx_directVersion_J_m((float)alpha, r_min_m, max_electron_range_m, Katz_dEdx_coeff_J_m);

      dEdx_MeV_cm2_g      =  dEdx_J_m / 100.0f / density_g_cm3 / MeV_to_J;

      float tmp           =  LET_MeV_cm2_g - dEdx_MeV_cm2_g;  // LET - dEdx (MeV_g_cm2)
      if(tmp > 0){
        tmp *=  density_g_cm3; // MeV / cm
        tmp *=  MeV_to_J;      // J / cm
        tmp *=  100.0f;        // J / m
        float core_kg_m   =  M_PI * gsl_pow_2(rdd_parameter[0]) * (float)density_kg_m3;  // kg / m
        d_max_Gy          =  tmp / core_kg_m;                                     // J / kg = Gy
      }else{
        d_max_Gy = 1e-11;
      } // end if tmp

    } else { // er_models other than ER_Waligorski or ER_Edmund
      dEdx_MeV_cm2_g      =  0.0f;
      d_max_Gy            =  0.0f;
    } //end if er_model

    // 4. set norm_constant_Gy (f1_parameters[5])
    norm_constant_Gy      =  (float)C_J_m;                      // Norm.constant k // TODO is Gy = J_m ?

    // 5. calculate single_impact_dose_Gy (f1_parameters[7])
    single_impact_dose_Gy =  LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;        // LET * fluence

    // 6. done after 3. and before 4.

    // Save parameters to f1_parameters table
    f1_parameters[1]      =  r_min_m;
  }// end RDD_Edmund


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Geiss
  if( *rdd_model == RDD_Geiss){

    // 1. set a0_m (f1_parameters[1])
    float a0_m            =  0.0;
    a0_m                  =  rdd_parameter[0];                                            // "r_min_m" = a0
    if (max_electron_range_m <= a0_m){
      a0_m                =  max_electron_range_m;
    }      // If r.max < a0, r.min = r.max, not a0

    // 2. calculate minimum dose d_min_Gy (f1_parameters[3])
    // 3. calculate maximum dose d_max_Gy (f1_parameters[4])
    // 4. set norm_constant_Gy (f1_parameters[5])
    float  tmp_cm         =  (float)(0.5f + log(max_electron_range_m / a0_m));
    tmp_cm               *=  2.0f * M_PI * gsl_pow_2(a0_m * m_to_cm);                     //Normalization to match with LET
    norm_constant_Gy      =  LET_MeV_cm2_g * MeV_g_to_J_kg / tmp_cm;                      // k = LET / tmp
    d_max_Gy              =  norm_constant_Gy;                                            // d_max_Gy = k
    d_min_Gy              =  norm_constant_Gy * gsl_pow_2( a0_m / max_electron_range_m);  // d_min_Gy

    // 5. calculate single_impact_dose_Gy (f1_parameters[7])
    single_impact_dose_Gy =  LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;   // single_impact_dose = LET * single_impact_fluence

    // 6. calculate dEdx_MeV_cm2_g (f1_parameters[8])
    dEdx_MeV_cm2_g        =  LET_MeV_cm2_g;                                               // dEdx = LET

    // Save parameters to f1_parameters table
    f1_parameters[1]      =  a0_m;
  } // end RDD_Geiss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Cucinotta
  if( *rdd_model == RDD_Cucinotta){ // TODO RDD does Cucinotta RDD model work only with Tabata ER ?

    // 1. set minimum r_min_m (f1_parameters[1])
    r_min_m               =  rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m             =  max_electron_range_m;
    }                  // If r.max < r.min, r.min = r.max

    // 2. calculate minimum dose d_min_Gy (f1_parameters[3])
    d_min_Gy              =  rdd_parameter[1];

    // 3. calculate maximum dose d_max_Gy (f1_parameters[4])
    // 4. set norm_constant_Gy (f1_parameters[5])
    double  LET_J_m       =  (double)LET_MeV_cm2_g * density_g_cm3; // [MeV / cm]
    LET_J_m              *=  100.0;        // [MeV / m]
    LET_J_m              *=  MeV_to_J;     // [J/m]

    norm_constant_Gy      =  (float)AT_RDD_Cucinotta_Cnorm((double)r_min_m, (double)max_electron_range_m, beta, density_kg_m3, LET_J_m, Katz_point_coeff_Gy);
    if( norm_constant_Gy == 0){
      printf("problem\n"); // TODO handle this situation
    }

    d_max_Gy              =  (float)AT_RDD_Cucinotta_Dpoint_Gy((double)r_min_m, (double)max_electron_range_m, beta, (double)norm_constant_Gy, Katz_point_coeff_Gy);

    // 5. calculate single_impact_dose_Gy (f1_parameters[7])
    single_impact_dose_Gy =  LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;        // LET * fluence

    // 6. calculate dEdx_MeV_cm2_g (f1_parameters[8])
    dEdx_MeV_cm2_g        =  LET_MeV_cm2_g;   // TODO move norm_constant_Gy to dEdx_MeV_cm2_g

    // Save parameters to f1_parameters table
    f1_parameters[1]      =  r_min_m;
  }// end RDD_Cucinotta


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzExtTarget
  if( *rdd_model == RDD_KatzExtTarget){ // TODO To be removed

    // 1. set minimum r_min_m (f1_parameters[1])
    r_min_m               =  rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m             =  max_electron_range_m;
    }                  // If r.max < r.min, r.min = r.max

    // 2. calculate minimum dose d_min_Gy (f1_parameters[3])
    d_min_Gy              =  rdd_parameter[2];

    // 3. calculate maximum dose d_max_Gy (f1_parameters[4])
    const long  n_tmp     =  1;
    const float a0_m      =  rdd_parameter[1];
    const float rdd_basic_parameter[] = {rdd_parameter[0], rdd_parameter[2]};

    AT_RDD_ExtendedTarget_Gy(n_tmp, &r_min_m, a0_m, *E_MeV_u, *particle_no, *material_no, RDD_KatzPoint, rdd_basic_parameter, *er_model, er_parameter, &d_max_Gy);

    // 4. set norm_constant_Gy (f1_parameters[5])
    norm_constant_Gy      =  0.0f;

    // 5. calculate single_impact_dose_Gy (f1_parameters[7])
    single_impact_dose_Gy =  LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;        // LET * fluence

    // 6. calculate dEdx_MeV_cm2_g (f1_parameters[8])
    dEdx_MeV_cm2_g        =  LET_MeV_cm2_g;   // TODO move norm_constant_Gy to dEdx_MeV_cm2_g

    // Save parameters to f1_parameters table
    f1_parameters[1]      =  r_min_m;
  }// end RDD_KatzExtTarget

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaExtTarget
  if( *rdd_model == RDD_CucinottaExtTarget){ // TODO To be removed

    // 1. set minimum r_min_m (f1_parameters[1])
    r_min_m               =  rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m             =  max_electron_range_m;
    }                  // If r.max < r.min, r.min = r.max

    // 2. calculate minimum dose d_min_Gy (f1_parameters[3])
    d_min_Gy              =  rdd_parameter[2];

    // 3. calculate maximum dose d_max_Gy (f1_parameters[4])
    const long  n_tmp     =  1;
    const float a0_m      =  rdd_parameter[1];
    const float rdd_basic_parameter[] = {rdd_parameter[0], rdd_parameter[2]};

    AT_RDD_ExtendedTarget_Gy(n_tmp, &r_min_m, a0_m, *E_MeV_u, *particle_no, *material_no, RDD_Cucinotta, rdd_basic_parameter, *er_model, er_parameter, &d_max_Gy);

    // 4. set norm_constant_Gy (f1_parameters[5])
    norm_constant_Gy      =  0.0f;

    // 5. calculate single_impact_dose_Gy (f1_parameters[7])
    single_impact_dose_Gy =  LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;        // LET * fluence

    // 6. calculate dEdx_MeV_cm2_g (f1_parameters[8])
    dEdx_MeV_cm2_g        =  LET_MeV_cm2_g;   // TODO move norm_constant_Gy to dEdx_MeV_cm2_g

    // Save parameters to f1_parameters table
    f1_parameters[1]      =  r_min_m;
  }// end RDD_CucinottaExtTarget


  // write data to output table (apart from f1_parameters[0] which sometimes is
  // r_min_m and sometimes a0_m )
  f1_parameters[0]  = LET_MeV_cm2_g;
  f1_parameters[2]  = max_electron_range_m;
  f1_parameters[3]  = d_min_Gy;
  f1_parameters[4]  = d_max_Gy;
  f1_parameters[5]  = norm_constant_Gy;
  f1_parameters[6]  = single_impact_fluence_cm2;
  f1_parameters[7]  = single_impact_dose_Gy;
  f1_parameters[8]  = dEdx_MeV_cm2_g;

}

void AT_D_RDD_Gy  ( const  long*  n,
    const float*  r_m,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const long*  particle_no,
    /* detector parameters */
    const long*  material_no,
    /* radial dose distribution model */
    const long*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const long*  er_model,
    const float*  er_parameter,
    float*  D_RDD_Gy)
{
  /********************************************************
   ********* CALCULATION BEFORE PARTICLE LOOP *************
   *******************************************************/

  // Get f1 parameters
  const long n_f1_parameters = 9;
  float*   f1_parameters     = (float*)calloc(n_f1_parameters, sizeof(float));

  AT_RDD_f1_parameters(      /* radiation field parameters */
      E_MeV_u,
      particle_no,
      /* detector parameters */
      material_no,
      /* radial dose distribution model */
      rdd_model,
      rdd_parameter,
      /* electron range model */
      er_model,
      er_parameter,
      /* calculated parameters */
      f1_parameters);


  // f1_parameters decoded
  const float LET_MeV_cm2_g             =  f1_parameters[0];
  const float max_electron_range_m      =  f1_parameters[2];
  const float d_min_Gy                  =  f1_parameters[3];
  const float d_max_Gy                  =  f1_parameters[4];
  const float norm_constant_Gy          =  f1_parameters[5];
  const float dEdx_MeV_cm2_g            =  f1_parameters[8];

  // Get material data
  double   density_g_cm3, electron_density_m3;
  AT_get_material_data(
              *material_no,
              &density_g_cm3,
              &electron_density_m3,
              NULL, NULL, NULL, NULL, NULL, NULL);
  double density_kg_m3      =  density_g_cm3 * 1000.0;

  double beta    =  AT_beta_from_E_single( (double)(*E_MeV_u) );

  // get alpha for some ER models
  double alpha               =  0.0;
  if( (*er_model == ER_Waligorski) || (*er_model == ER_Edmund) ){
    double wmax_MeV     =  AT_max_E_transfer_MeV_single((double)(*E_MeV_u));
    alpha                   =  1.667;
    if(wmax_MeV <= 1e-3){  // if wmax < 1keV
      alpha                 =  1.079;
    } // end if
  }

  double Katz_point_coeff_Gy  =  0.0;
  if( (*rdd_model == RDD_KatzPoint) || (*rdd_model == RDD_Site) || (*rdd_model == RDD_Edmund) || (*rdd_model == RDD_Cucinotta)){
    // TODO shall we move from J_m and m to more reasonable units ?
    // we have for water C_J_m = 1.22e-12 and r_m usually ~ 1e-8
    // calculations in C_J_um and r_um would be more precise
    const long  Z           =  AT_Z_from_particle_no_single(*particle_no);
    const double Z_eff      =  AT_effective_charge_from_beta_single(beta, Z);
    const double C_J_m      =  AT_RDD_Katz_C_J_m(electron_density_m3);
    Katz_point_coeff_Gy     =  AT_RDD_Katz_coeff_Gy(C_J_m, Z_eff, beta, density_kg_m3, (double)max_electron_range_m);
  }


  /********************************************************
   *************** LOOP OVER DISTANCE VECTOR **************
   *******************************************************/
  long     i;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Test
  if( *rdd_model == RDD_Test){
    // Loop over all r_m given
    for (i = 0; i < *n; i++){
      D_RDD_Gy[i]         =  0.0f;
      if (r_m[i] < max_electron_range_m){
        D_RDD_Gy[i]       =  norm_constant_Gy;
      } else {
        D_RDD_Gy[i]       =  0;
      }
    }
  }// end RDD_Test

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzPoint
  if( *rdd_model == RDD_KatzPoint){ // RDD formula will be determined by form of ER model
    const float r_min_m  =  f1_parameters[1];


    // use power law form of RDD for power law ER models (new Katz style)
    if( (*er_model == ER_Waligorski) || (*er_model == ER_Edmund) ){

      // Loop over all r_m given
      for (i = 0; i < *n; i++){
        D_RDD_Gy[i]       =  0.0f;
        if (r_m[i] >= r_min_m && r_m[i] <= max_electron_range_m){
          D_RDD_Gy[i]     =  (float)AT_RDD_Katz_PowerLawER_Dpoint_Gy((double)r_m[i], alpha, (double)max_electron_range_m, Katz_point_coeff_Gy);
          // D_RDD_Gy[i]        = fmaxf(D_RDD_Gy[i], d_min_Gy);          // Cut-off low doses // TODO is this cutoff necessary here ?
        } // end if
      } // end for

    }  else if (*er_model == ER_ButtsKatz){ // linear law form of RDD for other ER models (old Katz style)

      // Loop over all r_m given
      for (i = 0; i < *n; i++){
        D_RDD_Gy[i]       =  0.0f;
        if (r_m[i] >= r_min_m && r_m[i] <= max_electron_range_m){
          D_RDD_Gy[i]     = (float)AT_RDD_Katz_LinearER_Dpoint_Gy((double)r_m[i], (double)max_electron_range_m, Katz_point_coeff_Gy);
          // D_RDD_Gy[i]        = fmaxf(D_RDD_Gy[i], d_min_Gy);          // Cut-off low doses // TODO is this cutoff necessary here ?
        } // end if
      } // end for

    } else { // other not supported ER models

      // Loop over all r_m given
      for (i = 0; i < *n; i++){
        D_RDD_Gy[i]       =  0.0f;
      } // end for
    } // end if

  }// end RDD_KatzPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Site
  if( *rdd_model == RDD_Site){
    const double r_min_m   =  (double)f1_parameters[1];

    // convert LET_MeV_cm2_g to LET_J_m
    double  LET_J_m       =  (double)LET_MeV_cm2_g * density_g_cm3; // [MeV / cm]
    LET_J_m              *=  100.0;       // [MeV / m]
    LET_J_m              *=  MeV_to_J;     // [J/m]

    // convert dEdx_MeV_cm2_g to dEdx_J_m
    const double dEdx_J_m  = (double)dEdx_MeV_cm2_g * 100.0 * density_g_cm3 * MeV_to_J;

    if( (*er_model == ER_Waligorski) || (*er_model == ER_Edmund) ){

      // Loop over all r_m given
      for (i = 0; i < *n; i++){
        D_RDD_Gy[i]       =  0.0f;
        if (r_m[i] >= 0 && r_m[i] <= max_electron_range_m){
          D_RDD_Gy[i]     =  (float)AT_RDD_Katz_PowerLawER_DSite_Gy((double)r_m[i], r_min_m, (double)max_electron_range_m, density_kg_m3, alpha, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
        } // end if
      } // end for
    } else if (*er_model == ER_ButtsKatz){

      // Loop over all r_m given
      for (i = 0; i < *n; i++){
        D_RDD_Gy[i]       =  0.0f;
        if (r_m[i] >= 0 && r_m[i] <= max_electron_range_m){
          D_RDD_Gy[i]     =  (float)AT_RDD_Katz_LinearER_DSite_Gy((double)r_m[i], r_min_m, (double)max_electron_range_m, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
        } // end if
      } // end for
    } else { // not supported ER models

      // Loop over all r_m given
      for (i = 0; i < *n; i++){
        D_RDD_Gy[i]       =  0.0f;
      } // end for
    } //end if

  }// end RDD_Site

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Edmund
  if(*rdd_model == RDD_Edmund){
    const float r_min_m   =  f1_parameters[1];

    // This model will work only with ER_Waligorski and ER_Edmund
    if( (*er_model == ER_Waligorski) || (*er_model == ER_Edmund) ){
      // Loop over all r_m given
      for (i = 0; i < *n; i++){
        D_RDD_Gy[i]       =  0.0f;                                  // r < r_min_m (for RDD_KatzPoint) or r > r_max_m --> D = 0
        if (r_m[i] >= r_min_m && r_m[i] <= max_electron_range_m){          // in between r_min and r_max --> D = KatzPoint
          D_RDD_Gy[i]     =  (float)AT_RDD_Katz_PowerLawER_Dpoint_Gy((double)r_m[i], alpha, (double)max_electron_range_m, Katz_point_coeff_Gy);
          D_RDD_Gy[i]     =  fmaxf(D_RDD_Gy[i], d_min_Gy);          // Cut-off low doses
        } // end if
        if (r_m[i] < r_min_m){            // r < r_min_m --> D = d_max_Gy
          D_RDD_Gy[i]     =  d_max_Gy;
        } // end if
      } // end for
    } else {
      // Loop over all r_m given
      for (i = 0; i < *n; i++){
        D_RDD_Gy[i]       =  0.0f;
      } // end for
    } // end if

  }// end RDD_Site

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Geiss
  if( *rdd_model == RDD_Geiss){

    const float a0_m      =  f1_parameters[1];

    // Loop over all r_m given
    for (i = 0; i < *n; i++){
      D_RDD_Gy[i]         =  0.0f;
      if (r_m[i] < a0_m){
        D_RDD_Gy[i]       =  norm_constant_Gy;
      }
      if ((a0_m <= r_m[i]) && (r_m[i] <= max_electron_range_m)){
        D_RDD_Gy[i]       =  norm_constant_Gy * gsl_pow_2(a0_m / r_m[i]);
      }
    }
  }// end RDD_Geiss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Cucinotta
  if( *rdd_model == RDD_Cucinotta){

    const float r_min_m   =  f1_parameters[1];
    const double C_norm   =  (double)norm_constant_Gy;

    // Loop over all r_m given
    for (i = 0; i < *n; i++){
      D_RDD_Gy[i]         =  0.0f;
      if (r_m[i] < r_min_m){
        D_RDD_Gy[i]       =  0.0f;
      }
      if ((r_min_m <= r_m[i]) && (r_m[i] <= max_electron_range_m)){
        D_RDD_Gy[i]       =  (float)AT_RDD_Cucinotta_Dpoint_Gy((double)r_m[i], (double)max_electron_range_m, beta, C_norm, Katz_point_coeff_Gy);
      }
    }
  }// end RDD_Cucinotta

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzExtTarget TODO remove this dirty fix
  if( *rdd_model == RDD_KatzExtTarget){

    const float a0_m   =  rdd_parameter[1];

    const float rdd_basic_parameter[] = {rdd_parameter[0], rdd_parameter[2]};

    AT_RDD_ExtendedTarget_Gy(*n, r_m, a0_m, *E_MeV_u, *particle_no, *material_no, RDD_KatzPoint, rdd_basic_parameter, *er_model, er_parameter, D_RDD_Gy);
  }// end RDD_KatzExtTarget

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaExtTarget TODO remove this dirty fix
  if( *rdd_model == RDD_CucinottaExtTarget){

    const float a0_m   =  rdd_parameter[1];

    const float rdd_basic_parameter[] = {rdd_parameter[0], rdd_parameter[2]};

    AT_RDD_ExtendedTarget_Gy(*n, r_m, a0_m, *E_MeV_u, *particle_no, *material_no, RDD_Cucinotta, rdd_basic_parameter, *er_model, er_parameter, D_RDD_Gy);
  } // end RDD_CucinottaExtTarget


  free(f1_parameters);
}

float AT_D_RDD_Gy_solver( const float r , void * params ){
  AT_D_RDD_Gy_parameters* params_struct = (AT_D_RDD_Gy_parameters*)(params);
  *((*params_struct).n) = 1;
  params_struct->r_m = (float*)calloc(1,sizeof(float));
  (params_struct->r_m)[0] = r;
  AT_D_RDD_Gy(  params_struct->n,
          params_struct->r_m,
          params_struct->E_MeV_u,
          params_struct->particle_no,
          params_struct->material_no,
          params_struct->rdd_model,
          params_struct->rdd_parameter,
          params_struct->er_model,
          params_struct->er_parameter,
          params_struct->D_RDD_Gy);
  free(params_struct->r_m);
  return (params_struct->D_RDD_Gy)[0]-(params_struct->D0);
}

// TODO needs to be rewritten
void AT_r_RDD_m  ( const  long*  n,
    const float*  D_RDD_Gy,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const long*  particle_no,
    /* detector parameters */
    const long*  material_no,
    /* radial dose distribution model */
    const long*  rdd_model,       /* */
    const float*  rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
    /* electron range model */
    const long*  er_model,
    const float*  er_parameter,
    float*  r_RDD_m)
{
  long     i;
  // Get f1 parameters
  long     n_f1_parameters = 9;
  float*   f1_parameters   = (float*)calloc(n_f1_parameters, sizeof(float));
  AT_RDD_f1_parameters(      /* radiation field parameters */
    E_MeV_u,
    particle_no,
    /* detector parameters */
    material_no,
    /* radial dose distribution model */
    rdd_model,
    rdd_parameter,
    /* electron range model */
    er_model,
    er_parameter,
    /* calculated parameters */
    f1_parameters);

  // f1_parameters decoded
  //const float LET_MeV_cm2_g             =  f1_parameters[0]; //TODO not used
  const float max_electron_range_m      =  f1_parameters[2];
  const float d_min_Gy                  =  f1_parameters[3];
  const float d_max_Gy                  =  f1_parameters[4];
  const float norm_constant_Gy          =  f1_parameters[5];
  //const float dEdx_MeV_cm2_g            =  f1_parameters[8]; //TODO not used

  // Get material data
  double     density_g_cm3  =  AT_density_g_cm3_from_material_no(*material_no);
  double     density_kg_m3  =  density_g_cm3 * 1000.0;

  if( *rdd_model == RDD_Test){
    // Loop over all doses given
    for (i = 0; i < *n; i++){
      r_RDD_m[i]    =  0.0f;
      if (D_RDD_Gy[i] > 0.0f){
        r_RDD_m[i]    =  max_electron_range_m;
      }
    }
  }// end RDD_Test

  if( *rdd_model == RDD_Geiss){
    // Loop over all doses given
    for (i = 0; i < *n; i++){
      // if D is outside the definition, return -1
      r_RDD_m[i]    =  -1.0f;
      if ((d_min_Gy <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= d_max_Gy)){
        r_RDD_m[i]    =  f1_parameters[1] * (float)sqrt(norm_constant_Gy / D_RDD_Gy[i]);}
    }
  }// end RDD_Geiss

  long     n_tmp      = 1;
  if( (*rdd_model == RDD_KatzPoint) || (*rdd_model == RDD_Site) || (*rdd_model == RDD_KatzExtTarget) || (*rdd_model == RDD_Edmund)){

    const float r_min_m  =  f1_parameters[1];

    AT_D_RDD_Gy_parameters* params = (AT_D_RDD_Gy_parameters*)calloc(1,sizeof(AT_D_RDD_Gy_parameters));
    params->n             = (long*)calloc(1,sizeof(long));
    params->E_MeV_u       = (float*)E_MeV_u;
    params->particle_no   = (long*)particle_no;
    params->material_no   = (long*)material_no;
    params->rdd_model     = (long*)rdd_model;
    params->rdd_parameter = (float*)rdd_parameter;
    params->er_model      = (long*)er_model;
    params->er_parameter  = (float*)er_parameter;
    params->D_RDD_Gy      = (float*)calloc(1,sizeof(float));
    // Loop over all doses given
    float dev = 0.01f;
    float critical_d_Gy    = 0.0f;
    float critical_r_m     = 0.0f;
    float inv2_d_Gy        = 0.0f;
    float inv2_r_m         = 0.0f;
    if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund){
      critical_r_m   = rdd_parameter[0] * (1.0f + 1e-6f);
      inv2_r_m       = fmaxf(rdd_parameter[0], max_electron_range_m * dev);
    }
    if(*rdd_model == RDD_KatzExtTarget){
      critical_r_m   = rdd_parameter[1];
      inv2_r_m       = fmaxf(rdd_parameter[0], max_electron_range_m * dev);
    }
    if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund || *rdd_model == RDD_KatzExtTarget){
      AT_D_RDD_Gy  (  &n_tmp,                    // Use D(r) to find dose at jump of D_Site
                &critical_r_m,
                /* radiation field parameters */
                E_MeV_u,
                particle_no,
                /* detector parameters */
                material_no,
                /* radial dose distribution model */
                rdd_model,
                rdd_parameter,
                /* electron range model */
                er_model,
                er_parameter,
                &critical_d_Gy);
      AT_D_RDD_Gy  (  &n_tmp,                    // Use D(r) to find dose at jump of D_Site
                &inv2_r_m,
                /* radiation field parameters */
                E_MeV_u,
                particle_no,
                /* detector parameters */
                material_no,
                /* radial dose distribution model */
                rdd_model,
                rdd_parameter,
                /* electron range model */
                er_model,
                er_parameter,
                &inv2_d_Gy);
    }

///////////////////////////////////////////////////////////////////////////////////
    // Get beta, Z and Zeff
    const double beta   =  AT_beta_from_E_single( (double)(*E_MeV_u) );
    const long   Z      =  AT_Z_from_particle_no_single(*particle_no);
    const double Z_eff  =  AT_effective_charge_from_beta_single(beta, Z);

    const double Katz_point_coeff_Gy = AT_RDD_Katz_coeff_Gy((double)norm_constant_Gy, Z_eff, beta, density_kg_m3, (double)max_electron_range_m);

    for (i = 0; i < *n; i++){
      params->D0 = D_RDD_Gy[i];
      // if D is outside the definition, return -1
      r_RDD_m[i]        =  -1.0f;
      float  solver_accuracy  =  1e-13f;
      if ((d_min_Gy <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= d_max_Gy)){
        if(D_RDD_Gy[i] >= critical_d_Gy){
          if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund){
            r_RDD_m[i] = rdd_parameter[0];}
          if(*rdd_model == RDD_KatzExtTarget){
            r_RDD_m[i] = rdd_parameter[1];}            
        }
        if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund || *rdd_model == RDD_KatzExtTarget){
          if(D_RDD_Gy[i] < critical_d_Gy && D_RDD_Gy[i] >= inv2_d_Gy){
            r_RDD_m[i] = sqrt((float)Katz_point_coeff_Gy / D_RDD_Gy[i]) * max_electron_range_m;}
          if(D_RDD_Gy[i] < inv2_d_Gy){
            r_RDD_m[i] = zriddr(AT_D_RDD_Gy_solver,
                    (void*)params,
                    r_min_m,
                    max_electron_range_m,
                    solver_accuracy);
          }
        }else{
          r_RDD_m[i] = zriddr(AT_D_RDD_Gy_solver,
                  (void*)params,
                  r_min_m,
                  max_electron_range_m,
                  solver_accuracy);

        }
      }
    }
    free(params->n);
    free(params->D_RDD_Gy);
    free(params);
  }// end RDD_KatzPoint & RDD_Site

  free(f1_parameters);

}


double AT_D_RDD_Gy_int( double  r_m,
    void*   params){
  float D_Gy;
  long  n_tmp              = 1;
  AT_P_RDD_parameters* par = (AT_P_RDD_parameters*)params;
  float fr_m               = (float)r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  &n_tmp,
      &fr_m,
      par->E_MeV_u,
      par->particle_no,
      par->material_no,
      par->rdd_model,
      par->rdd_parameters,
      par->er_model,
      par->er_parameters,
      &D_Gy);

  return (2.0 * M_PI * r_m * (double)D_Gy);
}

double AT_sI_int( double  r_m,
    void*   params){
  double  P = AT_P_RDD(r_m, params);
  return (r_m * P);
}

double AT_P_RDD( double  r_m,
    void*   params)
{
  float  D_Gy;
  long   n_tmp              = 1;
  AT_P_RDD_parameters* par  = (AT_P_RDD_parameters*)params;
  float  fr_m               = (float)r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  &n_tmp,
      &fr_m,
      par->E_MeV_u,
      par->particle_no,
      par->material_no,
      par->rdd_model,
      par->rdd_parameters,
      par->er_model,
      par->er_parameters,
      &D_Gy);

  long gamma_model = GR_GeneralTarget;
  float P;
  AT_gamma_response(  n_tmp,
      &D_Gy,
      gamma_model,
      par->gamma_parameters,
      // return
      &P);
  return ((double)P);
}
