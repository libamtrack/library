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


long getRDDNo(const char* RDD_name){
  long i;
  for (i = 0; i < RDD_DATA_N; i++){
    if (strcmp(RDD_name, AT_RDD_Data.RDD_name[i]) == 0){
      return AT_RDD_Data.RDD_no[i];
    }
  }
  return -1;
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
  if( r_m > r_max_m )
    return 0.0;
  else
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
  if( r_m <= a0_m ){
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

////////////////////////////////////// GENERAL METHODS //////////////////////////////////////


inline double  AT_RDD_Test_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double norm_constant_Gy){

  if( r_m >= r_min_m && r_m <= r_max_m ){
    return norm_constant_Gy;
  } else {
    return 0.0;
  }
}


inline double  AT_RDD_KatzPoint_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const long   er_model,
    const double alpha,
    const double Katz_point_coeff_Gy){

  if (r_m >= r_min_m && r_m <= r_max_m){
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // use power law form of RDD for power law ER models (new Katz style)
      return AT_RDD_Katz_PowerLawER_Dpoint_Gy(r_m, alpha, r_max_m, Katz_point_coeff_Gy);
    } else if (er_model == ER_ButtsKatz){ // linear law form of RDD for other ER models (old Katz style)
      return AT_RDD_Katz_LinearER_Dpoint_Gy(r_m, r_max_m, Katz_point_coeff_Gy);
    }
  }
  return 0.;
}


inline double  AT_RDD_KatzSite_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double a0_m,
    const long er_model,
    const double alpha,
    const double density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy){

  if (r_m >= r_min_m && r_m <= r_max_m){
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){
      return AT_RDD_Katz_PowerLawER_DSite_Gy(r_m, a0_m, r_max_m, density_kg_m3, alpha, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
    } else if (er_model == ER_ButtsKatz){
      return AT_RDD_Katz_LinearER_DSite_Gy(r_m, a0_m, r_max_m, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
    }
  }
  return 0.;
}


inline double  AT_RDD_Edmund_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy){

  if (r_m >= r_min_m && r_m <= r_max_m){
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){
      if( r_m > a0_m ){
        return AT_RDD_Katz_PowerLawER_Dpoint_Gy(r_m, alpha, r_max_m, Katz_point_coeff_Gy);
      } else {

        double d_max_Gy;
        double dEdx_MeV_cm2_g  =  dEdx_J_m / (density_kg_m3 / 10.) / MeV_to_J;
        double LET_MeV_cm2_g   =  LET_J_m / (density_kg_m3 / 10.) / MeV_to_J;
        double tmp             =  LET_MeV_cm2_g - dEdx_MeV_cm2_g;  // LET - dEdx (MeV_g_cm2)
        if(tmp > 0){
          tmp *=  (density_kg_m3 / 10.) * MeV_to_J; // [MeV/cm] -> [J/m]
          double core_kg_m     =  M_PI * gsl_pow_2(a0_m) * density_kg_m3;  // kg / m
          d_max_Gy             =  tmp / core_kg_m;                                     // J / kg = Gy
        }else{
          d_max_Gy = 1e-11;
        } // end if tmp
        return d_max_Gy;

      }
    }
  }
  return 0.;
}


inline double  AT_RDD_Geiss_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double a0_m,
    const long   norm_constant_Gy){

  if ( (r_min_m <= r_m) && (r_m <= a0_m)){
    return norm_constant_Gy;
  }
  if ((a0_m < r_m) && (r_m <= r_max_m)){
    return norm_constant_Gy * gsl_pow_2(a0_m / r_m);
  }
  return 0.;
}

inline double   AT_RDD_CucinottaPoint_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy){

  if( r_m >= r_min_m && r_m <= r_max_m ){
    // D(r)    = Dexc(r) + Ddelta(r)
    return AT_RDD_Cucinotta_Dexc_Gy(r_m, r_max_m, beta, C_norm, Katz_point_coeff_Gy) + AT_RDD_Cucinotta_Ddelta_Gy(r_m, r_max_m, beta, Katz_point_coeff_Gy);
  } else {
    return 0.0;
  }
}


////////////////////////////////////// REST //////////////////////////////////////

double AT_RDD_r_min_m(
    const double  max_electron_range_m,
    /* radial dose distribution model */
    const long    rdd_model,
    const float   rdd_parameter[]){

  double r_min_m = 0.0;
  if( rdd_model == RDD_KatzPoint || rdd_model == RDD_CucinottaPoint || rdd_model == RDD_KatzExtTarget || rdd_model == RDD_CucinottaExtTarget){
    r_min_m = (double)rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m             =  max_electron_range_m;
    }                  // If r.max < r.min, r.min = r.max
  }
  return r_min_m;
}


double AT_RDD_a0_m(
    const double  max_electron_range_m,
    const long    rdd_model,
    const float   rdd_parameter[]){

  double a0_m = 0.0;
  if( rdd_model == RDD_Geiss || rdd_model == RDD_KatzSite || rdd_model == RDD_Edmund ){
    a0_m = (double)rdd_parameter[0];
  } else if ( rdd_model == RDD_KatzExtTarget || rdd_model == RDD_CucinottaExtTarget ) {
    a0_m = (double)rdd_parameter[1];
  }
  if (max_electron_range_m <= a0_m){
    a0_m             =  max_electron_range_m;
  }                  // If r.max < r.min, r.min = r.max
  return a0_m;
}


double AT_RDD_precalculated_constant_Gy(
    const double  max_electron_range_m,
    const double  LET_MeV_cm2_g,
    const double  E_MeV_u,
    const long    particle_no,
    /* detector parameters */
    const long    material_no,
    /* radial dose distribution model */
    const long    rdd_model,
    const float   rdd_parameter[],
    /* electron range model */
    const long    er_model){

  double precalculated_constant_Gy  = 0.0;

  double beta   =  AT_beta_from_E_single( E_MeV_u );

  // Get density
  double   density_g_cm3, electron_density_m3;
  AT_get_material_data(
      material_no,
      &density_g_cm3,
      &electron_density_m3,
      NULL, NULL, NULL, NULL, NULL, NULL);
  double density_kg_m3      =  density_g_cm3 * 1000.0;

  // For those models that have r_min cut-off
  double r_min_m = AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);

  // For those models that have a0
  double a0_m = AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);

  // get alpha for some ER models
  double alpha               =  0.0;
  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){
    double wmax_MeV     =  AT_max_E_transfer_MeV_single(E_MeV_u);
    alpha                   =  1.667;
    if(wmax_MeV <= 1e-3){  // if wmax < 1keV
      alpha                 =  1.079;
    } // end if
  }

  double C_J_m                =  0.0;
  double Katz_point_coeff_Gy  =  0.0;
  if( (rdd_model == RDD_KatzPoint) || (rdd_model == RDD_KatzSite) || (rdd_model == RDD_Edmund) || (rdd_model == RDD_CucinottaPoint) || (rdd_model == RDD_KatzExtTarget)){
    // TODO shall we move from J_m and m to more reasonable units ?
    // we have for water C_J_m = 1.22e-12 and r_m usually ~ 1e-8
    // calculations in C_J_um and r_um would be more precise
    const long Z            =  AT_Z_from_particle_no_single(particle_no);
    const double Z_eff      =  AT_effective_charge_from_beta_single(beta, Z);
    C_J_m                   =  AT_RDD_Katz_C_J_m(electron_density_m3);
    Katz_point_coeff_Gy     =  AT_RDD_Katz_coeff_Gy(C_J_m, Z_eff, beta, density_kg_m3, max_electron_range_m);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Test
  if( rdd_model == RDD_Test){
    const double single_impact_fluence_cm2 = AT_single_impact_fluence_cm2_single(E_MeV_u, material_no, er_model);
    precalculated_constant_Gy = LET_MeV_cm2_g * single_impact_fluence_cm2 * MeV_to_J * 1000.0;          // LET  / track area = Norm.constant k , TODO check units
  } // end RDD_Test


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzPoint
  if( rdd_model == RDD_KatzPoint){
    precalculated_constant_Gy   =  Katz_point_coeff_Gy;
  } // end RDD_KatzPoint


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzSite
  if( rdd_model == RDD_KatzSite){
    double  LET_J_m       =  LET_MeV_cm2_g * density_g_cm3; // [MeV / cm]
    LET_J_m              *=  100.0;        // [MeV / m]
    LET_J_m              *=  MeV_to_J;     // [J/m]

    double dEdx_J_m       =  0.0;
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // calculate dEdx_MeV_cm2_g from "new" Katz RDD
      dEdx_J_m            =  AT_RDD_Katz_PowerLawER_dEdx_J_m(a0_m, max_electron_range_m, density_kg_m3, alpha, Katz_point_coeff_Gy);
    } else if (er_model == ER_ButtsKatz){ // calculate dEdx_MeV_cm2_g from "old" Katz RDD
      dEdx_J_m            =  AT_RDD_Katz_LinearER_dEdx_J_m(a0_m, max_electron_range_m, density_kg_m3, Katz_point_coeff_Gy);
    } else {
      dEdx_J_m            =  0.0;
    }
    precalculated_constant_Gy      =  dEdx_J_m;
  } // end RDD_KatzSite

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Edmund
  //TODO rewrite in better way
  if( rdd_model == RDD_Edmund){
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // This model will work only with ER_Waligorski and ER_Edmund
      // Get dEdx by simple integration from r_min_m to r_max_m
      double dEdx_J_m     =  0.0;
      double Katz_dEdx_coeff_J_m  = (double)AT_RDD_Katz_dEdx_coeff_J_m((float)max_electron_range_m, (float)density_kg_m3, (float)Katz_point_coeff_Gy);
      dEdx_J_m            =  AT_RDD_Katz_PowerLawER_dEdx_directVersion_J_m((float)alpha, (float)a0_m, (float)max_electron_range_m, (float)Katz_dEdx_coeff_J_m);
      precalculated_constant_Gy    =  dEdx_J_m;
    }
  }// end RDD_Edmund


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Geiss
  if( rdd_model == RDD_Geiss){
    double  tmp_cm        =  0.5 + log(max_electron_range_m / a0_m);
    tmp_cm               *=  2.0 * M_PI * gsl_pow_2(a0_m * m_to_cm);                      // Normalization to match with LET
    precalculated_constant_Gy      =  LET_MeV_cm2_g * MeV_g_to_J_kg / tmp_cm;                      // k = LET / tmp
  } // end RDD_Geiss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaPoint
  if( rdd_model == RDD_CucinottaPoint){ // TODO does Cucinotta RDD model work only with Tabata ER ?
    double  LET_J_m       =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
    precalculated_constant_Gy      =  AT_RDD_Cucinotta_Cnorm(r_min_m, max_electron_range_m, beta, density_kg_m3, LET_J_m, Katz_point_coeff_Gy);
    if( precalculated_constant_Gy == 0){
      printf("problem in AT_RDD_precalculated_constant_Gy\n"); // TODO handle this situation
    }
  }// end RDD_CucinottaPoint


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzExtTarget
  if( rdd_model == RDD_KatzExtTarget){
    double Katz_plateau_Gy  =  0.0;
    double r_max_m          =  GSL_MIN(a0_m, max_electron_range_m);
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
      Katz_plateau_Gy  =  AT_RDD_Katz_PowerLawER_Daverage_Gy( r_min_m, r_max_m, max_electron_range_m, alpha, Katz_point_coeff_Gy );
    } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
      Katz_plateau_Gy  =  AT_RDD_Katz_LinearER_Daverage_Gy( r_min_m, r_max_m, max_electron_range_m, Katz_point_coeff_Gy );
    }
    precalculated_constant_Gy      =  Katz_plateau_Gy;
  }// end RDD_KatzExtTarget

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaExtTarget
  if( rdd_model == RDD_CucinottaExtTarget){
    // 4. TODO
    precalculated_constant_Gy      =  0.0;
  }// end RDD_CucinottaExtTarget


  return precalculated_constant_Gy;
}


double AT_RDD_d_min_Gy(
    const double  norm_constant_Gy,
    const double  max_electron_range_m,
    /* radial dose distribution model */
    const long    rdd_model,
    const float   rdd_parameter[]){

  double d_min_Gy  = 0.0;

  if( rdd_model == RDD_Test){
    d_min_Gy              =  norm_constant_Gy;
  } // end RDD_Test

  if( rdd_model == RDD_KatzPoint || rdd_model == RDD_KatzSite || rdd_model == RDD_Edmund || rdd_model == RDD_CucinottaPoint){
    d_min_Gy              =  (double)rdd_parameter[1];
  }

  if( rdd_model == RDD_Geiss){
    double a0_m           =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    d_min_Gy              =  AT_RDD_Geiss_Gy( max_electron_range_m, 0., max_electron_range_m, a0_m, norm_constant_Gy);
  } // end RDD_Geiss

  if( rdd_model == RDD_KatzExtTarget || rdd_model == RDD_CucinottaExtTarget){
    d_min_Gy              =  (double)rdd_parameter[2];
  }

  return d_min_Gy;
}


double AT_RDD_d_max_Gy(
    const double  E_MeV_u,
    const long    particle_no,
    /* detector parameters */
    const long    material_no,
    /* radial dose distribution model */
    const long    rdd_model,
    const float   rdd_parameter[],
    /* electron range model */
    const long    er_model,
    const float   er_parameter[]){

  double d_max_Gy  = 0.0;

  const double beta   =  AT_beta_from_E_single( E_MeV_u );

  const double max_electron_range_m = AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);

  // Get density
  double   density_g_cm3, electron_density_m3;
  AT_get_material_data(
              material_no,
              &density_g_cm3,
              &electron_density_m3,
              NULL, NULL, NULL, NULL, NULL, NULL);
  const double density_kg_m3      =  density_g_cm3 * 1000.0;

  const double LET_MeV_cm2_g = AT_LET_MeV_cm2_g_single(E_MeV_u, particle_no, material_no);

  const double r_min_m   =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);

  const double a0_m      =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);

  const double precalculated_constant_Gy = AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model);

  // get alpha for some ER models
  double alpha               =  0.0;
  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){
    double wmax_MeV     =  AT_max_E_transfer_MeV_single(E_MeV_u);
    alpha                   =  1.667;
    if(wmax_MeV <= 1e-3){  // if wmax < 1keV
      alpha                 =  1.079;
    } // end if
  }

  double C_J_m                =  0.0;
  double Katz_point_coeff_Gy  =  0.0;
  if( (rdd_model == RDD_KatzSite) || (rdd_model == RDD_Edmund) || (rdd_model == RDD_CucinottaPoint) || (rdd_model == RDD_KatzExtTarget) || (rdd_model == RDD_CucinottaExtTarget)){
    // TODO shall we move from J_m and m to more reasonable units ?
    // we have for water C_J_m = 1.22e-12 and r_m usually ~ 1e-8
    // calculations in C_J_um and r_um would be more precise
    const long Z            =  AT_Z_from_particle_no_single(particle_no);
    const double Z_eff      =  AT_effective_charge_from_beta_single(beta, Z);
    C_J_m                   =  AT_RDD_Katz_C_J_m(electron_density_m3);
    Katz_point_coeff_Gy     =  AT_RDD_Katz_coeff_Gy(C_J_m, Z_eff, beta, density_kg_m3, max_electron_range_m);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Test
  if( rdd_model == RDD_Test){
    d_max_Gy              = precalculated_constant_Gy;
  } // end RDD_Test

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzPoint
  if( rdd_model == RDD_KatzPoint){
    double Katz_point_coeff_Gy  =  precalculated_constant_Gy;
    d_max_Gy              =  AT_RDD_KatzPoint_Gy(r_min_m, r_min_m, max_electron_range_m, er_model, alpha, Katz_point_coeff_Gy);
  } // end RDD_KatzPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzSite
  if( rdd_model == RDD_KatzSite){
    double  LET_J_m       =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
    double dEdx_J_m       =  precalculated_constant_Gy;
    d_max_Gy              =  AT_RDD_KatzSite_Gy(0.0, 0.0, max_electron_range_m, a0_m, er_model, alpha, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
  } // end RDD_KatzSite

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Edmund
  if( rdd_model == RDD_Edmund){
    double  LET_J_m       =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
    double dEdx_J_m       =  precalculated_constant_Gy;
    d_max_Gy              =  AT_RDD_Edmund_Gy(0.0, 0.0, max_electron_range_m, a0_m, er_model, alpha, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
  }// end RDD_Edmund


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Geiss
  if( rdd_model == RDD_Geiss){
    d_max_Gy              =  precalculated_constant_Gy;
  } // end RDD_Geiss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaPoint
  if( rdd_model == RDD_CucinottaPoint){
    d_max_Gy              =  AT_RDD_CucinottaPoint_Gy(r_min_m, r_min_m, max_electron_range_m, beta, precalculated_constant_Gy, Katz_point_coeff_Gy);
  }// end RDD_CucinottaPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzExtTarget
  if( rdd_model == RDD_KatzExtTarget){
    double Katz_plateau_Gy     =  precalculated_constant_Gy;
    double Katz_point_r_min_m  =  r_min_m;
    d_max_Gy  =  AT_RDD_ExtendedTarget_KatzPoint_Gy( 0.0, a0_m, er_model, Katz_point_r_min_m, max_electron_range_m, alpha, Katz_plateau_Gy, Katz_point_coeff_Gy);
  }// end RDD_KatzExtTarget

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaExtTarget
  if( rdd_model == RDD_CucinottaExtTarget){
    const long   n_tmp     =  1;
    const float rdd_basic_parameter[] = {r_min_m, rdd_parameter[2]};
    const float r_min_m_float = (float)r_min_m;
    float d_max_Gy_float;
    AT_RDD_ExtendedTarget_Gy(n_tmp, &r_min_m_float, a0_m, E_MeV_u, particle_no, material_no, RDD_CucinottaPoint, rdd_basic_parameter, er_model, er_parameter, &d_max_Gy_float);
    d_max_Gy = (double)d_max_Gy_float;
  }// end RDD_CucinottaExtTarget

  return d_max_Gy;
}


void AT_RDD_f1_parameters(  /* radiation field parameters */
    const double  E_MeV_u,
    const long    particle_no,
    /* detector parameters */
    const long    material_no,
    /* radial dose distribution model */
    const long    rdd_model,
    const float   rdd_parameter[],
    /* electron range model */
    const long    er_model,
    const float   er_parameter[],
    /* calculated parameters */
    float         f1_parameters[])
{
  double LET_MeV_cm2_g              =  0.0;
  double max_electron_range_m       =  0.0;
  double r_min_m                    =  0.0;
  double single_impact_fluence_cm2  =  0.0;
  double single_impact_dose_Gy      =  0.0;
  double norm_constant_Gy           =  0.0;
  double d_min_Gy                   =  0.0;
  double d_max_Gy                   =  0.0;

  ///////////////////////////////////////////////////////////////////////////////
  // PARAMETER 0: Get the LET (same for all models)
  LET_MeV_cm2_g = AT_LET_MeV_cm2_g_single(E_MeV_u, particle_no, material_no);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 2: Get the maximum electron range (same for all RDD models)
  max_electron_range_m = AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 1: Get the r_min
  r_min_m   = AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 6: Get the single impact fluence (same for all RDD models)
  single_impact_fluence_cm2 = AT_single_impact_fluence_cm2_single(E_MeV_u, material_no, er_model);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 7: Get the single impact dose (same for all RDD models)
  single_impact_dose_Gy = AT_single_impact_dose_Gy_single(LET_MeV_cm2_g, single_impact_fluence_cm2);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 5: Get normalization constant
  norm_constant_Gy = AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 3: Get minimum dose d_min_Gy (f1_parameters[3])
  d_min_Gy = AT_RDD_d_min_Gy( norm_constant_Gy, max_electron_range_m, rdd_model, rdd_parameter);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 4: Get minimum dose d_min_Gy (f1_parameters[4])
  d_max_Gy = AT_RDD_d_max_Gy( E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, er_parameter);

  // write data to output table
  f1_parameters[0]  =  (float)LET_MeV_cm2_g;
  f1_parameters[1]  =  (float)r_min_m;
  f1_parameters[2]  =  (float)max_electron_range_m;
  f1_parameters[3]  =  (float)d_min_Gy;
  f1_parameters[4]  =  (float)d_max_Gy;
  f1_parameters[5]  =  (float)norm_constant_Gy;
  f1_parameters[6]  =  (float)single_impact_fluence_cm2;
  f1_parameters[7]  =  (float)single_impact_dose_Gy;
}

void AT_D_RDD_Gy( const long  n,
    const float   r_m[],
    /* radiation field parameters */
    const float   E_MeV_u,
    const long    particle_no,
    /* detector parameters */
    const long    material_no,
    /* radial dose distribution model */
    const long    rdd_model,
    const float   rdd_parameter[],
    /* electron range model */
    const long    er_model,
    const float   er_parameter[],
    float         D_RDD_Gy[])
{
  /********************************************************
   ********* CALCULATION BEFORE PARTICLE LOOP *************
   *******************************************************/

  const double LET_MeV_cm2_g          =  AT_LET_MeV_cm2_g_single(E_MeV_u, particle_no, material_no);

  const double max_electron_range_m   =  AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);

  const double r_min_m                =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);

  const double precalculated_constant_Gy  =  AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model);

  const double d_min_Gy               =  AT_RDD_d_min_Gy( precalculated_constant_Gy, max_electron_range_m, rdd_model, rdd_parameter);

  const double a0_m                   =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);

  // Get material data
  double   density_g_cm3, electron_density_m3;
  AT_get_material_data(
              material_no,
              &density_g_cm3,
              &electron_density_m3,
              NULL, NULL, NULL, NULL, NULL, NULL);
  double density_kg_m3      =  density_g_cm3 * 1000.0;

  double beta    =  AT_beta_from_E_single( (double)E_MeV_u );

  // get alpha for some ER models
  double alpha               =  0.0;
  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){
    double wmax_MeV     =  AT_max_E_transfer_MeV_single((double)E_MeV_u);
    alpha                   =  1.667;
    if(wmax_MeV <= 1e-3){  // if wmax < 1keV
      alpha                 =  1.079;
    } // end if
  }

  double Katz_point_coeff_Gy  =  0.0;
  if( (rdd_model == RDD_KatzPoint) || (rdd_model == RDD_KatzSite) || (rdd_model == RDD_Edmund) || (rdd_model == RDD_CucinottaPoint)){
    // TODO shall we move from J_m and m to more reasonable units ?
    // we have for water C_J_m = 1.22e-12 and r_m usually ~ 1e-8
    // calculations in C_J_um and r_um would be more precise
    const long  Z           =  AT_Z_from_particle_no_single(particle_no);
    const double Z_eff      =  AT_effective_charge_from_beta_single(beta, Z);
    const double C_J_m      =  AT_RDD_Katz_C_J_m(electron_density_m3);
    Katz_point_coeff_Gy     =  AT_RDD_Katz_coeff_Gy(C_J_m, Z_eff, beta, density_kg_m3, max_electron_range_m);
  }


  /********************************************************
   *************** LOOP OVER DISTANCE VECTOR **************
   *******************************************************/
  long     i;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Test
  if( rdd_model == RDD_Test){
    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]         =  (float)AT_RDD_Test_Gy((double)r_m[i], 0., max_electron_range_m, precalculated_constant_Gy);
      // Cut-off low doses at low doses not necessary here
    }
  }// end RDD_Test

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzPoint
  if( rdd_model == RDD_KatzPoint){ // RDD formula will be determined by form of ER model
    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]     =  (float)AT_RDD_KatzPoint_Gy((double)r_m[i], r_min_m, max_electron_range_m, er_model, alpha, Katz_point_coeff_Gy);
      D_RDD_Gy[i]     =  fmaxf(D_RDD_Gy[i], (double)d_min_Gy);          // Cut-off low doses, necessary in SPIFF
    } // end for
  }// end RDD_KatzPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzSite
  if( rdd_model == RDD_KatzSite){
    const double LET_J_m   =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J;     // convert LET_MeV_cm2_g to LET_J_m
    const double dEdx_J_m  = precalculated_constant_Gy;                                      // take dEdx averaged on outer shell from norm_constant

    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]     =  (float)AT_RDD_KatzSite_Gy((double)r_m[i], 0.0, max_electron_range_m, a0_m, er_model, alpha, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
      D_RDD_Gy[i]     =  fmaxf(D_RDD_Gy[i], (double)d_min_Gy);          // Cut-off low doses, necessary in SPIFF
    } // end for
  }// end RDD_KatzSite

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Edmund
  if(rdd_model == RDD_Edmund){
    const double LET_J_m   =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J;     // convert LET_MeV_cm2_g to LET_J_m
    const double dEdx_J_m  = precalculated_constant_Gy;                                      // take dEdx averaged on outer shell from norm_constant

    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]     =  (float)AT_RDD_Edmund_Gy((double)r_m[i], 0.0, max_electron_range_m, a0_m, er_model, alpha, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
      D_RDD_Gy[i]     =  fmaxf(D_RDD_Gy[i], (float)d_min_Gy);          // Cut-off low doses, necessary in SPIFF
    } // end for

  }// end RDD_KatzSite

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Geiss
  if( rdd_model == RDD_Geiss){
    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]     =  (float)AT_RDD_Geiss_Gy((double)r_m[i], 0., max_electron_range_m, a0_m, precalculated_constant_Gy);
      // Cut-off at low doses at low doses not necessary here
    }
  }// end RDD_Geiss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaPoint
  if( rdd_model == RDD_CucinottaPoint){
    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]     =  (float)AT_RDD_CucinottaPoint_Gy((double)r_m[i], r_min_m, max_electron_range_m, beta, precalculated_constant_Gy, Katz_point_coeff_Gy);
      D_RDD_Gy[i]     =  fmaxf(D_RDD_Gy[i], (double)d_min_Gy);          // Cut-off low doses, necessary in SPIFF
    }
  }// end RDD_CucinottaPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzExtTarget
  if( rdd_model == RDD_KatzExtTarget){
    const float rdd_basic_parameter[] = {(float)r_min_m, rdd_parameter[2]};
    AT_RDD_ExtendedTarget_Gy(n, r_m, a0_m, E_MeV_u, particle_no, material_no, RDD_KatzPoint, rdd_basic_parameter, er_model, er_parameter, D_RDD_Gy);
  }// end RDD_KatzExtTarget

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaExtTarget
  if( rdd_model == RDD_CucinottaExtTarget){
    const float rdd_basic_parameter[] = {(double)r_min_m, rdd_parameter[2]};
    AT_RDD_ExtendedTarget_Gy(n, r_m, a0_m, E_MeV_u, particle_no, material_no, RDD_CucinottaPoint, rdd_basic_parameter, er_model, er_parameter, D_RDD_Gy);
  } // end RDD_CucinottaExtTarget

}


inline double  AT_inverse_RDD_Test_m( const double D_Gy,
    const double r_max_m){
  if( D_Gy > 0)
    return r_max_m;
  return 0.0;
}

inline double  AT_inverse_RDD_Geiss_m( const double D_Gy,
    const double d_min_Gy,
    const double d_max_Gy,
    const double a0_m,
    const double norm_constant_Gy){
  if ((d_min_Gy <= D_Gy) && (D_Gy <= d_max_Gy)){
    return a0_m * sqrt(norm_constant_Gy / D_Gy);
  }
  if(  D_Gy < d_min_Gy ){ // TODO temporary fix to make SPIFF working
    return a0_m * sqrt(norm_constant_Gy / d_min_Gy);
  }
  return 0.0;
}

double  AT_inverse_RDD_KatzPoint_LinearER_m( const double D_Gy,
    const double r_min_m,
    const double r_max_m,
    const double Katz_point_coeff_Gy){
  double x = 0.5 * (1.0 + sqrt(1.0 + 4.0 * D_Gy / Katz_point_coeff_Gy));
  return r_max_m / x;
}

float AT_inverse_RDD_KatzPoint_PowerLawER_solver_function_Gy( const float r_m , void * params ){
  AT_inverse_RDD_KatzPoint_PowerLawER_parameters* RDD_parameters = (AT_inverse_RDD_KatzPoint_PowerLawER_parameters*)(params);

  const double  D_Gy                  =  RDD_parameters->D_Gy;
  const double  r_min_m               =  RDD_parameters->r_min_m;
  const double  r_max_m               =  RDD_parameters->r_max_m;
  const double  alpha                 =  RDD_parameters->alpha;
  const double  Katz_point_coeff_Gy   =  RDD_parameters->Katz_point_coeff_Gy;

  if( r_m < (double)r_min_m){
    return -(float)D_Gy;
  } else {
    return (float)AT_RDD_Katz_PowerLawER_Dpoint_Gy((double)r_m, alpha, r_max_m, Katz_point_coeff_Gy) - (float)D_Gy;
  }
}


double  AT_inverse_RDD_KatzPoint_m( const double D_Gy,
    const double r_min_m,
    const double r_max_m,
    const long   er_model,
    const double alpha,
    const double Katz_point_coeff_Gy){

  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
    float  solver_accuracy  =  1e-13f;

    AT_inverse_RDD_KatzPoint_PowerLawER_parameters RDD_parameters;

    RDD_parameters.D_Gy                 =  D_Gy;
    RDD_parameters.r_min_m              =  r_min_m;
    RDD_parameters.r_max_m              =  r_max_m;
    RDD_parameters.alpha                =  alpha;
    RDD_parameters.Katz_point_coeff_Gy  =  Katz_point_coeff_Gy;
    double r_m =  zriddr(AT_inverse_RDD_KatzPoint_PowerLawER_solver_function_Gy,
            (void*)(&RDD_parameters),
            (float)r_min_m,
            (float)r_max_m,
            solver_accuracy);

    return r_m;

  } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
    return AT_inverse_RDD_KatzPoint_LinearER_m(D_Gy, r_min_m, r_max_m, Katz_point_coeff_Gy);
  }

  return 0.0;
}

float AT_D_RDD_Gy_solver( const float r , void * params ){
  AT_D_RDD_Gy_parameters* params_struct = (AT_D_RDD_Gy_parameters*)(params);
  (*params_struct).n = 1;
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
void AT_r_RDD_m  ( const long  n,
    const float   D_RDD_Gy[],
    /* radiation field parameters */
    const float   E_MeV_u,
    const long    particle_no,
    /* detector parameters */
    const long    material_no,
    /* radial dose distribution model */
    const long    rdd_model,
    const float   rdd_parameter[],
    /* electron range model */
    const long    er_model,
    const float   er_parameter[],
    float         r_RDD_m[])
{
  long     i;

  const double LET_MeV_cm2_g  = AT_LET_MeV_cm2_g_single(E_MeV_u, particle_no, material_no);

  const double max_electron_range_m = AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);

  const double r_min_m        = AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);

  const double precalculated_constant_Gy = AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model);

  const double d_min_Gy       = AT_RDD_d_min_Gy( precalculated_constant_Gy, max_electron_range_m, rdd_model, rdd_parameter);

  const double d_max_Gy       = AT_RDD_d_max_Gy( E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, er_parameter);

  const double a0_m           =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);

  // Get material data
  const double density_g_cm3  =  AT_density_g_cm3_from_material_no(material_no);
  const double density_kg_m3  =  density_g_cm3 * 1000.0;

  const double electron_density_m3 = AT_electron_density_m3_from_material_no(material_no);

  ///////////////////////////////////////////////////////////////////////////////////
  // Get beta, Z and Zeff
  const double beta   =  AT_beta_from_E_single( (double)E_MeV_u );
  const long   Z      =  AT_Z_from_particle_no_single(particle_no);
  const double Z_eff  =  AT_effective_charge_from_beta_single(beta, Z);

  const double C_J_m      =  AT_RDD_Katz_C_J_m(electron_density_m3);
  const double Katz_point_coeff_Gy = AT_RDD_Katz_coeff_Gy(C_J_m, Z_eff, beta, density_kg_m3, max_electron_range_m);

  double alpha           =  0.0;
  if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
    alpha               =  1.667;
    double wmax_MeV     =  AT_max_E_transfer_MeV_single(E_MeV_u);
    if(wmax_MeV <= 1e-3){
      alpha             =  1.079;
    }
  }


  if( rdd_model == RDD_Test){
    // Loop over all doses given
    for (i = 0; i < n; i++){
      r_RDD_m[i]    =  (float)AT_inverse_RDD_Test_m( (double)D_RDD_Gy[i], max_electron_range_m);
    }
  }// end RDD_Test

  if( rdd_model == RDD_Geiss){
    // Loop over all doses given
    for (i = 0; i < n; i++){
      r_RDD_m[i]    =  (float)AT_inverse_RDD_Geiss_m( (double)D_RDD_Gy[i], d_min_Gy, d_max_Gy, a0_m, precalculated_constant_Gy);
    }
  }// end RDD_Geiss


  if( rdd_model == RDD_KatzPoint){
    // Loop over all doses given
    const double Katz_point_coeff_Gy = precalculated_constant_Gy;
    for (i = 0; i < n; i++){
      if( (D_RDD_Gy[i] >= d_min_Gy) && (D_RDD_Gy[i] <= d_max_Gy) ){
        r_RDD_m[i]    =  (float)AT_inverse_RDD_KatzPoint_m( (double)D_RDD_Gy[i], r_min_m, max_electron_range_m, er_model, alpha, Katz_point_coeff_Gy);
      } else {
        r_RDD_m[i]    =  0.0f;
      }
    }
  }// end RDD_KatzPoint

  if( rdd_model == RDD_KatzExtTarget){
    // Loop over all doses given
    double Katz_plateau_Gy  =  precalculated_constant_Gy;
    for (i = 0; i < n; i++){
      if( (D_RDD_Gy[i] >= d_min_Gy) && (D_RDD_Gy[i] <= d_max_Gy) ){
        r_RDD_m[i]    =  (float)AT_inverse_RDD_ExtendedTarget_KatzPoint_m( (double)D_RDD_Gy[i], r_min_m, max_electron_range_m, a0_m, er_model, alpha, Katz_plateau_Gy, Katz_point_coeff_Gy);
      } else {
        r_RDD_m[i]    =  0.0f;
      }
    }
  }// end RDD_KatzExtTarget


  long     n_tmp      = 1;
  if( (rdd_model == RDD_KatzSite) || (rdd_model == RDD_Edmund)){
    AT_D_RDD_Gy_parameters* params = (AT_D_RDD_Gy_parameters*)calloc(1,sizeof(AT_D_RDD_Gy_parameters));
    params->n             = 1;
    params->E_MeV_u       = (float)(E_MeV_u);
    params->particle_no   = (long)(particle_no);
    params->material_no   = (long)(material_no);
    params->rdd_model     = (long)(rdd_model);
    params->rdd_parameter = (float*)rdd_parameter;
    params->er_model      = (long)(er_model);
    params->er_parameter  = (float*)er_parameter;
    params->D_RDD_Gy      = (float*)calloc(1,sizeof(float));
    // Loop over all doses given
    float dev = 0.01f;
    float critical_d_Gy    = 0.0f;
    float critical_r_m     = 0.0f;
    float inv2_d_Gy        = 0.0f;
    float inv2_r_m         = 0.0f;
    if(rdd_model == RDD_KatzSite || rdd_model == RDD_Edmund){
      critical_r_m   = (float)a0_m * (1.0f + 1e-8f);
      inv2_r_m       = fmaxf((float)a0_m, (float)max_electron_range_m * dev);
    }
    if(rdd_model == RDD_KatzSite || rdd_model == RDD_Edmund){
      AT_D_RDD_Gy  (  n_tmp,                    // Use D(r) to find dose at jump of D_Site
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
      AT_D_RDD_Gy  (  n_tmp,                    // Use D(r) to find dose at jump of D_Site
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

    for (i = 0; i < n; i++){
      params->D0 = D_RDD_Gy[i];
      // if D is outside the definition, return -1
      r_RDD_m[i]        =  -1.0f;
      float  solver_accuracy  =  1e-13f;
      if (((float)d_min_Gy <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= (float)d_max_Gy)){
        if(D_RDD_Gy[i] >= critical_d_Gy){
          if(rdd_model == RDD_KatzSite || rdd_model == RDD_Edmund){
            r_RDD_m[i] = r_min_m;
          }
        }
          if(rdd_model == RDD_KatzSite || rdd_model == RDD_Edmund){
//          if(D_RDD_Gy[i] < critical_d_Gy && D_RDD_Gy[i] >= inv2_d_Gy){
//            r_RDD_m[i] = sqrt((float)Katz_point_coeff_Gy / D_RDD_Gy[i]) * (float)max_electron_range_m;}
//          if(D_RDD_Gy[i] < inv2_d_Gy){
          if(D_RDD_Gy[i] < critical_d_Gy){
            r_RDD_m[i] = zriddr(AT_D_RDD_Gy_solver,
                    (void*)params,
                    (float)a0_m,
                     (float)max_electron_range_m,
                    solver_accuracy);
          }
        }else{
          r_RDD_m[i] = zriddr(AT_D_RDD_Gy_solver,
                  (void*)params,
                  (float)r_min_m,
                  (float)max_electron_range_m,
                  solver_accuracy);

        }
      }
    }
    free(params->D_RDD_Gy);
    free(params);
  }// end RDD_KatzSite
}


double AT_D_RDD_Gy_int( double  r_m,
    void*   params){
  float D_Gy;
  long  n_tmp              = 1;
  AT_P_RDD_parameters* par = (AT_P_RDD_parameters*)params;
  float fr_m               = (float)r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  n_tmp,
      &fr_m,
      *(par->E_MeV_u),
      *(par->particle_no),
      *(par->material_no),
      *(par->rdd_model),
      par->rdd_parameters,
      *(par->er_model),
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
  AT_D_RDD_Gy  (  n_tmp,
      &fr_m,
      *(par->E_MeV_u),
      *(par->particle_no),
      *(par->material_no),
      *(par->rdd_model),
      par->rdd_parameters,
      *(par->er_model),
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
