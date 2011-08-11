/**
 * @brief Radial Dose Distribution models
 */

/*
 *    AT_RDD.c
 *    ========
 *
 *    Created on: 28.07.2009
 *    Creator: greilich
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

#include "AT_RDD.h"


long AT_RDD_index_from_RDD_number( const long RDD_number ){
  long  index           =  -1;
  long  number_of_RDDs  =  1;
  find_elements_int(  &RDD_number,
      number_of_RDDs,
      AT_RDD_Data.RDD_no,
      AT_RDD_Data.n,
      &index);   // TODO replace call to find_elements_int by call to simpler function which will find the index just for one argument
  return index;
}


int AT_RDD_name_from_number(const long RDD_no, char* RDD_name){
  long  index = AT_RDD_index_from_RDD_number( RDD_no );

  if( index != -1){
    strcpy(RDD_name, AT_RDD_Data.RDD_name[index]);
  } else {
    strcpy(RDD_name,"*** invalid choice of RDD ***");
    return -1;
  }
  return AT_Success;
}


long AT_RDD_number_from_name(const char* RDD_name){
  // find look-up index for material name in material data table
  long  match;
  const long n_tmp = 1;

  find_elements_char(  &RDD_name,
      n_tmp,
      AT_RDD_Data.RDD_name,
      AT_RDD_Data.n,
      &match);

  if( match != -1){
    return AT_RDD_Data.RDD_no[match];
  } else {
    return -1;
  }
}


int AT_RDD_number_of_parameters( const long RDD_model){
  long  index = AT_RDD_index_from_RDD_number( RDD_model );
  if( index == -1){
    printf("RDD no %ld not found\n", RDD_model);
    return 0;
  }
  return AT_RDD_Data.n_parameters[index];
}


////////////////////////////////////// REST //////////////////////////////////////

double AT_RDD_r_min_m(
    const double   max_electron_range_m,
    const long     rdd_model,
    const double   rdd_parameter[]){

  double r_min_m = 0.0;
  if( rdd_model == RDD_KatzPoint || rdd_model == RDD_CucinottaPoint || rdd_model == RDD_KatzExtTarget || rdd_model == RDD_CucinottaExtTarget){
    r_min_m = rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m      =  max_electron_range_m;
    }                  // If r.max < r.min, r.min = r.max
  }
  return r_min_m;
}


double AT_RDD_a0_m(
    const double   max_electron_range_m,
    const long     rdd_model,
    const double   rdd_parameter[]){

  double a0_m = 0.0;
  if( rdd_model == RDD_Geiss || rdd_model == RDD_KatzSite){
    a0_m = rdd_parameter[0];
  } else if ( rdd_model == RDD_KatzExtTarget || rdd_model == RDD_CucinottaExtTarget ) {
    a0_m = rdd_parameter[1];
  } else {
    a0_m = 0.0;
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
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model){

  double precalculated_constant_Gy  = 0.0;

  // Get density // TODO move it down
  double density_g_cm3  =  AT_density_g_cm3_from_material_no( material_no );
  double density_kg_m3  =  density_g_cm3 * 1000.0;

  double Katz_point_coeff_Gy  =  0.0; // TODO move it down
  if( (rdd_model == RDD_KatzPoint) || (rdd_model == RDD_KatzSite) || (rdd_model == RDD_CucinottaPoint) || (rdd_model == RDD_KatzExtTarget) || (rdd_model == RDD_CucinottaExtTarget)){
    Katz_point_coeff_Gy     =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Test
  if( rdd_model == RDD_Test){
    const double single_impact_fluence_cm2 = AT_single_impact_fluence_cm2_single(E_MeV_u, material_no, er_model);
    precalculated_constant_Gy = LET_MeV_cm2_g * single_impact_fluence_cm2 * MeV_to_J * 1000.0;          // LET  / track area = Norm.constant k
  } // end RDD_Test


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzPoint
  if( rdd_model == RDD_KatzPoint){
    precalculated_constant_Gy   =  Katz_point_coeff_Gy;
  } // end RDD_KatzPoint


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzSite
  if( rdd_model == RDD_KatzSite){
    double       dEdx_J_m =  0.0;
    const double a0_m     = AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // calculate dEdx_MeV_cm2_g from "new" Katz RDD
      const double alpha  =  AT_ER_PowerLaw_alpha(E_MeV_u);
      dEdx_J_m            =  AT_RDD_Katz_PowerLawER_dEdx_J_m(a0_m, max_electron_range_m, density_kg_m3, alpha, Katz_point_coeff_Gy);
    } else if (er_model == ER_ButtsKatz){ // calculate dEdx_MeV_cm2_g from "old" Katz RDD
      dEdx_J_m            =  AT_RDD_Katz_LinearER_dEdx_J_m(a0_m, max_electron_range_m, density_kg_m3, Katz_point_coeff_Gy);
    } else {
      dEdx_J_m            =  0.0;
    }
    precalculated_constant_Gy      =  dEdx_J_m;
  } // end RDD_KatzSite


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Geiss
  if( rdd_model == RDD_Geiss){
    const double  a0_m    =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);   // TODO rename tmp_cm to something more reasonable
    double  tmp_cm        =  0.5 + log(max_electron_range_m / a0_m);
    tmp_cm               *=  2.0 * M_PI * gsl_pow_2(a0_m * m_to_cm);                      // Normalization to match with LET
    precalculated_constant_Gy      =  LET_MeV_cm2_g * MeV_g_to_J_kg / tmp_cm;             // k = LET / tmp
  } // end RDD_Geiss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaPoint
  if( rdd_model == RDD_CucinottaPoint){
    const double  r_min_m  = AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double  beta     =  AT_beta_from_E_single( E_MeV_u );
    const double  LET_J_m  =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
    precalculated_constant_Gy      =  AT_RDD_Cucinotta_Cnorm(r_min_m, max_electron_range_m, beta, density_kg_m3, LET_J_m, Katz_point_coeff_Gy);
    if( precalculated_constant_Gy == 0){
      printf("problem in AT_RDD_precalculated_constant_Gy\n"); // TODO handle this situation
    }
  }// end RDD_CucinottaPoint


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzExtTarget
  if( rdd_model == RDD_KatzExtTarget){
    double Katz_plateau_Gy  =  0.0;
    const double r_min_m    =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double a0_m       =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double r_max_m    =  GSL_MIN(a0_m, max_electron_range_m);
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
      double alpha     =  AT_ER_PowerLaw_alpha( E_MeV_u );
      Katz_plateau_Gy  =  AT_RDD_Katz_PowerLawER_Daverage_Gy( r_min_m, r_max_m, max_electron_range_m, alpha, Katz_point_coeff_Gy );
    } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
      Katz_plateau_Gy  =  AT_RDD_Katz_LinearER_Daverage_Gy( r_min_m, r_max_m, max_electron_range_m, Katz_point_coeff_Gy );
    }
    precalculated_constant_Gy      =  Katz_plateau_Gy;
  }// end RDD_KatzExtTarget

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaExtTarget
  if( rdd_model == RDD_CucinottaExtTarget){
    const double  r_min_m          =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double  beta             =  AT_beta_from_E_single( E_MeV_u );
    const double  LET_J_m          =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
    precalculated_constant_Gy      =  AT_RDD_Cucinotta_Cnorm( r_min_m, max_electron_range_m, beta, density_kg_m3, LET_J_m, Katz_point_coeff_Gy);
  }// end RDD_CucinottaExtTarget


  return precalculated_constant_Gy;
}


double AT_RDD_d_min_Gy(
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const double  precalculated_constant_Gy){

  double d_min_Gy  = 0.0;

  if( rdd_model == RDD_Test){
    d_min_Gy              =  precalculated_constant_Gy;
  } // end RDD_Test

  if( rdd_model == RDD_KatzPoint ){
    d_min_Gy              =  rdd_parameter[1];
  } // end RDD_Test

  if( rdd_model == RDD_KatzSite || rdd_model == RDD_CucinottaPoint){
    d_min_Gy              =  rdd_parameter[1];
  }

  if( rdd_model == RDD_Geiss){
    double max_electron_range_m = AT_max_electron_range_m(E_MeV_u, material_no, er_model);
    const double a0_m     =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    d_min_Gy              =  AT_RDD_Geiss_Gy( max_electron_range_m, 0., max_electron_range_m, a0_m, precalculated_constant_Gy);
  } // end RDD_Geiss

  if( rdd_model == RDD_KatzExtTarget || rdd_model == RDD_CucinottaExtTarget){
    d_min_Gy              =  rdd_parameter[2];
  }

  return d_min_Gy;
}


double AT_RDD_d_max_Gy(
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no){

  // TODO implement handling of non-compatible er and rdd

  double d_max_Gy  = 0.0;

  const double max_electron_range_m      =  AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);
  const double LET_MeV_cm2_g             =  AT_Stopping_Power_MeV_cm2_g_single( stopping_power_source_no, E_MeV_u, particle_no, material_no);
  const double precalculated_constant_Gy =  AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model);

  // Get density // TODO move it down
  double density_g_cm3  =  AT_density_g_cm3_from_material_no( material_no );
  double density_kg_m3  =  density_g_cm3 * 1000.0;

  double Katz_point_coeff_Gy  =  0.0; // TODO move it down
  if( (rdd_model == RDD_KatzSite) || (rdd_model == RDD_CucinottaPoint) || (rdd_model == RDD_KatzExtTarget) || (rdd_model == RDD_CucinottaExtTarget)){
    Katz_point_coeff_Gy     =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Test
  if( rdd_model == RDD_Test){
    d_max_Gy              = precalculated_constant_Gy;
  } // end RDD_Test

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzPoint
  if( rdd_model == RDD_KatzPoint){
    const double Katz_point_coeff_Gy  =  precalculated_constant_Gy;
    const double alpha                =  AT_ER_PowerLaw_alpha( E_MeV_u );
    const double r_min_m              =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
    d_max_Gy              =  AT_RDD_KatzPoint_Gy(r_min_m, r_min_m, max_electron_range_m, er_model, alpha, Katz_point_coeff_Gy);
  } // end RDD_KatzPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzSite
  if( rdd_model == RDD_KatzSite){
    const double  a0_m      =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double  LET_J_m   =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
    const double  dEdx_J_m  =  precalculated_constant_Gy;
    const double  alpha     =  AT_ER_PowerLaw_alpha( E_MeV_u );
    d_max_Gy                =  AT_RDD_KatzSite_Gy(0.0, 0.0, max_electron_range_m, a0_m, er_model, alpha, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
  } // end RDD_KatzSite

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Geiss
  if( rdd_model == RDD_Geiss){
    d_max_Gy              =  precalculated_constant_Gy;
  } // end RDD_Geiss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaPoint
  if( rdd_model == RDD_CucinottaPoint){
    const double beta     =  AT_beta_from_E_single( E_MeV_u );
    const double r_min_m  =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
    d_max_Gy              =  AT_RDD_CucinottaPoint_Gy(r_min_m, r_min_m, max_electron_range_m, beta, precalculated_constant_Gy, Katz_point_coeff_Gy);
  }// end RDD_CucinottaPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzExtTarget
  if( rdd_model == RDD_KatzExtTarget){
    const double a0_m                =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double Katz_plateau_Gy     =  precalculated_constant_Gy;
    const double Katz_point_r_min_m  =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double alpha               =  AT_ER_PowerLaw_alpha( E_MeV_u );
    d_max_Gy  =  AT_RDD_ExtendedTarget_KatzPoint_Gy( 0.0, a0_m, er_model, Katz_point_r_min_m, max_electron_range_m, alpha, Katz_plateau_Gy, Katz_point_coeff_Gy);
  }// end RDD_KatzExtTarget

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaExtTarget
  if( rdd_model == RDD_CucinottaExtTarget){
    const double beta                  =  AT_beta_from_E_single( E_MeV_u );
    const double a0_m                  =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double C_norm                =  precalculated_constant_Gy;
    const double Katz_point_r_min_m    =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);

    const double r_max_m               =  GSL_MIN(a0_m, max_electron_range_m);
    double Cucinotta_plateau_Gy        =  AT_RDD_Cucinotta_Ddelta_average_Gy( Katz_point_r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
    Cucinotta_plateau_Gy              +=  C_norm * AT_RDD_Cucinotta_Dexc_average_Gy( Katz_point_r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
    d_max_Gy   =  AT_RDD_ExtendedTarget_CucinottaPoint_Gy( 0.0, a0_m, Katz_point_r_min_m, max_electron_range_m, beta, Katz_point_coeff_Gy, C_norm, Cucinotta_plateau_Gy);
  }// end RDD_CucinottaExtTarget

  return d_max_Gy;
}


void AT_RDD_f1_parameters_single_field(
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no,
    double        f1_parameters[])
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
  LET_MeV_cm2_g = AT_Stopping_Power_MeV_cm2_g_single( stopping_power_source_no, E_MeV_u, particle_no, material_no);

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
  d_min_Gy = AT_RDD_d_min_Gy( E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, norm_constant_Gy);

  ////////////////////////////////////////////////////////////////////////////////
  // PARAMETER 4: Get minimum dose d_min_Gy (f1_parameters[4])
  d_max_Gy = AT_RDD_d_max_Gy( E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, stopping_power_source_no);

  // write data to output table
  f1_parameters[0]  =  LET_MeV_cm2_g;
  f1_parameters[1]  =  r_min_m;
  f1_parameters[2]  =  max_electron_range_m;
  f1_parameters[3]  =  d_min_Gy;
  f1_parameters[4]  =  d_max_Gy;
  f1_parameters[5]  =  norm_constant_Gy;
  f1_parameters[6]  =  single_impact_fluence_cm2;
  f1_parameters[7]  =  single_impact_dose_Gy;
}


void AT_RDD_f1_parameters_mixed_field(
    const long    n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no,
    double        f1_parameters[]){

  long  i;
  for (i = 0; i < n; i++){
    /* get RDD parameters for all particles and energies */
    AT_RDD_f1_parameters_single_field(  E_MeV_u[i],
        particle_no[i],
        material_no,
        rdd_model,
        rdd_parameter,
        er_model,
        stopping_power_source_no,
        &f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH]);
  }
}


int AT_D_RDD_Gy( const long  n,
    const double  r_m[],
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no,
    double        D_RDD_Gy[])
{
  /********************************************************
   ********* CALCULATION BEFORE PARTICLE LOOP *************
   *******************************************************/

  const double LET_MeV_cm2_g          =  AT_Stopping_Power_MeV_cm2_g_single( stopping_power_source_no, E_MeV_u, particle_no, material_no);

  const double max_electron_range_m   =  AT_max_electron_range_m( E_MeV_u, (int)material_no, (int)er_model);

  const double precalculated_constant_Gy  =  AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model);

  const double d_min_Gy               =  AT_RDD_d_min_Gy( E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, precalculated_constant_Gy );

  // Get density // TODO move it down
  double density_g_cm3  =  AT_density_g_cm3_from_material_no( material_no );
  double density_kg_m3  =  density_g_cm3 * 1000.0;

  double Katz_point_coeff_Gy  =  0.0; // TODO move it down
  if( (rdd_model == RDD_KatzSite) || (rdd_model == RDD_CucinottaPoint) || (rdd_model == RDD_KatzExtTarget) || (rdd_model == RDD_CucinottaExtTarget)){
    Katz_point_coeff_Gy     =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);
  }

  long     i;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Test
  if( rdd_model == RDD_Test){
    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]         =  AT_RDD_Test_Gy(r_m[i], 0., max_electron_range_m, precalculated_constant_Gy);
      // Cut-off low doses at low doses not necessary here
    }
  }// end RDD_Test

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzPoint
  if( rdd_model == RDD_KatzPoint){ // RDD formula will be determined by form of ER model
    const double Katz_point_coeff_Gy  =  precalculated_constant_Gy;
    // Compatible ER models
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) || (er_model == ER_ButtsKatz) ){
      const double alpha              =  AT_ER_PowerLaw_alpha( E_MeV_u );
      const double r_min_m            =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
      // Loop over all r_m given
      for (i = 0; i < n; i++){
        D_RDD_Gy[i]     =  AT_RDD_KatzPoint_Gy(r_m[i], r_min_m, max_electron_range_m, er_model, alpha, Katz_point_coeff_Gy);
        if( D_RDD_Gy[i] > 0.0)
          D_RDD_Gy[i]   =  GSL_MAX(D_RDD_Gy[i], d_min_Gy);          // Cut-off low doses, necessary in CPPSC
        // TODO maybe this cutoff can be moved to CPPSC implementation, in a place where it is needed
      } // end for
    } else {
      for (i = 0; i < n; i++){
        D_RDD_Gy[i]     =  0.0;
      } // end for
      return 1;
    }
  }// end RDD_KatzPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzSite
  if( rdd_model == RDD_KatzSite){
    // Compatible ER models
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) || (er_model == ER_ButtsKatz) ){
      const double LET_J_m   =  LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J;     // convert LET_MeV_cm2_g to LET_J_m
      const double dEdx_J_m  =  precalculated_constant_Gy;                            // take dEdx averaged on outer shell from norm_constant
      const double alpha     =  AT_ER_PowerLaw_alpha( E_MeV_u );
      const double a0_m      =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);

      // Loop over all r_m given
      for (i = 0; i < n; i++){
        D_RDD_Gy[i]     =  AT_RDD_KatzSite_Gy(r_m[i], 0.0, max_electron_range_m, a0_m, er_model, alpha, density_kg_m3, LET_J_m, dEdx_J_m, Katz_point_coeff_Gy);
        if( D_RDD_Gy[i] > 0.0)
          D_RDD_Gy[i]   =  GSL_MAX(D_RDD_Gy[i], d_min_Gy);          // Cut-off low doses, necessary in CPPSC
      } // end for
    } else {
      for (i = 0; i < n; i++){
        D_RDD_Gy[i]     =  0.0;
      } // end for
      return 1;
    }
  }// end RDD_KatzSite

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_Geiss
  if( rdd_model == RDD_Geiss){
    const double a0_m =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]     =  AT_RDD_Geiss_Gy(r_m[i], 0., max_electron_range_m, a0_m, precalculated_constant_Gy);
      // Cut-off at low doses at low doses not necessary here
    }
  }// end RDD_Geiss

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaPoint
  if( rdd_model == RDD_CucinottaPoint){
    const double beta =  AT_beta_from_E_single( E_MeV_u );
    const double r_min_m  =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
    // Loop over all r_m given
    for (i = 0; i < n; i++){
      D_RDD_Gy[i]     =  AT_RDD_CucinottaPoint_Gy(r_m[i], r_min_m, max_electron_range_m, beta, precalculated_constant_Gy, Katz_point_coeff_Gy);
      if( D_RDD_Gy[i] > 0.0)
        D_RDD_Gy[i]   =  GSL_MAX(D_RDD_Gy[i], d_min_Gy);          // Cut-off low doses, necessary in CPPSC
    }
  }// end RDD_CucinottaPoint

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_KatzExtTarget
  if( rdd_model == RDD_KatzExtTarget){
    // Compatible ER models
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) || (er_model == ER_ButtsKatz) ){
      const double Katz_plateau_Gy     =  precalculated_constant_Gy;
      const double Katz_point_r_min_m  =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);
      const double alpha               =  AT_ER_PowerLaw_alpha( E_MeV_u );
      const double a0_m                =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
      for (i = 0; i < n; i++){
        D_RDD_Gy[i]     =  AT_RDD_ExtendedTarget_KatzPoint_Gy(r_m[i], a0_m, er_model, Katz_point_r_min_m, max_electron_range_m, alpha, Katz_plateau_Gy, Katz_point_coeff_Gy);
        if( D_RDD_Gy[i] > 0.0)
          D_RDD_Gy[i]   =  GSL_MAX(D_RDD_Gy[i], d_min_Gy);          // Cut-off low doses, necessary in CPPSC
      }
    } else {
      for (i = 0; i < n; i++){
        D_RDD_Gy[i]     =  0.0;
      } // end for
      return 1;
    }
  }// end RDD_KatzExtTarget

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RDD_CucinottaExtTarget
  if( rdd_model == RDD_CucinottaExtTarget){
    const double a0_m                  =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    const double beta                  =  AT_beta_from_E_single( E_MeV_u );
    const double C_norm                =  precalculated_constant_Gy;
    const double Katz_point_r_min_m    =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);

    const double r_max_m               =  GSL_MIN(a0_m, max_electron_range_m);
    double Cucinotta_plateau_Gy        =  AT_RDD_Cucinotta_Ddelta_average_Gy( Katz_point_r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
    Cucinotta_plateau_Gy              +=  C_norm * AT_RDD_Cucinotta_Dexc_average_Gy( Katz_point_r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);

    for (i = 0; i < n; i++){
      D_RDD_Gy[i]     =  AT_RDD_ExtendedTarget_CucinottaPoint_Gy( r_m[i], a0_m, Katz_point_r_min_m, max_electron_range_m, beta, Katz_point_coeff_Gy, C_norm, Cucinotta_plateau_Gy);
      if( D_RDD_Gy[i] > 0.0)
        D_RDD_Gy[i]   =  GSL_MAX(D_RDD_Gy[i], d_min_Gy);          // Cut-off low doses, necessary in CPPSC
    }
  } // end RDD_CucinottaExtTarget

  return 0;
}


int AT_r_RDD_m  ( const long  n,
    const double  D_RDD_Gy[],
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no,
    double        r_RDD_m[])
{
  long     i;

  const double LET_MeV_cm2_g              =  AT_Stopping_Power_MeV_cm2_g_single( stopping_power_source_no, E_MeV_u, particle_no, material_no);

  const double max_electron_range_m       =  AT_max_electron_range_m(E_MeV_u, (int)material_no, (int)er_model);

  const double r_min_m                    =  AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameter);

  const double precalculated_constant_Gy  =  AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model);

  const double d_min_Gy                   =  AT_RDD_d_min_Gy(E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, precalculated_constant_Gy );

  const double d_max_Gy                   =  AT_RDD_d_max_Gy(E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, stopping_power_source_no);

  if( rdd_model == RDD_Test){
    // Loop over all doses given
    for (i = 0; i < n; i++){
      r_RDD_m[i]    =  AT_inverse_RDD_Test_m( D_RDD_Gy[i], max_electron_range_m);
    }
  }// end RDD_Test

  if( rdd_model == RDD_KatzPoint){
    // Compatible ER models
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) || (er_model == ER_ButtsKatz) ){
      // Loop over all doses given
      const double Katz_point_coeff_Gy  =  precalculated_constant_Gy;
      const double alpha                =  AT_ER_PowerLaw_alpha( E_MeV_u );
      for (i = 0; i < n; i++){
        if( (D_RDD_Gy[i] >= d_min_Gy) && (D_RDD_Gy[i] <= d_max_Gy) ){
          r_RDD_m[i]    =  AT_inverse_RDD_KatzPoint_m( D_RDD_Gy[i], r_min_m, max_electron_range_m, er_model, alpha, Katz_point_coeff_Gy);
        } else {
          r_RDD_m[i]    =  0.0;
        }
      }
    } else {
      for (i = 0; i < n; i++){
          r_RDD_m[i]    =  0.0;
      }
      return 1;
    }
  }// end RDD_KatzPoint

  if( rdd_model == RDD_Geiss){
    const double a0_m    =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);
    // Loop over all doses given
    for (i = 0; i < n; i++){
      r_RDD_m[i]    =  AT_inverse_RDD_Geiss_m( D_RDD_Gy[i], d_min_Gy, d_max_Gy, a0_m, precalculated_constant_Gy);
    }
  }// end RDD_Geiss

  if( rdd_model == RDD_KatzSite){
    // Compatible ER models
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) || (er_model == ER_ButtsKatz) ){
      // TODO it is not good to use r = 0.0 for doses outside the range [dmin,dmax]
      // in the case of RDD_KatzSite RDD is well defined for r = 0.0 : D(r) = dmax
      const double Katz_point_coeff_Gy  =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);
      const double alpha                =  AT_ER_PowerLaw_alpha( E_MeV_u );
      const double a0_m                 =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);

      // Loop over all doses given
      for (i = 0; i < n; i++){
        if( (D_RDD_Gy[i] >= d_min_Gy) && (D_RDD_Gy[i] <= d_max_Gy) ){
          r_RDD_m[i]    =  AT_inverse_RDD_KatzSite_m( D_RDD_Gy[i], r_min_m, max_electron_range_m, a0_m, er_model, alpha, d_max_Gy, Katz_point_coeff_Gy);
        } else {
          r_RDD_m[i]    =  0.0;
        }
      }
    } else {
      for (i = 0; i < n; i++){
        r_RDD_m[i]    =  0.0;
      }
      return 1;
    }
  }// end RDD_KatzSite

  if( rdd_model == RDD_CucinottaPoint){
    // Loop over all doses given
    const double Katz_point_coeff_Gy  =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);
    const double C_norm  =  precalculated_constant_Gy;
    const double beta    =  AT_beta_from_E_single( E_MeV_u );

    for (i = 0; i < n; i++){
      if( (D_RDD_Gy[i] >= d_min_Gy) && (D_RDD_Gy[i] <= d_max_Gy) ){
        r_RDD_m[i]    =  AT_inverse_RDD_Cucinotta_m( D_RDD_Gy[i], r_min_m, max_electron_range_m, er_model, beta, C_norm, Katz_point_coeff_Gy);
      } else {
        r_RDD_m[i]    =  0.0;
      }
    }
  }// end RDD_CucinottaPoint

  if( rdd_model == RDD_KatzExtTarget){
    // Compatible ER models
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) || (er_model == ER_ButtsKatz) ){
      // Loop over all doses given
      const double Katz_point_coeff_Gy  =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);
      const double Katz_plateau_Gy  =  precalculated_constant_Gy;
      const double alpha            =  AT_ER_PowerLaw_alpha( E_MeV_u );
      const double a0_m             =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);

      for (i = 0; i < n; i++){
        if( (D_RDD_Gy[i] >= d_min_Gy) && (D_RDD_Gy[i] <= d_max_Gy) ){
          r_RDD_m[i]    =  AT_inverse_RDD_ExtendedTarget_KatzPoint_m( D_RDD_Gy[i], r_min_m, max_electron_range_m, a0_m, er_model, alpha, Katz_plateau_Gy, Katz_point_coeff_Gy);
        } else {
          r_RDD_m[i]    =  0.0;
        }
      }
    } else {
      for (i = 0; i < n; i++){
        r_RDD_m[i]    =  0.0;
      }
      return 1;
    }
  }// end RDD_KatzExtTarget

  if( rdd_model == RDD_CucinottaExtTarget){
    // Loop over all doses given
    const double Katz_point_coeff_Gy   =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, material_no, er_model);
    const double beta    =  AT_beta_from_E_single( E_MeV_u );
    const double C_norm  =  precalculated_constant_Gy;
    const double a0_m    =  AT_RDD_a0_m(max_electron_range_m, rdd_model, rdd_parameter);

    const double r_max_m               =  GSL_MIN(a0_m, max_electron_range_m);
    double Cucinotta_plateau_Gy        =  AT_RDD_Cucinotta_Ddelta_average_Gy( r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);
    Cucinotta_plateau_Gy              +=  C_norm * AT_RDD_Cucinotta_Dexc_average_Gy( r_min_m, r_max_m, max_electron_range_m, beta, Katz_point_coeff_Gy);

    // TODO around d_max it is hard to find using ziddler algorithm inverse RDD =>
    // calculation in that case should be done manually
    for (i = 0; i < n; i++){
      if( (D_RDD_Gy[i] >= d_min_Gy) && (D_RDD_Gy[i] < (1.0-1e-5)*d_max_Gy) ){
        r_RDD_m[i]    =  AT_inverse_RDD_ExtendedTarget_CucinottaPoint_m( D_RDD_Gy[i], a0_m, r_min_m, max_electron_range_m, beta, Katz_point_coeff_Gy, C_norm, Cucinotta_plateau_Gy);
      } else {
        r_RDD_m[i]    =  0.0;
      }
    }
  }// end RDD_CucinottaExtTarget

  return 0;
}
