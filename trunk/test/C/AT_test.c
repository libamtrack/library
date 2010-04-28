/**
 * @file
 * @brief Dummy file to enable debugging, to be changed by the user anyway they like.
 */

/*
 *    AT_test.c
 *    ===================
 *
 *    Created on: 2009-06-08
 *    Author: grzanka
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

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "AmTrack.h"
#include "AT_Wrapper_R.h"

void test_AT_GSM(  double E_MeV_u[])
{

  long          n                       = 1;
  long          particle_no[]           = {PARTICLE_PROTON_NUMBER}; //proton
  long          material_no             = 1;    //alanine

  double        fluence_cm2[]           = {-3.0};
  long          RDD_model               = 3;    //RDD_site
  double        RDD_parameters[]        = {1e-9 ,1e-10, 0.0};
  long          ER_model                = 4;  //ButtsKatz
  long          gamma_model             = 1; // single hit, single target
  double        gamma_parameters[]      = {1, 10.5e-2, 1 ,1, 0};
  long          N2                      = 20;
  long          N_runs                  = 1;
  long          nX                      = 100;
  double        voxel_size_m            = 1e-9;
  double        fluence_factor          = 1.0;
  bool          write_output            = false;
  bool          lethal_events_mode      = false;

  double        results[10];

//      const long*  N_runs,
//      const double*  voxel_size_m,

  AT_run_GSM_method(n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      RDD_model,
      RDD_parameters,
      ER_model,
      gamma_model,
      gamma_parameters,
      N_runs,
      N2,
      fluence_factor,
      write_output,
      nX,
      voxel_size_m,
      lethal_events_mode,
      results);
  printf("GSM: %.3e %.3e\n", E_MeV_u[0],results[0]);
  printf("Dose = %g, dose check = %g\n", fluence_cm2[0], results[1]);
}

void test_AT_SPIFF(  double E_MeV_u[])
{

  long          n                       = 1;
  long          particle_no[]           = {PARTICLE_PROTON_NUMBER}; //proton
  long          material_no             = 1;    //alanine

  double        fluence_cm2[]           = {-10.0};
  long          RDD_model               = RDD_KatzSite;    //RDD_site
  double        RDD_parameters[]        = {5e-8 ,1e-10, 0};
  long          ER_model                = ER_Geiss;  //ButtsKatz
  long          gamma_model             = 2; // single hit, single target
  double        gamma_parameters[]      = {1, 10.0, 1 ,1, 0};
  long          N2                      = 30;
  double        fluence_factor          = 1.0;
  bool          write_output            = false;
  bool          shrink_tails            = true;
  double        shrink_tails_under      = 1e-40;
  bool          adjust_N2               = true;
  bool          lethal_events_mode      = false;

  double        results[10];

  AT_run_SPIFF_method(n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      RDD_model,
      RDD_parameters,
      ER_model,
      gamma_model,
      gamma_parameters,
      N2,
      fluence_factor,
      write_output,
      shrink_tails,
      shrink_tails_under,
      adjust_N2,
      lethal_events_mode,
      results);
  printf("SPIFF: %.3e %.3e\n", E_MeV_u[0],results[0]);
  printf("Dose = %g, dose check = %g\n", fluence_cm2[0], results[1]);
}


void test_AT_IGK(double E_MeV_u[])
{
        long    n =1 ;
//      double  E_MeV_u[] ={100.0};
        long    particle_no[] = {PARTICLE_PROTON_NUMBER};
        double  fluence_cm2[] = {-10.0};
        long    material_no = 1;
        long    RDD_model = 3;
        double  RDD_parameters[] = {5e-8,1e-10, 0.0};
        long    ER_model = 3;
        long    gamma_model = 2;
        double  gamma_parameters[] = {1.0 , 1.0,1 ,1};
        double  saturation_cross_section_factor = 1.0;
        double  results[10];
 AT_run_IGK_method(
        n,
        E_MeV_u,
        particle_no,
        fluence_cm2,
        material_no,
        RDD_model,
        RDD_parameters,
        ER_model,
        gamma_model,
        gamma_parameters,
        saturation_cross_section_factor,
        results);

 printf("IGK: %.3e %.3e %.3e %.3e\n", E_MeV_u[0],results[0], results[2], results[3]);

}

int main(){
//  long test_pn[] = {1001, 2004, 6012, 8016, 92238};
//  double test_E_MeV_u[] = {100,100,100,100,100};
//  long material_no = 1;
//  long test_A[] = {0,0,0,0,0};
//  long test_Z[] = {0,0,0,0,0};
//  double test_w[] = {0,0,0,0,0};
//  double test_LET[] = {0,0,0,0,0};
//
//  long n_tmp = 5;
//  AT_A_from_particle_no(       n_tmp,
//    test_pn,
//    test_A);
//  AT_Z_from_particle_no(       n_tmp,
//    test_pn,
//    test_Z);
//  AT_atomic_weight_from_particle_no(       n_tmp,
//    test_pn,
//    test_w);
//
//  AT_LET_MeV_cm2_g(     n_tmp,
//      test_E_MeV_u,
//      test_pn,
//      material_no,
//      test_LET);
//  double E_MeV_u[] = {energy};
//  test_AT_IGK(E_MeV_u);
//  test_AT_SPIFF(E_MeV_u);
//  test_AT_GSM(E_MeV_u);
//  while (energy >= 0.1){
//    E_MeV_u[0] = energy;
//    //test_AT_IGK(E_MeV_u);
//    test_AT_SPIFF(E_MeV_u);
//    energy = energy *0.9;
//    }


//  const int n                   =       1;
//  const float E_MeV_u[]         =       {3};
//  const int particle_no[]       =       {1001};
//  const float fluence_cm2[]     =       {-0.01f};
//  const int material_no         =       1;
//  const int rdd_model           =       3;
//  const float rdd_parameters[]  =       {5e-8f};
//  const int er_model            =       3;
//  const int N2                  =       20;
//  int n_bins_f1                 =       0;
//  float f1_parameters[9];
//
//  AT_SC_get_f1_array_size_R(
//      &n,
//      E_MeV_u,
//      particle_no,
//      &material_no,
//      &rdd_model,
//      rdd_parameters,
//      &er_model,
//      &N2,
//      &n_bins_f1,
//      f1_parameters);
//
//  float norm_fluence[]           =      {0};
//  float dose_contribution_Gy[]   =      {0};
//  float f_parameters[7];
//  float* f1_d_Gy                  =      (float*)malloc(n_bins_f1*sizeof(float));
//  float* f1_dd_Gy                 =      (float*)malloc(n_bins_f1*sizeof(float));
//  float* f1                       =      (float*)malloc(n_bins_f1*sizeof(float));
//
//  AT_SC_get_f1_R(
//      &n,
//      E_MeV_u,
//      particle_no,
//      fluence_cm2,
//      &material_no,
//      &rdd_model,
//      rdd_parameters,
//      &er_model,
//      &N2,
//      &n_bins_f1,
//      f1_parameters,
//      norm_fluence,
//      dose_contribution_Gy,
//      f_parameters,
//      f1_d_Gy,
//      f1_dd_Gy,
//      f1);
//
//  free(f1_d_Gy);
//  free(f1_dd_Gy);
//  free(f1);

//  const long n = 9;
//  const long particle_no[] = {1001,2003,2004,3005,6012,8016,92238,137137,512256};
//  char particle_name[n][PARTICLE_NAME_NCHAR];
//
//  AT_particle_name_from_particle_no( n,
//      particle_no,
//      particle_name);
//
//  long i;
//  for( i = 0 ; i < n; i++){
//    printf("particle_name [ %ld ] = %s \n", particle_no[i], particle_name[i]);
//  }

//  const int number_of_bins = 10;
//  const float d_Gy[] = {1,2,3,4,5,6,7,8,9,10};
//  const float dd_Gy[] = {1,1,1,1,1,1,1,1,1,1};
//  const float f[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
//  const int gamma_model = 2;
//  const float gamma_parameter[] = {1,10,1,1,0};
//  const int lethal_events_mode = 0;
//  float S[] = {0,0,0,0,0,0,0,0,0,0};
//  float f0, S_HCP, S_gamma, efficiency;
//
//   AT_SC_get_gamma_response_R(  &number_of_bins,
//      d_Gy,
//      dd_Gy,
//      f,
//      &f0,
//      &gamma_model,
//      gamma_parameter,
//      &lethal_events_mode,
//      S,
//      &S_HCP,
//      &S_gamma,
//      &efficiency);


//const long n = 3;
//const double E_MeV_u[] = {1,10,100};
//const long particle_no[] = {1001,6012,1001};
//const double fluence_cm2[] = {1e7, 1e6, 1e8};
//const long material_no = 1;
//double f_parameters[7];
//const long er_model = 4;
//
//
//AT_RDD_f_parameters( n,
//    E_MeV_u,
//   particle_no,
//    fluence_cm2,
//    material_no,
//    er_model,
//    f_parameters);

//const int number_of_materials = 3;
//const int material_no[] = {1,2,3};
//float dummy1[] = {0,0,0};
//float dummy2[] = {0,0,0};
//float dummy3[] = {0,0,0};
//float dummy4[] = {0,0,0};
//float dummy5[] = {0,0,0};
//float dummy6[] = {0,0,0};
//float dummy7[] = {0,0,0};
//float dummy8[] = {0,0,0};
//
// AT_get_materials_data_R( &number_of_materials,
//    material_no,
//    dummy1,
//    dummy2,
//    dummy3,
//    dummy4,
//    dummy5,
//    dummy6,
//    dummy7,
//    dummy8);


const int particle_no[] = {1001,1001,1001,6012,6012,6012};
int A[] = {0,0,0,0,0,0};


const int n = 6;
AT_A_from_particle_no_R( &n,
    particle_no,
    A);

return 0;

};
