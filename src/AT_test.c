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

void test_AT_GSM(  float* E_MeV_u)
{

  long          n                       = 1;
  long          particle_no[]           = {1001}; //carbon
  long          material_no             = 1;    //alanine

  float         fluence_cm2[]           = {-2.0f};
  long          RDD_model               = 4;    //RDD_site
  float         RDD_parameters[]        = {2.5e-8 ,1e-10, 0.0};
  long          ER_model                = 2;  //ButtsKatz
  float         ER_parameters[]         = {0.0f};
  long          gamma_model             = 1; // single hit, single target
  float         gamma_parameters[]      = {1, 10.5e-4, 1 ,1, 0};
  long          N2                      = 40;
  long          N_runs                  = 30;
  long          nX                      = 500;
  float         voxel_size_m            = 1e-9;
  float         fluence_factor          = 1.0f;
  bool          write_output            = false;
  bool          lethal_events_mode      = false;

  float         results[10];

//      const long*  N_runs,
//      const float*  voxel_size_m,

  AT_run_GSM_method(&n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      &material_no,
      &RDD_model,
      RDD_parameters,
      &ER_model,
      ER_parameters,
      &gamma_model,
      gamma_parameters,
      &N_runs,
      &N2,
      &fluence_factor,
      &write_output,
      &nX,
      &voxel_size_m,
      &lethal_events_mode,
      results);
  printf("GSM: %.3e %.3e\n", E_MeV_u[0],results[0]);
  printf("Dose = %g, dose check = %g\n", fluence_cm2[0], results[1]);
}

void test_AT_SPIFF(  float* E_MeV_u)
{

  long          n                       = 1;
  long          particle_no[]           = {1001}; //carbon
  long          material_no             = 1;    //alanine

  float         fluence_cm2[]           = {-3.0f};
  long          RDD_model               = 3;    //RDD_site
  float         RDD_parameters[]        = {5e-8 ,1e-10, 0.0};
  long          ER_model                = 3;  //ButtsKatz
  float         ER_parameters[]         = {0.0f};
  long          gamma_model             = 2; // single hit, single target
  float         gamma_parameters[]      = {1, 10.0, 1 ,1, 0};
  long          N2                      = 10;
  float         fluence_factor          = 1.0f;
  bool          write_output            = false;
  bool          shrink_tails            = true;
  float         shrink_tails_under      = 1e-30;
  bool          adjust_N2               = true;
  bool          lethal_events_mode      = false;

  float         results[10];

  AT_run_SPIFF_method(n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      material_no,
      RDD_model,
      RDD_parameters,
      ER_model,
      ER_parameters,
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


void test_AT_IGK(float* E_MeV_u)
{
        long   n =1 ;
//      float  E_MeV_u[] ={100.0};
        long   particle_no[] = {6012};
        float  fluence_cm2[] ={-1.0};
        long   material_no =5;
        long   RDD_model =3;
        float  RDD_parameters[] ={5e-8,1e-10, 0.0};
        long   ER_model= 3;
        float  ER_parameters[] ={0.0};
        long   gamma_model = 2;
        float  gamma_parameters[] = {1.0 , 10.0,1 ,1};
        float  saturation_cross_section_factor = 1.0;
        float  results[10];
 AT_run_IGK_method(
        &n,
        E_MeV_u,
        particle_no,
        fluence_cm2,
        &material_no,
        &RDD_model,
        RDD_parameters,
        &ER_model,
        ER_parameters,
        &gamma_model,
        gamma_parameters,
        &saturation_cross_section_factor,
        results);

 printf("IGK: %.3e %.3e %.3e %.3e\n", E_MeV_u[0],results[0], results[2], results[3]);

}

int main(){
//  long test_pn[] = {1001, 2004, 6012, 8016, 92238};
//  float test_E_MeV_u[] = {100,100,100,100,100};
//  long material_no = 1;
//  long test_A[] = {0,0,0,0,0};
//  long test_Z[] = {0,0,0,0,0};
//  float test_w[] = {0,0,0,0,0};
//  float test_LET[] = {0,0,0,0,0};
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

  double energy = 40.0;
  float E_MeV_u[] = {300.0};
//  test_AT_GSM(E_MeV_u);
  test_AT_SPIFF(E_MeV_u);
//  while (energy >= 0.1){
//    E_MeV_u[0] = energy;
//    //test_AT_IGK(E_MeV_u);
//    test_AT_SPIFF(E_MeV_u);
//    energy = energy *0.9;
//    }

  return 0;

};
