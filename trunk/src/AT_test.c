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

void test_AT_GSM(  double E_MeV_u[])
{

  long          n                       = 1;
  long          particle_no[]           = {1001}; //carbon
  long          material_no             = 1;    //alanine

  double        fluence_cm2[]           = {-3.0};
  long          RDD_model               = 4;    //RDD_site
  double        RDD_parameters[]        = {2.5e-8 ,1e-10, 0.0};
  long          ER_model                = 2;  //ButtsKatz
  double        ER_parameters[]         = {0.0};
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
      ER_parameters,
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
  long          particle_no[]           = {1001}; //carbon
  long          material_no             = 1;    //alanine

  double        fluence_cm2[]           = {-3.0};
  long          RDD_model               = 3;    //RDD_site
  double        RDD_parameters[]        = {5e-8 ,1e-10, 0.0};
  long          ER_model                = 3;  //ButtsKatz
  double        ER_parameters[]         = {0.0};
  long          gamma_model             = 2; // single hit, single target
  double        gamma_parameters[]      = {1, 10.0, 1 ,1, 0};
  long          N2                      = 10;
  double        fluence_factor          = 1.0;
  bool          write_output            = false;
  bool          shrink_tails            = true;
  double        shrink_tails_under      = 1e-30;
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


void test_AT_IGK(double E_MeV_u[])
{
        long    n =1 ;
//      double  E_MeV_u[] ={100.0};
        long    particle_no[] = {1001};
        double  fluence_cm2[] = {-10.0};
        long    material_no = 1;
        long    RDD_model = 3;
        double  RDD_parameters[] = {5e-8,1e-10, 0.0};
        long    ER_model = 3;
        double  ER_parameters[] = {0.0};
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
        ER_parameters,
        gamma_model,
        gamma_parameters,
        saturation_cross_section_factor,
        results);

 printf("IGK: %.3e %.3e %.3e %.3e\n", E_MeV_u[0],results[0], results[2], results[3]);

}

//int main(){
////  long test_pn[] = {1001, 2004, 6012, 8016, 92238};
////  double test_E_MeV_u[] = {100,100,100,100,100};
////  long material_no = 1;
////  long test_A[] = {0,0,0,0,0};
////  long test_Z[] = {0,0,0,0,0};
////  double test_w[] = {0,0,0,0,0};
////  double test_LET[] = {0,0,0,0,0};
////
////  long n_tmp = 5;
////  AT_A_from_particle_no(       n_tmp,
////    test_pn,
////    test_A);
////  AT_Z_from_particle_no(       n_tmp,
////    test_pn,
////    test_Z);
////  AT_atomic_weight_from_particle_no(       n_tmp,
////    test_pn,
////    test_w);
////
////  AT_LET_MeV_cm2_g(     n_tmp,
////      test_E_MeV_u,
////      test_pn,
////      material_no,
////      test_LET);
//
//  double energy = 10.0;
//  double E_MeV_u[] = {energy};
////  test_AT_IGK(E_MeV_u);
//  test_AT_SPIFF(E_MeV_u);
////  test_AT_GSM(E_MeV_u);
////  while (energy >= 0.1){
////    E_MeV_u[0] = energy;
////    //test_AT_IGK(E_MeV_u);
////    test_AT_SPIFF(E_MeV_u);
////    energy = energy *0.9;
////    }
//
//  return 0;
//
//};
