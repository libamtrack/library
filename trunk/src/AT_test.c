/**
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

//#include "AT_Constants.h"
//#include "AT_RDD.h"
//#include "AT_SuccessiveConvolutions.h"

//void AT_SC_get_f1_array_sizeS(  long*  n,
//    float*  E_MeV_u,
//    long*  particle_no,
//    char**  material_name,
//    float*  parameter,
//    long*  N2,
//    // from here: return values
//    long*  n_bins_f1,
//    bool*  debug);
//
//void  AT_SC_get_f_array_size(  float*  u,
//    float*  fluence_factor,
//    long*  N2,
//    long*  n_bins_f1,
//    float*  f1_d_Gy,
//    float*  f1_dd_Gy,
//    float*  f1,
//    // from here: return values
//    long*  n_bins_f,
//    float*  u_start,
//    long*  n_convolutions);
//
//void AT_SC_get_f1S(  long*  n,
//    float*  E_MeV_u,
//    long*  particle_no,
//    float*  fluence_cm2,
//    char**  material_name,
//    float*  parameter,
//    long*  N2,
//    long*  n_bins_f1,
//    // return values
//    float*  norm_fluence,
//    float*  LET_MeV_cm2_g,
//    float*  r_min_m,
//    float*  r_max_m,
//    float*  d_min_Gy,
//    float*  d_max_Gy,
//    float*  k_Gy,
//    float*  single_impact_fluence_cm2,
//    float*  single_impact_dose_Gy,
//    float*  dose_contribution_Gy,
//    float*  total_fluence_cm2,
//    float*  total_dose_Gy,
//    float*  ave_E_MeV,
//    float*  dw_E_MeV,
//    float*  ave_LET_MeV_cm2_g,
//    float*  dw_LET_MeV_cm2_g,
//    float*  u,
//    float*  f1_d_Gy,
//    float*  f1_dd_Gy,
//    float*  f1,
//    bool*  debug);
//
//void  AT_SC_get_f_start(      float*  u_start,
//    long*  n_bins_f1,
//    long*  N2,
//    float*  f1_d_Gy,
//    float*  f1_dd_Gy,
//    float*  f1,
//    long*  n_bins_f,
//    // from here: return values
//    float*  f_d_Gy,
//    float*  f_dd_Gy,
//    float*  f_start);
//
//void AT_SuccessiveConvolutions(  float*  u,
//    long*  n_bins_f,
//    long*  N2,
//    // input + return values
//    long*  n_bins_f_used,
//    float*  f_d_Gy,
//    float*  f_dd_Gy,
//    float*  f,
//    // return values
//    float*  f0,
//    float*  fdd,
//    float*  dfdd,
//    float*  d,
//    bool*  write_output,
//    bool*  shrink_tails,
//    float*  shrink_tails_under,
//    bool*  adjust_N2);
//
//void AT_SC_Loop(  long*  n,
//    float*  E_MeV_u,
//    long*  particle_no,
//    float*  fluence_cm2,
//    long*  slab_no,
//    char*  material_name,
//    float*  parameter,
//    long*  N2,
//    long*  n_slabs,
//    long*  n_gamma_parameter,
//    long*  gamma_model,
//    float*  gamma_parameter,
//    long*  verbose_level,
//    char*  output_fileName,
//    // return values
//    float*  u,
//    float*  total_d_Gy,
//    float*  d,
//    float*  S_HCP,
//    float*  S_gamma,
//    float*  efficiency,
//    float*  S_HCP_total,
//    float*  S_gamma_total,
//    float*  efficiency_total);
//void testRDD(){
//
//    long   n           = 5;
//    float   r_m[]         = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
//
//    float  E_MeV_u        = 100;
//    long  particle_no      = 1;
//    long   material_no       = 2;
//
//    long  rdd_model      = 1;
//  long  n_rdd_parameter    = 0;
//    float   rdd_parameter[]    = {1e-11f, 1e-11f, 0.0f};
//
//  long  er_model      = 2;
//  long  n_er_parameter    = 0;
//    float   er_parameter    = 0.0f;
//
//    float  D_RDD_Gy[]       = { 0, 0, 0, 0, 0 };
//  float  r_RDD_m_back[]    = { 0, 0, 0, 0, 0 };
//
//  int i;
//
//  printf("begin model %d\n", rdd_model);
//  AT_D_RDD_Gy(   &n,
//          r_m,
//          &E_MeV_u,  &particle_no,     &material_no,
//          &rdd_model, &n_rdd_parameter,   rdd_parameter,
//          &er_model,   &n_er_parameter,   &er_parameter,
//          D_RDD_Gy);
//
//  for( i = 0 ; i < n ; i++){
//    printf("end, D_RRD_Gy[%g] = %g, r_RRD_m_back = %g\n", r_m[i], D_RDD_Gy[i], r_RDD_m_back[i]);}
//
//  rdd_model      = 2;
//  n_rdd_parameter    = 2;
//  printf("begin model %d\n", rdd_model);
//  AT_D_RDD_Gy(   &n,
//          r_m,
//          &E_MeV_u,  &particle_no,     &material_no,
//          &rdd_model, &n_rdd_parameter,   rdd_parameter,
//          &er_model,   &n_er_parameter,   &er_parameter,
//          D_RDD_Gy);
//  for( i = 0 ; i < n ; i++){
//    printf("end, D_RRD_Gy[%g] = %g, r_RRD_m_back = %g\n", r_m[i], D_RDD_Gy[i], r_RDD_m_back[i]);}
//
//  rdd_model      = 3;
//  n_rdd_parameter    = 1;
//  rdd_parameter[0]  = 5e-8;
//  printf("begin model %d\n", rdd_model);
//  AT_D_RDD_Gy(   &n,
//          r_m,
//          &E_MeV_u,  &particle_no,     &material_no,
//          &rdd_model, &n_rdd_parameter,   rdd_parameter,
//          &er_model,   &n_er_parameter,   &er_parameter,
//          D_RDD_Gy);
//  for( i = 0 ; i < n ; i++){
//    printf("end, D_RRD_Gy[%g] = %g, r_RRD_m_back = %g\n", r_m[i], D_RDD_Gy[i], r_RDD_m_back[i]);
//  }
//
//  rdd_model      = 4;
//  n_rdd_parameter    = 2;
//  rdd_parameter[0]  = 5e-8;
//  rdd_parameter[1]  = 1e-11;
//  printf("begin model %d\n", rdd_model);
//  AT_D_RDD_Gy(   &n,
//          r_m,
//          &E_MeV_u,  &particle_no,     &material_no,
//          &rdd_model, &n_rdd_parameter,   rdd_parameter,
//          &er_model,   &n_er_parameter,   &er_parameter,
//          D_RDD_Gy);
//  for( i = 0 ; i < n ; i++){
//    printf("end, D_RRD_Gy[%g] = %g, r_RRD_m_back = %g\n", r_m[i], D_RDD_Gy[i], r_RDD_m_back[i]);
//  }
//
//  rdd_model      = 5;
//  n_rdd_parameter    = 3;
//  rdd_parameter[0]  = 1e-11;
//  rdd_parameter[1]  = 5e-8;
//  rdd_parameter[2]  = 1e-11;
//  printf("begin model %d\n", rdd_model);
//  AT_D_RDD_Gy(   &n,
//          r_m,
//          &E_MeV_u,  &particle_no,     &material_no,
//          &rdd_model, &n_rdd_parameter,   rdd_parameter,
//          &er_model,   &n_er_parameter,   &er_parameter,
//          D_RDD_Gy);
//  for( i = 0 ; i < n ; i++){
//    printf("end, D_RRD_Gy[%g] = %g, r_RRD_m_back = %g\n", r_m[i], D_RDD_Gy[i], r_RDD_m_back[i]);
//  }
//
//  n            = 1;
//  float  p_E_MeV_u[]    = {100}; //{100, 10, 1};
//  long  p_particle_no[]  = {1}; //{1, 1, 1};
//  rdd_model        = 4;
//    float  rdd_parameter2[]= {5e-8f, 1e-11f};
//    n_rdd_parameter      = 2;
//  long  N2        = 10;
//  bool  debug      = false;
//  long  n_bins_f1;
//  float*  f1_parameters  = (float*)calloc(n * 9, sizeof(float));
//
//  AT_SC_get_f1_array_size(  /* radiation field parameters */
//                &n,
//                p_E_MeV_u,
//                p_particle_no,
//                /* detector parameters */
//                &material_no,
//                &rdd_model,
//                &n_rdd_parameter,
//                &rdd_parameter2,
//                /* electron range model */
//                &er_model,
//                &n_er_parameter,
//                &er_parameter,
//                /* algorith parameters*/
//                &N2,
//                &debug,
//                // from here: return values
//                &n_bins_f1,
//                f1_parameters);
//
//  float  p_fluence_cm2[]      =  {1e8}; //{1e8, 1e8, 1e8};
//  float*  norm_fluence      =   (float*)calloc(n, sizeof(float));
//  float*  dose_contribution_Gy  =   (float*)calloc(n, sizeof(float));
//  float*  f_parameters      =   (float*)calloc(7, sizeof(float));
//  float*  f1_d_Gy          =   (float*)calloc(n_bins_f1, sizeof(float));
//  float*  f1_dd_Gy        =   (float*)calloc(n_bins_f1, sizeof(float));
//  float*  f1            =   (float*)calloc(n_bins_f1, sizeof(float));
//
//  AT_SC_get_f1(  /* radiation field parameters */
//              &n,
//              p_E_MeV_u,
//              p_particle_no,
//              p_fluence_cm2,
//              /* detector parameters */
//              &material_no,
//              &rdd_model,
//              &n_rdd_parameter,
//              &rdd_parameter,
//              /* electron range model */
//              &er_model,
//              &n_er_parameter,
//              &er_parameter,
//              /* algorith parameters*/
//              &N2,
//              &n_bins_f1,
//              /* f1 parameters*/
//              f1_parameters,
//              // from here: return values
//              norm_fluence,
//              dose_contribution_Gy,
//              f_parameters,
//              f1_d_Gy,
//              f1_dd_Gy,
//              f1);
//
//  float  fluence_factor = 1.0f;
//  long  n_bins_f;
//  float  u_start;
//  long  n_convolutions;
//  AT_SC_get_f_array_size(  &f_parameters[0],
//      &fluence_factor,
//      &N2,
//      &n_bins_f1,
//      f1_d_Gy,
//      f1_dd_Gy,
//      f1,
//      // from here: return values
//      &n_bins_f,
//      &u_start,
//      &n_convolutions);
//
//  float*  f_d_Gy          =   (float*)calloc(n_bins_f, sizeof(float));
//  float*  f_dd_Gy          =   (float*)calloc(n_bins_f, sizeof(float));
//  float*  f_start          =   (float*)calloc(n_bins_f, sizeof(float));
//
//  AT_SC_get_f_start(  &f_parameters[0],
//      &n_bins_f1,
//      &N2,
//      f1_d_Gy,
//      f1_dd_Gy,
//      f1,
//      &n_bins_f,
//      // from here: return values
//      f_d_Gy,
//      f_dd_Gy,
//      f_start);
//
//  bool   write_output = true;
//  bool  shrink_tails       = true;
//  float  shrink_tails_under    = 1e-30f;
//  bool  adjust_N2        = true;
//  long  n_bins_f_used      = n_bins_f;
//
//  float  f0;
//  float*  fdd            =   (float*)calloc(n_bins_f, sizeof(float));
//  float*  dfdd          =   (float*)calloc(n_bins_f, sizeof(float));
//  float  d;
//
//  AT_SuccessiveConvolutions(  &f_parameters[0],
//      &n_bins_f,
//      &N2,
//      // input + return values
//      &n_bins_f_used,
//      f_d_Gy,
//      f_dd_Gy,
//      f_start,
//      // return values
//      &f0,
//      fdd,
//      dfdd,
//      &d,
//      &write_output,
//      &shrink_tails,
//      &shrink_tails_under,
//      &adjust_N2);
//
//
//  free(f1_parameters);
//  free(norm_fluence);
//  free(dose_contribution_Gy);
//  free(f_parameters);
//  free(f1_d_Gy);
//  free(f1_dd_Gy);
//  free(f1);
//  free(f_d_Gy);
//  free(f_dd_Gy);
//  free(f_start);
//  free(fdd);
//  free(dfdd);
//}
//
////
////void test_AT_SC_1(){
////
////  /****************** STEP 1 **************/
////
////  // INPUT :
////  long n = 1;
////  float E_MeV_u[] = {60.};
////  float fluence_cm2[] = {1e8};
////  long particle_index[] = {1};
////  char material_name[50] = "Water, Liquid";
////  char ** mn = (char**)calloc(1,sizeof(char*));
////  *mn = material_name;
////  float rdd_parameter[] = {1e-8};
////  long N2 = 10;
////  bool debug = false;
////
////  // OUTPUT :
////  long n_bins_f1;
////
////  AT_SC_get_f1_array_sizeS(&n,E_MeV_u,particle_index,mn,rdd_parameter,&N2,\
////      &n_bins_f1,\
////      &debug);
////
////  printf("AT_SC_get_f1_array_sizeS finished:\n");
////  printf(" - n_bins_f1 = %ld\n", n_bins_f1);
////
////  /****************** STEP 2 **************/
////
////  // OUTPUT :
////  float*  norm_fluence = (float*)calloc(n,sizeof(float));
////  float*  LET_MeV_cm2_g = (float*)calloc(n,sizeof(float));
////  float*  r_min_m = (float*)calloc(n,sizeof(float));
////  float*  r_max_m = (float*)calloc(n,sizeof(float));
////  float*  d_min_Gy = (float*)calloc(n,sizeof(float));
////  float*  d_max_Gy = (float*)calloc(n,sizeof(float));
////  float*  k_Gy = (float*)calloc(n,sizeof(float));
////  float*  single_impact_fluence_cm2 = (float*)calloc(n,sizeof(float));
////  float*  single_impact_dose_Gy = (float*)calloc(n,sizeof(float));
////  float*  dose_contribution_Gy = (float*)calloc(n,sizeof(float));
////  float  total_fluence_cm2;
////  float  total_dose_Gy;
////  float  ave_E_MeV;
////  float  dw_E_MeV;
////  float  ave_LET_MeV_cm2_g;
////  float  dw_LET_MeV_cm2_g;
////  float  u;
////  float*  f1_d_Gy = (float*)calloc(n_bins_f1,sizeof(float));
////  float*  f1_dd_Gy = (float*)calloc(n_bins_f1,sizeof(float));
////  float*  f1 = (float*)calloc(n_bins_f1,sizeof(float));
////
////
////  AT_SC_get_f1S(&n,E_MeV_u,particle_index,fluence_cm2,mn,rdd_parameter,&N2,&n_bins_f1,\
////      norm_fluence,LET_MeV_cm2_g,r_min_m,r_max_m,d_min_Gy,d_max_Gy,k_Gy,\
////      single_impact_fluence_cm2, single_impact_dose_Gy,dose_contribution_Gy,\
////      &total_fluence_cm2,&total_dose_Gy,&ave_E_MeV,&dw_E_MeV,&ave_LET_MeV_cm2_g,\
////      &dw_LET_MeV_cm2_g,&u,f1_d_Gy,f1_dd_Gy,f1,\
////      &debug);
////
////  printf("AT_SC_get_f1S finished:\n");
////  printf(" - norm_fluence[0] = %g\n", norm_fluence[0]);
////  printf(" - LET_MeV_cm2_g[0] = %g\n", LET_MeV_cm2_g[0]);
////  printf(" - r_min_m[0] = %g\n", r_min_m[0]);
////  printf(" - r_max_m[0] = %g\n", r_max_m[0]);
////  printf(" - d_min_Gy[0] = %g\n", d_min_Gy[0]);
////  printf(" - d_max_Gy[0] = %g\n", d_max_Gy[0]);
////  printf(" - k_Gy[0] = %g\n", k_Gy[0]);
////  printf(" - single_impact_fluence_cm2[0] = %g\n", single_impact_fluence_cm2[0]);
////  printf(" - single_impact_dose_Gy[0] = %g\n", single_impact_dose_Gy[0]);
////  printf(" - dose_contribution_Gy[0] = %g\n", dose_contribution_Gy[0]);
////  printf(" - total_fluence_cm2 = %g\n", total_fluence_cm2);
////  printf(" - total_dose_Gy = %g\n", total_dose_Gy);
////  printf(" - ave_E_MeV = %g\n", ave_E_MeV);
////  printf(" - dw_E_MeV = %g\n", dw_E_MeV);
////  printf(" - ave_LET_MeV_cm2_g = %g\n", ave_LET_MeV_cm2_g);
////  printf(" - dw_LET_MeV_cm2_g = %g\n", dw_LET_MeV_cm2_g);
////  printf(" - u = %g\n", u);
////  printf(" - f1_d_Gy[] = %g,%g,%g,...\n", f1_d_Gy[0], f1_d_Gy[1], f1_d_Gy[2]);
////  printf(" - f1_dd_Gy[] = %g,%g,%g,...\n", f1_dd_Gy[0], f1_dd_Gy[1], f1_dd_Gy[2]);
////  printf(" - f1[] = %g,%g,%g,...\n", f1[0], f1[1], f1[2]);
////
////  /****************** STEP 3 **************/
////
////  // INPUT :
////  float fluence_factor = 1.0;
////  float u_1 = u;
////  // OUTPUT :
////  long n_bins_f;
////  float u_start;
////  long n_convolutions;
////
////  AT_SC_get_f_array_size(\
////      &u_1,&fluence_factor,&N2,&n_bins_f1,\
////      f1_d_Gy,f1_dd_Gy,f1,\
////      &n_bins_f,&u_start,&n_convolutions);
////
////  printf("AT_SC_get_f_array_size: finished\n");
////  printf(" - n_bins_f = %ld\n", n_bins_f);
////  printf(" - u_start = %g\n", u_start);
////  printf(" - n_convolutions = %ld\n", n_convolutions);
////
////
////  /****************** STEP 4 **************/
////
////  // OUTPUT:
////  float * f_d_Gy = (float*)calloc(n_bins_f,sizeof(float));
////  float * f_dd_Gy = (float*)calloc(n_bins_f,sizeof(float));
////  float * f_start = (float*)calloc(n_bins_f,sizeof(float));
////
////  AT_SC_get_f_start(\
////      &u_start,&n_bins_f1,&N2,f1_d_Gy,f1_dd_Gy,f1,&n_bins_f,\
////      f_d_Gy,f_dd_Gy,f_start);
////
////  printf("AT_SC_get_f_start: finished\n");
////  printf(" - f_d_Gy[] = %g,%g,%g,...\n", f_d_Gy[0], f_d_Gy[1], f_d_Gy[2]);
////  printf(" - f_dd_Gy[] = %g,%g,%g,...\n", f_dd_Gy[0], f_dd_Gy[1], f_dd_Gy[2]);
////  printf(" - f_start[] = %g,%g,%g,...\n", f_start[0], f_start[1], f_start[2]);
////
////  /****************** STEP 5 **************/
////
////  // INPUT:
////  bool write_output = true;
////  bool shrink_tails = false;
////  float shrink_tails_under = 1e-30;
////  bool adjust_N2 = true;
////
////  long n_bins_f_used = n_bins_f1;
////  // OUTPUT:
////  float f0 = 0.0;
////  float * fdd = (float*)calloc(n_bins_f,sizeof(float));
////  float * dfdd = (float*)calloc(n_bins_f,sizeof(float));
////  float d = 0.0;
////
////  AT_SuccessiveConvolutions(&u,&n_bins_f,&N2,&n_bins_f_used,\
////      f_d_Gy,f_dd_Gy,f_start,\
////      &f0,fdd,dfdd,&d,\
////      &write_output,&shrink_tails,&shrink_tails_under,&adjust_N2);
////
////  printf("AT_SuccessiveConvolutions: finished\n");
////  printf(" - f0 = %g\n", f0);
////  printf(" - fdd[] = %g,%g,%g,...\n", fdd[0], fdd[1], fdd[2]);
////  printf(" - dfdd[] = %g,%g,%g,...\n", dfdd[0], dfdd[1], dfdd[2]);
////  printf(" * n_bins_f_used = %ld\n", n_bins_f_used);
////  printf(" * f_d_Gy[] = %g,%g,%g,...\n", f_d_Gy[0], f_d_Gy[1], f_d_Gy[2]);
////  printf(" * f_dd_Gy[] = %g,%g,%g,...\n", f_dd_Gy[0], f_dd_Gy[1], f_dd_Gy[2]);
////  printf(" * f_start[] = %g,%g,%g,...\n", f_start[0], f_start[1], f_start[2]);
////
////  free(norm_fluence);
////  free(LET_MeV_cm2_g);
////  free(r_min_m);
////  free(r_max_m);
////  free(d_min_Gy);
////  free(d_max_Gy);
////  free(k_Gy);
////  free(single_impact_fluence_cm2);
////  free(single_impact_dose_Gy);
////  free(dose_contribution_Gy);
////  free(f1_d_Gy);
////  free(f1_dd_Gy);
////  free(f1);
////  free(f_d_Gy);
////  free(f_dd_Gy);
////  free(f_start);
////  free(fdd);
////  free(dfdd);
////}
////
////void test_AT_SC_2(){
////
////  // INPUT :
////  long n = 1;
////  float E_MeV_u[] = {60.};
////  long particle_no[] = {1};
////  float fluence_cm2[] = {1e8};
////  long slab_no[] = {1};
////  char material_name[50] = "Water, Liquid";
////  //char ** mn = (char**)calloc(1,sizeof(char*));
////  //*mn = material_name;
////  float rdd_parameter[] = {1e-8};
////  long N2 = 1;
////  long n_slabs = 1;
////  long n_gamma_parameter = 2;
////  long gamma_model = 0;
////  float gamma_parameter[] = {1., 1.};
////  long verbose_level = 2;
////  char output_fileName[50] = "SC.dat";
////  float*  u = (float*)calloc(n_slabs,sizeof(float));
////  float*  total_d_Gy = (float*)calloc(n_slabs,sizeof(float));
////  float*  d = (float*)calloc(n_slabs,sizeof(float));
////  float*  S_HCP = (float*)calloc(n_slabs,sizeof(float));
////  float*  S_gamma = (float*)calloc(n_slabs,sizeof(float));
////  float*  efficiency = (float*)calloc(n_slabs,sizeof(float));
////  float*  S_HCP_total = (float*)calloc(n_slabs,sizeof(float));
////  float*  S_gamma_total = (float*)calloc(n_slabs,sizeof(float));
////  float  efficiency_total;
////
////  AT_SC_Loop(&n, E_MeV_u, particle_no, fluence_cm2, slab_no, material_name,\
////      rdd_parameter,&N2,&n_slabs,&n_gamma_parameter,\
////      &gamma_model,gamma_parameter,&verbose_level,output_fileName,\
////      u,total_d_Gy,d,S_HCP,S_gamma,efficiency,S_HCP_total,\
////      S_gamma_total,&efficiency_total);
////
////  free(u);
////  free(total_d_Gy);
////  free(d);
////  free(S_HCP);
////  free(S_gamma);
////  free(efficiency);
////  free(S_HCP_total);
////  free(S_gamma_total);
////}

void test_AT_SPISS(){

  // INPUT :
  long 		n 					= 2;
  float 	E_MeV_u[] 			= {10.,50};
  long 		particle_no[] 		= {1, 1};
  float 	fluence_cm2[] 		= {-0.025f, -0.005f};
  long		material_no 		= 1;			// Water

  long		RDD_model			= 3;			// Gei�
  float 	RDD_parameters[] 	= {5e-8};
  long		ER_model			= 4;			// Gei�
  float		ER_parameters[]		= {0.0f};
  long		GR_model			= 4;			// Exp-sat
  float		GR_parameters[]		= {1, 10};

  long		n_runs				= 1000000;
  long 		N2 					= 40;
  float		fluence_factor		= 1.0f;
  int		write_output		= 1;
  long		importance_sampling	= 0;

  float		results[10];

  AT_efficiency_SPISS(	&n,
						E_MeV_u,
						particle_no,
						fluence_cm2,
						&material_no,
						&RDD_model,
						RDD_parameters,
						&ER_model,
						ER_parameters,
						&GR_model,
						GR_parameters,
						&n_runs,
						&N2,
						&fluence_factor,
						&write_output,
						&importance_sampling,
						results);

}

void test_AT_grid(){

  // INPUT :
  long 		n 					= 1;
  float 	E_MeV_u[] 			= {10.,50};
  long 		particle_no[] 		= {1, 1};
  float 	fluence_cm2[] 		= {-0.25f, -0.005f};
  long		material_no 		= 1;			// Water

  long		RDD_model			= 3;			// Gei�
  float 	RDD_parameters[] 	= {5e-8};
  long		ER_model			= 4;			// Gei�
  float		ER_parameters[]		= {0.0f};
  long		GR_model			= 4;			// Exp-sat
  float		GR_parameters[]		= {1, 10};

  long		n_runs				= 10;
  long 		N2 					= 20;
  float		fluence_factor		= 1.0f;
  bool		write_output		= true;

  long		nX					= 100;
  float		grid_size_m			= 1e-6;
  bool		lethal_events_mode	= false;

  float		results[10];

  AT_efficiency_grid(	&n,
						E_MeV_u,
						particle_no,
						fluence_cm2,
						&material_no,
						&RDD_model,
						RDD_parameters,
						&ER_model,
						ER_parameters,
						&GR_model,
						GR_parameters,
						&n_runs,
						&N2,
						&fluence_factor,
						&write_output,
						&nX,
						&grid_size_m,
						&lethal_events_mode,
						results);
  }

int main(){

  //testRDD();
  //test_AT_SC_1();
  //test_AT_SC_2();

//  test_AT_SPISS();
  test_AT_grid();

  return 0;
};
