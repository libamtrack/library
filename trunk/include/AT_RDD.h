#ifndef AT_RDD_H_
#define AT_RDD_H_

/**
 *    AT_RDD.h
 *    ========
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

#include "AT_Constants.h"
#include "AT_Data.h"
#include "AT_Utils.h"
#include "AT_Functions.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>

extern int indent_counter;
extern char isp[];
extern FILE * debf;

void       AT_D_RDD_Gy(      long*  n,
                    float*  r_m,
                    /* radiation field parameters */
                    float*  E_MeV_u,
                    long*  particle_no,
                    /* detector parameters */
                    long*  material_no,
                    /* radial dose distribution model */
                    long*  rdd_model,
                    float*  rdd_parameter,
                    /* electron range model */
                    long*  er_model,
                    float*  er_parameter,
                    float*  D_RDD_Gy);

void       AT_r_RDD_m  (      long*  n,
                    float*  D_RDD_Gy,
                    /* radiation field parameters */
                    float*  E_MeV_u,
                    long*  particle_no,
                    /* detector parameters */
                    long*  material_no,
                    /* radial dose distribution model */
                    long*  rdd_model,       /* */
                    float*  rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
                    /* electron range model */
                    long*  er_model,
                    float*  er_parameter,
                    float*  r_RDD_m);

void       AT_RDD_f1_parameters(  /* radiation field parameters */
                    float*  E_MeV_u,
                    long*  particle_no,
                    /* detector parameters */
                    long*  material_no,
                    /* radial dose distribution model */
                    long*  rdd_model,
                    float*  rdd_parameter,
                    /* electron range model */
                    long*  er_model,
                    float*  er_parameter,
                    /* calculated parameters */
                    float * f1_parameters);

float       AT_D_RDD_Gy_solver(   float r , void * params );

inline float   AT_RDD_Katz_point_kernel(    float* x, float* alpha);
void       AT_RDD_Katz_point_kernelS(    int *n, float* x, float* alpha, float* f);
inline float   AT_RDD_Katz_point_coeff_Gy(  float* C_J_m,float* Z_eff, float* beta, float* alpha, float* density_kg_m3, float* r_max_m);
inline float   AT_RDD_Katz_point_Gy(      float* r_m, float* alpha, float* r_max_m, float* Katz_point_coeff_Gy);
void       AT_RDD_Katz_point_GyS(      int *n, float* r_m, float* alpha, float* r_max_m,float* Katz_point_coeff_Gy, float * D);

inline float   AT_RDD_Katz_dEdx_kernel(    float* x, float* alpha);
void       AT_RDD_Katz_dEdx_kernelS(    int *n, float* x, float* alpha, float* f);
double       AT_RDD_Katz_dEdx_integrand(  double x, void * params);
inline float   AT_RDD_Katz_dEdx_coeff_J_m(  float* r_max_m, float* density_kg_m3, float* Katz_point_coeff_Gy);
float       AT_RDD_Katz_dEdx_J_m(      float* alpha, float* r_min_m, float* r_max_m, float* Katz_dEdx_coeff_J_m);

inline float   AT_RDD_Katz_site_Gy(      float* r_m, float* alpha, float* r_min_m, float* r_max_m, float* LET_J_m, float* density_kg_m3, float* Katz_dEdx_J_m, float* Katz_point_coeff_Gy);
void       AT_RDD_Katz_site_GyS(      int *n, float* r_m, float* alpha, float* r_min_m, float* r_max_m, float* LET_J_m, float* density_kg_m3, float* Katz_dEdx_J_m, float* Katz_point_coeff_Gy, float * D_Gy);

float      geometryFunctionPhi(      float* r0_m, float* a0_m, float* r_m);

inline float   AT_RDD_Katz_ext_kernel_Gy(    float* t_m, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy);
void       AT_RDD_Katz_ext_kernel_GyS(  int *n, float* t_m, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy);
double       AT_RDD_Katz_ext_integrand_Gy(  double t_m, void * params);
inline float   AT_RDD_Katz_ext_Gy(      float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy);
void       AT_RDD_Katz_ext_GyS(      int *n, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy);



#endif // AT_RDD_H_
