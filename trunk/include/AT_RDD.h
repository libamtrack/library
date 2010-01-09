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

#include <string.h>

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

#include "AT_Constants.h"
#include "AT_DataLET.h"
#include "AT_PhysicsRoutines.h"


extern int indent_counter;
extern char isp[];
extern FILE * debf;


///////////////////////////////////////////////////////////////////////
// RDD DATA

enum RDDModels{
  RDD_Test                 = 1,      /* no parameters */
      RDD_KatzPoint        = 2,      /* parameters: 0 - r_min [m] (lower integration limit), 1 - d_min_Gy (lower dose cut-off) */
      RDD_Geiss            = 3,      /* parameters: 0 - a0 [m] (core diameter) */
      RDD_Site             = 4,      /* parameters: 0 - a0 [m] (core diameter), 1 - d_min_Gy (lower dose cut-off)  */ // after Edmund et al., 2007, but modified with dose-cut off
      RDD_ExtTarget        = 5       /* parameters: 0 - r_min [m] (core diameter), 1 - a0 [m] (target diameter), 2 - D_min [Gy] (cut-off dose) */ //as defined in Edmund et al. , 2007
};

#define RDD_DATA_N    5

typedef struct {
  long    n;
  long    RDD_no[RDD_DATA_N];
  long    n_parameters[RDD_DATA_N];
  char*   parameter_name[RDD_DATA_N][3];
  float   parameter_default[RDD_DATA_N][3];
  char*   RDD_name[RDD_DATA_N];
} rdd_data;

static const rdd_data AT_RDD_Data = {
    RDD_DATA_N,
    {  RDD_Test,          RDD_KatzPoint,          RDD_Geiss,        RDD_Site, RDD_ExtTarget},
    {  0, 2, 1, 2, 3},
    {  {"","",""},{"r_min_m", "d_min_Gy",""},{"a0_m","",""},{"a0_m","d_min_Gy",""},{"r_min_m","a0_m","D_min_Gy"}},
    {  {0,0,0}, {1e-10, 1e-10,0}, {5e-8,0,0}, {5e-8,1e-10,0}, {1e-10, 5e-8, 1e-10}},
    {  "Simple step test function",  "Katz' point target RDD [Katz et al., 1972]",  "Geiss' RDD [Geiss et al., 1998]",    "Site RDD, as defined in [Edmund et al., 2007]", "Katz' extended target, as defined in [Edmund et al., 2007]"}
};

void   getRDDName(  long* RDD_no, char* RDD_name);
void   getRDDNo(char* RDD_name, long* RDD_no);

void AT_D_RDD_Gy(      long*  n,
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

void AT_r_RDD_m  (      long*  n,
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

void AT_RDD_f1_parameters(  /* radiation field parameters */
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

float          AT_D_RDD_Gy_solver(          float r , void * params );

inline float   AT_RDD_Katz_point_kernel(    float* x, float* alpha);
void           AT_RDD_Katz_point_kernelS(   int *n, float* x, float* alpha, float* f);
inline float   AT_RDD_Katz_point_coeff_Gy(  float* C_J_m,float* Z_eff, float* beta, float* alpha, float* density_kg_m3, float* r_max_m);
inline float   AT_RDD_Katz_point_Gy(        float* r_m, float* alpha, float* r_max_m, float* Katz_point_coeff_Gy);
void           AT_RDD_Katz_point_GyS(       int *n, float* r_m, float* alpha, float* r_max_m,float* Katz_point_coeff_Gy, float * D);

inline float   AT_RDD_Katz_dEdx_kernel(     float* x, float* alpha);
void           AT_RDD_Katz_dEdx_kernelS(    int *n, float* x, float* alpha, float* f);
double         AT_RDD_Katz_dEdx_integrand(  double x, void * params);
inline float   AT_RDD_Katz_dEdx_coeff_J_m(  float* r_max_m, float* density_kg_m3, float* Katz_point_coeff_Gy);
float          AT_RDD_Katz_dEdx_J_m(        float* alpha, float* r_min_m, float* r_max_m, float* Katz_dEdx_coeff_J_m);

inline float   AT_RDD_Katz_site_Gy(         float* r_m, float* alpha, float* r_min_m, float* r_max_m, float* LET_J_m, float* density_kg_m3, float* Katz_dEdx_J_m, float* Katz_point_coeff_Gy);
void           AT_RDD_Katz_site_GyS(        int *n, float* r_m, float* alpha, float* r_min_m, float* r_max_m, float* LET_J_m, float* density_kg_m3, float* Katz_dEdx_J_m, float* Katz_point_coeff_Gy, float * D_Gy);

float          geometryFunctionPhi(         float* r0_m, float* a0_m, float* r_m);

inline float   AT_RDD_Katz_ext_kernel_Gy(   float* t_m, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy);
void           AT_RDD_Katz_ext_kernel_GyS(  int *n, float* t_m, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy);
double         AT_RDD_Katz_ext_integrand_Gy(double t_m, void * params);
inline float   AT_RDD_Katz_ext_Gy(          float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy);
void           AT_RDD_Katz_ext_GyS(         int *n, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy);


#endif // AT_RDD_H_
