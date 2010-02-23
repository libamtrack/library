#ifndef AT_RDD_H_
#define AT_RDD_H_

/**
 * @file
 * @brief Radial Dose Distribution models
 */

/*
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
#include "AT_GammaResponse.h"

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
      RDD_ExtTarget        = 5,      /* parameters: 0 - r_min [m] (core diameter), 1 - a0 [m] (target diameter), 2 - D_min [Gy] (cut-off dose) */ //as defined in Edmund et al. , 2007
      RDD_Edmund           = 6,      /* parameters: 0 - a0 [m] (core diameter), 1 - d_min_Gy (lower dose cut-off)  */ // after Edmund et al., 2007, but modified with dose-cut off
      RDD_Cucinotta        = 7       /* parameters: TODO  */
};

#define RDD_DATA_N    7

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
    {  RDD_Test,                     RDD_KatzPoint,                                RDD_Geiss,                         RDD_Site,                                        RDD_ExtTarget,                                                RDD_Edmund,                     RDD_Cucinotta},
    {  0,                            2,                                            1,                                 2,                                               3,                                                            2,                              1},
    {  {"","",""},                   {"r_min_m", "d_min_Gy",""},                   {"a0_m","",""},                    {"a0_m","d_min_Gy",""},                          {"r_min_m","a0_m","D_min_Gy"},                                {"a0_m","d_min_Gy",""},         {"r_min_m","",""}},
    {  {0,0,0},                      {1e-10, 1e-10,0},                             {5e-8,0,0},                        {5e-8,1e-10,0},                                  {1e-10, 5e-8, 1e-10},                                         {5e-8,1e-10,0},                 {5e-11,0,0}},
    {  "Simple step test function",  "Katz' point target RDD [Katz et al., 1972]", "Geiss' RDD [Geiss et al., 1998]", "Site RDD, as defined in [Edmund et al., 2007]", "Katz' extended target, as defined in [Edmund et al., 2007]", "Edmund, as defined in [TODO]", "Cucinotta, as defined in [Cucinotta et al. ]"}
};

/**
* Returns name of the radial dose distribution model from index
*
* @param  RDD_no   radial dose distribution model index
* @param  RDD_name string containing radial dose distribution model name (output)
*/
void getRDDName( const long* RDD_no,
    char* RDD_name);

void getRDDNo( const char* RDD_name,
    long* RDD_no);

void AT_D_RDD_Gy( const long*  n,
    const float*  r_m,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const long*   particle_no,
    /* detector parameters */
    const long*   material_no,
    /* radial dose distribution model */
    const long*   rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const long*   er_model,
    const float*  er_parameter,
    float*        D_RDD_Gy);

void AT_r_RDD_m  ( const long*  n,
    const float*  D_RDD_Gy,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const long*   particle_no,
    /* detector parameters */
    const long*   material_no,
    /* radial dose distribution model */
    const long*   rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const long*   er_model,
    const float*  er_parameter,
    float*        r_RDD_m);

/**
 * @param E_MeV_u
 * @param particle_no
 * @param material_no
 * @param rdd_model
 * @param rdd_parameter
 * @param er_model
 * @param er_parameter
 * @param f1_parameters (output)\n
 *     0 - LET_MeV_cm2_g \n
 *     1 - r_min_m \n
 *     2 - r_max_m \n
 *     3 - d_min_Gy \n
 *     4 - d_max_Gy \n
 *     5 - k             (norm. constant) \n
 *     6 - single_impact_fluence_cm2 \n
 *     7 - single_impact_dose_Gy \n
 *     8 - dEdx_MeV_cm2_g
 */
void AT_RDD_f1_parameters(  /* radiation field parameters */
    const float* E_MeV_u,
    const long*  particle_no,
    /* detector parameters */
    const long*  material_no,
    /* radial dose distribution model */
    const long*  rdd_model,
    const float* rdd_parameter,
    /* electron range model */
    const long*  er_model,
    const float* er_parameter,
    /* calculated parameters */
    float * f1_parameters);

//TODO rewrite Katz functions with const parameters

float          AT_D_RDD_Gy_solver(          const float r , void * params );
double         AT_P_RDD(                    double  r_m, void* params);
double         AT_sI_int(                   double  r_m, void* params);
double         AT_D_RDD_Gy_int(             double  r_m, void* params);

inline float   AT_RDD_Katz_point_kernel(    const float* x, const float* alpha);
void           AT_RDD_Katz_point_kernelS(   int *n, float* x, float* alpha, float* f);
inline float   AT_RDD_Katz_point_coeff_Gy(  const float* C_J_m,const float* Z_eff, const float* beta, const float* alpha, const float* density_kg_m3, const float* r_max_m);
inline float   AT_RDD_Katz_point_Gy(        const float* r_m, const float* alpha, const float* r_max_m, const float* Katz_point_coeff_Gy);
void           AT_RDD_Katz_point_GyS(       int *n, float* r_m, float* alpha, float* r_max_m,float* Katz_point_coeff_Gy, float * D);

inline float   AT_RDD_Katz_dEdx_kernel(     const float* x, const float* alpha);
void           AT_RDD_Katz_dEdx_kernelS(    int *n, float* x, float* alpha, float* f);
double         AT_RDD_Katz_dEdx_integrand(  double x, void * params);
inline float   AT_RDD_Katz_dEdx_coeff_J_m(  const float* r_max_m, const float* density_kg_m3, const float* Katz_point_coeff_Gy);
float          AT_RDD_Katz_dEdx_J_m(        const float* alpha, const float* r_min_m, const float* r_max_m, const float* Katz_dEdx_coeff_J_m);

inline float   AT_RDD_Katz_site_Gy(         const float* r_m, const float* alpha, const float* r_min_m, const float* r_max_m, const float* LET_J_m, const float* density_kg_m3, const float* Katz_dEdx_J_m, const float* Katz_point_coeff_Gy);
void           AT_RDD_Katz_site_GyS(        int *n, float* r_m, float* alpha, float* r_min_m, float* r_max_m, float* LET_J_m, float* density_kg_m3, float* Katz_dEdx_J_m, float* Katz_point_coeff_Gy, float * D_Gy);

float          geometryFunctionPhi(         const float* r0_m, const float* a0_m, const float* r_m);

inline float   AT_RDD_Katz_ext_kernel_Gy(   const float* t_m, const float *r_m, const float* a0_m, const float* alpha, const float* r_min_m, const float* r_max_m, const float* Katz_point_coeff_Gy);
void           AT_RDD_Katz_ext_kernel_GyS(  int *n, float* t_m, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy);
double         AT_RDD_Katz_ext_integrand_Gy(double t_m, void * params);
inline float   AT_RDD_Katz_ext_Gy(          const float *r_m, const float* a0_m, const float* alpha, const float* r_min_m, const float* r_max_m, const float* Katz_point_coeff_Gy);
void           AT_RDD_Katz_ext_GyS(         int *n, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy);

typedef struct {
  /* radiation field parameters */
  float*  E_MeV_u;
  long*   particle_no;
  /* detector parameters */
  long*   material_no;
  /* radial dose distribution model */
  long*   rdd_model;
  float*  rdd_parameters;
  /* electron range model */
  long*   er_model;
  float*  er_parameters;
  /* gamma response parameter*/
  float  gamma_parameters[5];
} AT_P_RDD_parameters;

typedef struct {
  long*  n;
  float*  r_m;
  /* radiation field parameters */
  float*  E_MeV_u;
  long*   particle_no;
  /* detector parameters */
  long*   material_no;
  /* radial dose distribution model */
  long*   rdd_model;
  float*  rdd_parameter;
  /* electron range model */
  long*   er_model;
  float*  er_parameter;
  /* calculated parameters */
  float*  D_RDD_Gy;
  float   D0;
} AT_D_RDD_Gy_parameters;

#endif // AT_RDD_H_
