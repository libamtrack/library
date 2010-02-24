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


/**
 * RDD code numbers
 */
enum RDDModels {
  RDD_Test                 = 1,      /**< no parameters */
      RDD_KatzPoint        = 2,      /**< parameters: 0 - r_min [m] (lower integration limit), 1 - d_min_Gy (lower dose cut-off) */
      RDD_Geiss            = 3,      /**< parameters: 0 - a0 [m] (core diameter) */
      RDD_Site             = 4,      /**< parameters: 0 - a0 [m] (core diameter), 1 - d_min_Gy (lower dose cut-off) \n after Edmund et al., 2007, but modified with dose-cut off  */
      RDD_KatzExtTarget    = 5,      /**< parameters: 0 - r_min [m] (core diameter), 1 - a0 [m] (target diameter), 2 - D_min [Gy] (cut-off dose) */
      RDD_Edmund           = 6,      /**< parameters: 0 - a0 [m] (core diameter), 1 - d_min_Gy (lower dose cut-off) \n after Edmund et al., 2007, but modified with dose-cut off */
      RDD_Cucinotta        = 7       /**< parameters: 0 - r_min [m] (lower integration limit)  */
};

#define RDD_DATA_N    7

/**
 * RDD data
 */
typedef struct {
  long    n;                                   /** TODO */
  long    RDD_no[RDD_DATA_N];                  /** TODO */
  long    n_parameters[RDD_DATA_N];            /** TODO */
  char*   parameter_name[RDD_DATA_N][3];       /** TODO */
  float   parameter_default[RDD_DATA_N][3];    /** TODO */
  char*   RDD_name[RDD_DATA_N];                /** TODO */
} rdd_data;

static const rdd_data AT_RDD_Data = {
    RDD_DATA_N,
    {  RDD_Test,                     RDD_KatzPoint,                                      RDD_Geiss,                         RDD_Site,                                        RDD_KatzExtTarget,                                            RDD_Edmund,                      RDD_Cucinotta},
    {  0,                            2,                                                  1,                                 2,                                               3,                                                             2,                               1},
    {  {"","",""},                   {"r_min_m", "d_min_Gy",""},                         {"a0_m","",""},                    {"a0_m","d_min_Gy",""},                          {"r_min_m","d_min_Gy","a_0_m"},                                {"a0_m","d_min_Gy",""},          {"r_min_m","",""}},
    {  {0,0,0},                      {1e-10, 1e-10,0},                                   {5e-8,0,0},                        {5e-8,1e-10,0},                                  {1e-10, 1e-10,5e-8},                                           {5e-8,1e-10,0},                  {5e-11,0,0}},
    {  "Simple step test function",  "Katz' point target RDD",                           "Geiss' RDD [Geiss et al., 1998]", "Site RDD, as defined in [Edmund et al., 2007]", "Katz' extended target RDD", "Edmund, as defined in [TODO]", "Cucinotta, as defined in [Cucinotta et al. 1997]"}
};

/**
 * Returns name of the radial dose distribution model from index
 *
 * @param[in]   RDD_no   radial dose distribution model index
 * @param[out]  RDD_name string containing radial dose distribution model name
 */
void getRDDName( const long* RDD_no,
    char* RDD_name);

/**
 * Returns number of the radial dose distribution model from its name
 *
 * @param[in]   RDD_name  string containing radial dose distribution model name
 * @param[out]  RDD_no    radial dose distribution model index
 */
void getRDDNo( const char* RDD_name,
    long* RDD_no);

/**
 * Returns RDD as a function of distance r_m
 *
 * @param[in]   n
 * @param[in]   r_m            distance [m]
 * @param[in]   E_MeV_u
 * @param[in]   particle_no
 * @param[in]   material_no
 * @param[in]   rdd_model
 * @param[in]   rdd_parameter
 * @param[in]   er_model
 * @param[in]   er_parameter
 * @param[out]  D_RDD_Gy       dose [Gy]
 */
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

/**
 * Returns distance as a function of dose
 *
 * @param[in]   n
 * @param[in]   D_RDD_Gy            dose [Gy]
 * @param[in]   E_MeV_u
 * @param[in]   particle_no
 * @param[in]   material_no
 * @param[in]   rdd_model
 * @param[in]   rdd_parameter
 * @param[in]   er_model
 * @param[in]   er_parameter
 * @param[out]  r_RDD_m             distance [m]
 */
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
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @param[in] rdd_model
 * @param[in] rdd_parameter
 * @param[in] er_model
 * @param[in] er_parameter
 * @param[out] f1_parameters
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


/**
 * Calculates C constant given by equation: C = 2 pi N e^4 / ( m c^2 (4 pi eps_0)^2)
 * For water: C = 1.36662e-12 [J/m] = 8.53 [MeV/m]
 *
 * @param[in] electron_density_m3 electron density of given material [1/m^3]
 * @return constant C [J/m]
 */
inline float AT_RDD_Katz_C_J_m( const float electron_density_m3);


/**
 * Calculates coefficient
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] C_J_m                    constant C [J/m]
 * @param[in] Z_eff                    effective ion charge Zeff
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @return coeff [Gy]    calculated coefficient
 */
inline float   AT_RDD_Katz_coeff_Gy(  const float C_J_m,
    const float Z_eff,
    const float beta,
    const float material_density_kg_m3,
    const float r_max_m);

/**
 * Calculates power-law ER kernel of Katz point RDD
 *
 * kernel(r)  =  1/x^2 * 1/alpha * (1 - x)^(1/alpha)              [here: x = r/rmax]
 *
 * @param[in] x                        dimensionless x = r/rmax
 * @param[in] alpha                    parameter of ER model
 * @return kernel    calculated kernel
 */
inline float AT_RDD_Katz_PowerLawER_kernel(    const float x,
    const float alpha);


/**
 * Calculates "old" Katz RDD (derived from linear ER model):
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/r * (1/r - 1/rmax)
 *
 * @param[in] r_m                      distance r [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Z_eff                    effective ion charge Zeff
 * @param[in] C_J_m                    constant C [J/m]
 * @return D(r) [Gy] radial dose distribution at distance r
 */
inline float AT_RDD_Katz_LinearER_versionB_Gy(const float r_m,
    const float r_max_m,
    const float material_density_kg_m3,
    const float beta,
    const float Z_eff,
    const float C_J_m);


/**
 * Calculates "new" Katz RDD (derived from power-law (on wmax) ER model):
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/r^2 * 1/alpha * (1 - r/rmax)^(1/alpha)
 *
 * Version A : using pre-calculated constant in following manner:
 *
 * D(r) = coeff * kernel(r)
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 * kernel(r)  =  1/x^2 * 1/alpha * (1 - x)^(1/alpha)              [here: x = r/rmax]
 *
 * @param[in] r_m                      distance r [m]
 * @param[in] alpha                    parameter of ER model
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] radial dose distribution at distance r
 */
inline float   AT_RDD_Katz_PowerLawER_versionA_Gy(        const float r_m,
    const float alpha,
    const float r_max_m,
    const float Katz_point_coeff_Gy);


/**
 * Calculates "new" Katz RDD (derived from power-law (on wmax) ER model):
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/r^2 * 1/alpha * (1 - r/rmax)^(1/alpha)
 *
 * Version B : Direct calculation
 *
 * @param[in] r_m                      distance r [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Z_eff                    effective ion charge Zeff
 * @param[in] C_J_m                    constant C [J/m]
 * @param[in] alpha                    parameter of ER model
 * @return D(r) [Gy] radial dose distribution at distance r
 */
inline float    AT_RDD_Katz_PowerLawER_versionB_Gy(const float r_m,
    const float r_max_m,
    const float material_density_kg_m3,
    const float beta,
    const float Z_eff,
    const float C_J_m,
    const float alpha);


/**
 * TODO
 */
inline float   AT_RDD_Katz_dEdx_coeff_J_m(  const float r_max_m,
    const float material_density_kg_m3,
    const float Katz_point_coeff_Gy);

/**
 * TODO
 */
float          AT_RDD_Katz_PowerLawER_dEdx_versionA_J_m(        const float alpha,
    const float r_min_m,
    const float r_max_m,
    const float Katz_dEdx_coeff_J_m);

/**
 * TODO
 */
inline float   AT_RDD_Katz_PowerLawER_site_versionA_Gy(         const float r_m,
    const float alpha,
    const float r_min_m,
    const float r_max_m,
    const float LET_J_m,
    const float material_density_kg_m3,
    const float Katz_dEdx_J_m,
    const float Katz_point_coeff_Gy);

/**
 * TODO
 */
float          geometryFunctionPhi(         const float r0_m,
    const float a0_m,
    const float r_m);

/**
 * TODO
 */
inline float   AT_RDD_Katz_ext_kernel_Gy(   const float t_m,
    const float r_m,
    const float a0_m,
    const float alpha,
    const float r_min_m,
    const float r_max_m,
    const float Katz_point_coeff_Gy);

/**
 * TODO
 */
double         AT_RDD_Katz_ext_integrand_Gy(double t_m,
    void * params);

/**
 * TODO
 */
inline float   AT_RDD_Katz_ext_Gy(          const float r_m,
    const float a0_m,
    const float alpha,
    const float r_min_m,
    const float r_max_m,
    const float Katz_point_coeff_Gy);

/**
 * TODO
 */
float          AT_D_RDD_Gy_solver(          const float r ,
    void * params );

/**
 * TODO
 */
double         AT_P_RDD(                    double  r_m,
    void* params);

/**
 * TODO
 */
double         AT_sI_int(                   double  r_m,
    void* params);

/**
 * TODO
 */
double         AT_D_RDD_Gy_int(             double  r_m,
    void* params);


/**
 * TODO
 */
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


/**
 * TODO
 */
typedef struct {
  long*  n;
  float*  r_m;
  /* radiation field parameters */
  float*  E_MeV_u;          /**< energy per nucleon */
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
