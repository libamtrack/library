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
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_expint.h>
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
#include "AT_RDD_ExtendedTarget.h"


//////////////////////////////////////////////////////// DATA STRUCTURES ////////////////////////////////////////////////////////

/**
 * RDD code numbers
 */
enum RDDModels {
  RDD_Test                   = 1,      /**< no parameters */
      RDD_KatzPoint          = 2,      /**< parameters: 0 - r_min [m] (lower integration limit), 1 - d_min_Gy (lower dose cut-off) */
      RDD_Geiss              = 3,      /**< parameters: 0 - a0 [m] (core diameter) */
      RDD_KatzSite           = 4,      /**< parameters: 0 - a0 [m] (core diameter), 1 - d_min_Gy (lower dose cut-off) \n after Edmund et al., 2007, but modified with dose-cut off  */
      RDD_Edmund             = 5,      /**< parameters: 0 - a0 [m] (core diameter), 1 - d_min_Gy (lower dose cut-off) \n after Edmund et al., 2007, but modified with dose-cut off */
      RDD_CucinottaPoint     = 6,      /**< parameters: 0 - r_min [m] (lower integration limit),1 - d_min_Gy (lower dose cut-off)   */
      RDD_KatzExtTarget      = 7,      /**< TODO: This is a dummy, replace */
      RDD_CucinottaExtTarget = 8,      /**< TODO: This is a dummy, replace */
};

/**
 * Total number of RDD models
 */
#define RDD_DATA_N    8

/**
 * RDD data
 */
typedef struct {
  long    n;                                   /** total number of RDD models */
  long    RDD_no[RDD_DATA_N];                  /** code number of RDD model*/
  long    n_parameters[RDD_DATA_N];            /** number of model parameters */
  char*   parameter_name[RDD_DATA_N][3];       /** list of names of model parameters*/
  float   parameter_default[RDD_DATA_N][3];    /** default values of model parameters */
  char*   RDD_name[RDD_DATA_N];                /** model names */
} rdd_data;

/**
 * Default model parameters and names
 */
static const rdd_data AT_RDD_Data = {
    RDD_DATA_N,
    {  RDD_Test,                     RDD_KatzPoint,                                      RDD_Geiss,                         RDD_KatzSite,                                    RDD_Edmund,                      RDD_CucinottaPoint,                                 RDD_KatzExtTarget,              RDD_CucinottaExtTarget},
    {  0,                            2,                                                  1,                                 2,                                               2,                               2,                                                  3,                              3},
    {  {"","",""},                   {"r_min_m", "d_min_Gy",""},                         {"a0_m","",""},                    {"a0_m","d_min_Gy",""},                          {"a0_m","d_min_Gy",""},          {"r_min_m","d_min_Gy",""},                          {"r_min_m","a0_m","d_min_Gy"},  {"r_min_m","a0_m","d_min_Gy"}},
    {  {0,0,0},                      {1e-10, 1e-10,0},                                   {5e-8,0,0},                        {5e-8,1e-10,0},                                  {5e-8,1e-10,0},                  {5e-11,1e-10,0},                                    {1e-10,1e-8,1e-10},             {5e-11,1e-8,1e-10}},
    {  "Simple step test function",  "Katz' point target RDD",                           "Geiss' RDD [Geiss et al., 1998]", "Site RDD, as defined in [Edmund et al., 2007]", "Edmund, as defined in [TODO]", "Cucinotta, as defined in [Cucinotta et al. 1997]",  "Katz Extended Target",         "Cucinotta Extended Target"}
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
 * @return      RDD_no    radial dose distribution model index
 */
long getRDDNo( const char* RDD_name );

//////////////////////////////////////////////////////// GENERAL FUNCTIONS ////////////////////////////////////////////////////////

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
    float         D_RDD_Gy[]);

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
    float         r_RDD_m[]);


//////////////////////////////////////////////////////// HELPER FUNCTIONS ////////////////////////////////////////////////////////

/**
 * TODO
 * @param max_electron_range_m
 * @param rdd_model
 * @param rdd_parameter
 * @return
 */
double AT_RDD_r_min_m(
    const double  max_electron_range_m,
    const long    rdd_model,
    const float   rdd_parameter[]);


/**
 * TODO
 * @param max_electron_range_m
 * @param rdd_model
 * @param rdd_parameter
 * @return a0
 */
double AT_RDD_a0_m(
    const double  max_electron_range_m,
    const long    rdd_model,
    const float   rdd_parameter[]);

/**
 * TODO
 * @param max_electron_range_m
 * @param LET_MeV_cm2_g
 * @param E_MeV_u
 * @param particle_no
 * @param material_no
 * @param rdd_model
 * @param rdd_parameter
 * @param er_model
 * @return
 */
double AT_RDD_norm_constant_Gy(
    const double  max_electron_range_m,
    const double  LET_MeV_cm2_g,
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const float   rdd_parameter[],
    const long    er_model);

/**
 * TODO
 * @param norm_constant_Gy
 * @param max_electron_range_m
 * @param rdd_model
 * @param rdd_parameter
 * @return
 */
double AT_RDD_d_min_Gy(
    const double  norm_constant_Gy,
    const double  max_electron_range_m,
    const long    rdd_model,
    const float   rdd_parameter[]);


/**
 * TODO
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @param[in] rdd_model
 * @param[in] rdd_parameter
 * @param[in] er_model
 * @param[in] er_parameter
 * @return    d_max_Gy
 */
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
    const float   er_parameter[]);


/**
 * TODO
 * @param[in]  E_MeV_u
 * @param[in]  particle_no
 * @param[in]  material_no
 * @param[in]  rdd_model
 * @param[in]  rdd_parameter
 * @param[in]  er_model
 * @param[in]  er_parameter
 * @param[out] f1_parameters
 *     0 - LET_MeV_cm2_g \n
 *     1 - r_min_m \n
 *     2 - r_max_m \n
 *     3 - d_min_Gy \n
 *     4 - d_max_Gy \n
 *     5 - normalization constant [Gy] \n
 *     6 - single_impact_fluence_cm2 \n
 *     7 - single_impact_dose_Gy \n
 */
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
    float         f1_parameters[]);

//////////////////////////////////////////////////////// KATZ RDD ////////////////////////////////////////////////////////

/**
 * Calculates C constant given by equation:
 *
 * C = 2 pi N e^4 / ( m c^2 (4 pi eps_0)^2)
 *
 * thus:
 *
 * C = 2 pi N e^4 / ( m c^2 16 pi^2 eps_0^2)
 *
 * in simpler form:
 *
 * C = N e^4 / ( 8 m pi (c * eps_0)^2)
 *
 * For water: C = 1.36662e-12 [J/m] = 8.53 [MeV/m]
 *
 * @param[in] electron_density_m3 electron density of given material [1/m^3]
 * @return constant C [J/m]
 */
inline double AT_RDD_Katz_C_J_m( const double electron_density_m3);


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
 * @return coeff [Gy]                  calculated coefficient
 */
inline double AT_RDD_Katz_coeff_Gy(const double C_J_m,
    const double Z_eff,
    const double beta,
    const double density_kg_m3,
    const double r_max_m);


/**
 * Calculates "old" Katz RDD (derived from linear (on wmax) ER model):
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho *  1/r * (1/r - 1/rmax)
 *
 * or:
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/rmax^2 * rmax/r * (rmax/r - 1.)
 *
 * thus:
 *
 * D(r) = coeff * rmax/r * (rmax/r - 1.)
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r_m                      distance r [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] radial dose distribution at distance r
 */
inline double   AT_RDD_Katz_LinearER_Dpoint_Gy(        const double r_m,
    const double r_max_m,
    const double Katz_point_coeff_Gy);



/**
 * Calculates "new" Katz RDD (derived from power-law (on wmax) ER model):
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/r^2 * 1/alpha * (1 - r/rmax)^(1/alpha)
 *
 * thus:
 *
 * D(r) = coeff * 1/alpha * (rmax/r)^2 * (1 - r/rmax)^(1/alpha)
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r_m                      distance r [m]
 * @param[in] alpha                    parameter of ER model
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] radial dose distribution at distance r
 */
inline double   AT_RDD_Katz_PowerLawER_Dpoint_Gy(        const double r_m,
    const double alpha,
    const double r_max_m,
    const double Katz_point_coeff_Gy);


/**
 * Calculates average dose kernel for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * kernel = 2/(x2^2 - x1^2) * \int_x1^x2 1/x^2 * 1/alpha * (1 - x)^(1/alpha) x dx =
 *        = 2/(x2^2 - x1^2) * \int_x1^x2 1/x * 1/alpha * (1 - x)^(1/alpha) dx
 *
 * now we use the information that:
 *
 * \int 1/x * 1/alpha * (1 - x)^(1/alpha) dx = (1-x)^(1/alpha) ((x-1)/x)^(-1/alpha) _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x)+constant
 *
 * thus:
 *
 * kernel =  2/(x2^2 - x1^2) * (F2 - F1)
 *
 * where:
 *
 * F1 = (1-x1)^(1/alpha) ((x1-1)/x1)^(-1/alpha) _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x1)
 * F2 = (1-x2)^(1/alpha) ((x2-1)/x2)^(-1/alpha) _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x2)
 *
 * here _2F_1 is the special hypergeometric function
 *
 * @param[in] x1                     inner radius x1 (lower integration limit)
 * @param[in] x2                     outer radius x2 (upper integration limit)
 * @param[in] alpha                  parameter of ER model
 * @return kernel                    calculated kernel
 */
double   AT_RDD_Katz_PowerLawER_DaverageKernel(  const double x1,
    const double x2,
    const double alpha);


/**
 * Calculates approximate value of average dose kernel for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * kernel = 2/(x2^2 - x1^2) * \int_x1^x2 1/x^2 * 1/alpha * (1 - x)^(1/alpha) x dx =
 *        = 2/(x2^2 - x1^2) * \int_x1^x2 1/x * 1/alpha * (1 - x)^(1/alpha) dx
 *
 * now we use the series expansion:
 *
 * 1/x * 1/alpha * (1 - x)^(1/alpha) = 1 / (x alpha) - 1 / alpha^2  + (1/alpha  - 1) x / (2 alpha^2) + O(x^2)
 *
 * to calculate the integral:
 *
 * \int 1/x * 1/alpha * (1 - x)^(1/alpha) dx \approx x / alpha^2 ( (x / 4alpha) * (1/alpha - 1) - 1 ) + log(x) + C
 *
 * thus:
 *
 * kernel =  2/(x2^2 - x1^2) * (F2 - F1)
 *
 * where:
 *
 * F1 = x1 / alpha^2 ( (x1 / 4alpha) * (1/alpha - 1) - 1 ) + log(x1)
 * F2 = x2 / alpha^2 ( (x2 / 4alpha) * (1/alpha - 1) - 1 ) + log(x2)
 *
 * @param[in] x1                     inner radius x1 (lower integration limit)
 * @param[in] x2                     outer radius x2 (upper integration limit)
 * @param[in] alpha                  parameter of ER model
 * @return kernel                    calculated kernel
 */
double   AT_RDD_Katz_PowerLawER_DaverageKernel_approx(  const double x1,
    const double x2,
    const double alpha);


/**
 * Calculates average dose for "old" Katz RDD (derived from linear (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr
 *
 * thus:
 *
 * Dav(r1,r2) = coeff / (pi r2^2 - pi r1^2)  *  \int_r1^r2 rmax/r * (rmax/r - 1.) 2 pi r dr
 * Dav(r1,r2) = 2 * coeff * rmax^2/ (r2^2 - r1^2)  *  \int_r1^r2 (1/r - 1/rmax) dr
 * Dav(r1,r2) = 2 * coeff / ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 (1/r - 1/rmax) dr
 *
 * Let us calculate integral:
 *
 * \int_r1^r2 (1/r - 1/rmax) dr = log(r) - r/rmax |_r1^r2 = log(r2/r1) - (r2 - r1)/rmax
 *
 * in other words:
 *
 * Dav(r1,r2) = 2 * coeff * ( log(r2/r1) - (r2 - r1)/rmax ) / ((r2/rmax)^2 - (r1/rmax)^2)
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
inline double   AT_RDD_Katz_LinearER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double r_max_m,
    const double Katz_point_coeff_Gy);


/**
 * Calculates average dose for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr
 *
 * thus:
 *
 * D(r) = coeff * kernel(r)
 *
 * Dav(r1,r2) = coeff * 2 / (r2^2 - r1^2) * \int_r1^r2 kernel(r) r dr
 *
 * substituting x1 = r1/rmax , x2 = r2/rmax we will have:
 *
 * Dav(r1,r2) = coeff * 2 / (x2^2 - x1^2) * \int_x1^x2 kernel(x) x dx
 *
 * in other words:
 *
 * Dav(r1,r2) = coeff * kernel_av( x1, x2 )
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] alpha                    parameter of ER model
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double   AT_RDD_Katz_PowerLawER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double r_max_m,
    const double alpha,
    const double Katz_point_coeff_Gy);



/**
 * Calculates energy delivered to shell between radius a_0 and r_max
 * for "old" Katz RDD (derived from linear (on wmax) ER model).
 *
 * dEdx = rho \int_a0^rmax  D(r) 2 pi r dr =
 *      = rho * (pi rmax^2 - pi a0^2) D_av(a0,rmax)
 *
 * because:
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double   AT_RDD_Katz_LinearER_dEdx_J_m(  const double a0_m,
    const double r_max_m,
    const double material_density_kg_m3,
    const double Katz_point_coeff_Gy);

/**
 * Calculates energy delivered to shell between radius a_0 and r_max
 * for "new" Katz RDD (derived from power-law (on wmax) ER model).
 *
 * dEdx = rho \int_a0^rmax D(r) 2 pi r dr =
 *      = rho * (pi rmax^2 - pi a0^2) D_av(a0,rmax)
 *
 * because:
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] alpha                    parameter of ER model
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double   AT_RDD_Katz_PowerLawER_dEdx_J_m(  const double a0_m,
    const double r_max_m,
    const double material_density_kg_m3,
    const double alpha,
    const double Katz_point_coeff_Gy);

/**
 * Calculates Site RDD, which is LET-normalized
 * for "old" Katz RDD (derived from linear (on wmax) ER model).
 *
 * Dsite(r) = 1 / (rho pi a0^2) * (LET - dEdx)   for r < a0
 * Dsite(r) = D(r)                               for r >= a0
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] dEdx_J_m                 dEdx [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double   AT_RDD_Katz_LinearER_DSite_Gy( const double r_m,
    const double a0_m,
    const double r_max_m,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy);

/**
 * Calculates Site RDD, which is LET-normalized
 * for "new" Katz RDD (derived from power-law (on wmax) ER model).
 *
 * Dsite(r) = 1 / (rho pi a0^2) * (LET - dEdx)   for r < a0
 * Dsite(r) = D(r)                               for r >= a0
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] alpha                    parameter of ER model
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] dEdx_J_m                 dEdx [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double   AT_RDD_Katz_PowerLawER_DSite_Gy( const double r_m,
    const double a0_m,
    const double r_max_m,
    const double material_density_kg_m3,
    const double alpha,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy);

/**
 * TODO
 */
inline float   AT_RDD_Katz_dEdx_coeff_J_m(  const float r_max_m,
    const float material_density_kg_m3,
    const float Katz_point_coeff_Gy);

/**
 * TODO
 */
float          AT_RDD_Katz_PowerLawER_dEdx_directVersion_J_m(        const float alpha,
    const float r_min_m,
    const float r_max_m,
    const float Katz_dEdx_coeff_J_m);


//////////////////////////////////////////////////////// CUCINOTTA RDD ////////////////////////////////////////////////////////

/**
 * Calculates short range modification function fS(r)
 * for Cucinotta RDD
 *
 * fS(r) = 1.0/( r0/r + 0.6 + 1.7 beta + 1.1 beta^2)           [here r0 = 10^(-9) [m]]
 *
 * @param[in] r_m                      distance [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @return fS(r)
 */
inline double   AT_RDD_Cucinotta_f_shortRange( const double r_m,
    const double beta);

/**
 * Calculates long range modification function fL(r)
 * for Cucinotta RDD
 *
 * fL(r) = exp( -(r/(0.37rmax))^2 )
 *
 * @param[in] r_m                      distance [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @return fL(r)
 */
inline double   AT_RDD_Cucinotta_f_longRange( const double r_m,
    const double r_max_m);


/**
 * Calculates radial component D_delta
 * for Cucinotta RDD
 *
 * Ddelta(r) = C z^2 / beta^2 1/rho fS(r) fL(r) /r^2
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * Ddelta(r) = coeff * fS(r) * fL(r) * rmax^2/r^2
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r_m                      distance [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return Ddelta(r) [Gy]
 */
inline double   AT_RDD_Cucinotta_Ddelta_Gy( const double r_m,
    const double r_max_m,
    const double beta,
    const double Katz_point_coeff_Gy);

/**
 * Integrand in Cucinotta Ddelta calculation:
 *
 *  fS(r) * fL(r) /r;
 * @param[in] r_m                      distance [m]
 * @param[in] params                   vector of parameters (table of doubles)
 * @return integrand  fS(r) * fL(r) /r
 */
double AT_RDD_Cucinotta_Ddelta_average_integrand_m(  double r_m,
    void * params);

/**
 * Calculates average dose for Cucinotta delta RDD.
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 Ddelta(r) 2 pi r dr
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * Ddelta(r) = coeff * fS(r) * fL(r) * rmax^2/r^2
 *
 * Thus
 *
 * Dav(r1,r2) = 2 * coeff/ ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 fS(r) * fL(r) * 1/r dr
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double  AT_RDD_Cucinotta_Ddelta_average_Gy(  const double r1_m,
    const double r2_m,
    const double r_max_m,
    const double beta,
    const double Katz_point_coeff_Gy);


/**
 * Calculates average dose for Cucinotta excitation RDD with Cnorm = 1
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 Dexc(r) 2 pi r dr
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * Dexc(r) = Cnorm * coeff * exp( - r / 2d ) * (rmax/r)^2            [where d = (beta/2) * (hbar * c / wr) and wr = 13eV ]
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 * Cnorm      =  1
 *
 * Thus
 *
 * Dav(r1,r2) = 2 coeff/ ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 exp( - r / 2d ) * 1/r dr
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double   AT_RDD_Cucinotta_Dexc_average_Gy(  const double r1_m,
    const double r2_m,
    const double r_max_m,
    const double beta,
    const double Katz_point_coeff_Gy);


/**
 * Calculates normalization constant
 * for Cucinotta RDD
 *
 * We should have:
 *
 * LET = 2 pi rho \int_rmin^rmax D(r) r dr
 *
 * Thus:
 *
 * LET = 2 pi rho \int_rmin^rmax Ddelta(r) r dr + 2 pi rho \int_rmin^rmax Dexc(r) r dr
 *
 * and
 *
 * LET / (2 pi rho) = (pi rmax^2 - pi rmin^2) Ddelta_average(rmin,rmax) + (pi rmax^2 - pi rmin^2) * Cnorm * Dexc_average(rmin,rmax)
 *
 * so
 *
 * LET / (2 pi rho * ((pi rmax^2 - pi rmin^2)) )  = Ddelta_average(rmin,rmax) + Cnorm * Dexc_average(rmin,rmax)
 *
 * finally:
 *
 * Cnorm = (LET / (2 pi rho * ((pi rmax^2 - pi rmin^2)) ) - Ddelta_average(rmin,rmax) ) / Dexc_average(rmin,rmax)
 *
 * @param[in] r_min_m                  minimum radius cut-off distance [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return C norm
 */
double   AT_RDD_Cucinotta_Cnorm( const double r_min_m,
    const double r_max_m,
    const double beta,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double Katz_point_coeff_Gy);

/**
 * Calculates excitation component D_exc
 * for Cucinotta RDD
 *
 * Dexc(r) = C exp( - r / 2d ) / r^2            [where d = (beta/2) * (hbar * c / wr) and wr = 13eV ]
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * Dexc(r) = Cnorm * coeff * exp( - r / 2d ) * (rmax/r)^2
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * Cnorm      =  Cnorm = (LET / (2 pi rho * ((pi rmax^2 - pi rmin^2)) ) - Ddelta_average(rmin,rmax) ) / Dexc_average(rmin,rmax)
 *
 * @param[in] r_m                      distance [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return Dexc(r) [Gy]
 */
double   AT_RDD_Cucinotta_Dexc_Gy( const double r_m,
    const double r_max_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy);


//////////////////////////////////////////////////////// RDD MODELS IMPLEMENTATION, SINGLE PARTICLE ////////////////////////////////////////////////////////

/**
 * TODO
 * @param r_m
 * @param r_min_m
 * @param r_max_m
 * @param norm_constant_Gy
 * @return
 */
inline double  AT_RDD_Test_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double norm_constant_Gy);

/**
 * TODO
 * @param r_m
 * @param r_min_m
 * @param r_max_m
 * @param er_model
 * @param alpha
 * @param Katz_point_coeff_Gy
 * @return
 */
inline double  AT_RDD_KatzPoint_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const long er_model,
    const double alpha,
    const double Katz_point_coeff_Gy);


/**
 * TODO
 * @param r_m
 * @param r_min_m
 * @param r_max_m
 * @param a0_m
 * @param er_model
 * @param alpha
 * @param density_kg_m3
 * @param LET_J_m
 * @param dEdx_J_m
 * @param Katz_point_coeff_Gy
 * @return
 */
inline double  AT_RDD_KatzSite_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double a0_m,
    const long er_model,
    const double alpha,
    const double density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy);


/**
 * TODO
 * @param r_m
 * @param r_min_m
 * @param r_max_m
 * @param a0_m
 * @param er_model
 * @param alpha
 * @param density_kg_m3
 * @param LET_J_m
 * @param dEdx_J_m
 * @param Katz_point_coeff_Gy
 * @return
 */
inline double  AT_RDD_Edmund_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy);


/**
 * TODO
 * @param r_m
 * @param r_min_m
 * @param r_max_m
 * @param a0_m
 * @param norm_constant_Gy
 * @return
 */
inline double  AT_RDD_Geiss_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double a0_m,
    const long   norm_constant_Gy);


/**
 * Calculates Cucinotta point RDD
 *
 * D(r)    = Dexc(r) + Ddelta(r)
 *
 * Ddelta(r) = C z^2 / beta^2 1/rho fS(r) fL(r) /r^2
 *
 * Dexc(r) = C exp( - r / 2d ) / r^2
 *
 * @param[in] r_m                      distance [m]
 * @param[in] r_min_m                  minimum cut-off [m]
 * @param[in] r_max_m                  delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy]
 */
inline double   AT_RDD_CucinottaPoint_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy);

//////////////////////////////////////////////////////// OTHERS ////////////////////////////////////////////////////////

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
  long    n;
  float*  r_m;
  /* radiation field parameters */
  float   E_MeV_u;          /**< energy per nucleon */
  long    particle_no;
  /* detector parameters */
  long    material_no;
  /* radial dose distribution model */
  long    rdd_model;
  float*  rdd_parameter;
  /* electron range model */
  long    er_model;
  float*  er_parameter;
  /* calculated parameters */
  float*  D_RDD_Gy;
  float   D0;
} AT_D_RDD_Gy_parameters;

#endif // AT_RDD_H_
