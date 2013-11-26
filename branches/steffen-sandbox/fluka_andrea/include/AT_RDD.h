#ifndef AT_RDD_H_
#define AT_RDD_H_

/**
 * @brief Radial Dose Distribution models
 */

/*
 *    AT_RDD.h
 *    ========
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
#include "AT_DataStoppingPower.h"
#include "AT_PhysicsRoutines.h"
#include "AT_RDD_ExtendedTarget.h"
#include "AT_RDD_Simple.h"
#include "AT_RDD_ShellAveraged.h"

//////////////////////////////////////////////////////// DATA STRUCTURES ////////////////////////////////////////////////////////

/**
 * RDD code numbers
 */
enum RDDModels {
  RDD_Test                   = 1,      /**< no parameters */
      RDD_KatzPoint          = 2,      /**< parameters: 0 - r_min [m] (lower integration limit), 1 - d_min_Gy (lower dose cut-off) */
      RDD_Geiss              = 3,      /**< parameters: 0 - a0 [m] (core diameter/target size) */
      RDD_KatzSite           = 4,      /**< parameters: 0 - a0 [m] (core diameter/target size), 1 - d_min_Gy (lower dose cut-off) \n after Edmund et al., 2007, but modified with dose-cut off  */
      RDD_CucinottaPoint     = 5,      /**< parameters: 0 - r_min [m] (lower integration limit),1 - d_min_Gy (lower dose cut-off)   */
      RDD_KatzExtTarget      = 6,      /**< parameters: 0 - KatzPoint_r_min [m] (lower integration limit in KatzPoint RDD), 1 - a0 [m] (core diameter/target size), 2 - d_min_Gy (lower dose cut-off) */
      RDD_CucinottaExtTarget = 7,      /**< parameters: 0 - CucinottaPoint_r_min [m] (lower integration limit in CucinottaPoint RDD), 1 - a0 [m] (core diameter/target size), 2 - d_min_Gy (lower dose cut-off) */
};


/**
 * Total number of RDD models
 */
#define RDD_DATA_N    7


/**
 * Maximum number of RDD model parameters
 */
#define RDD_MAX_NUMBER_OF_PARAMETERS    3


/**
 * Length of the single particle component characteristics,
 * length of f1_parameters array is AT_SC_F1_PARAMETERS_SINGLE_LENGTH * number of particle components in the field
 */
#define AT_SC_F1_PARAMETERS_SINGLE_LENGTH 8


/**
 * @struct AT_rdd_data_struct
 * RDD data
 */
typedef struct {
  const long    n;                                                              /** total number of RDD models */
  const long    RDD_no[RDD_DATA_N];                                             /** code number of RDD model */
  const int     n_parameters[RDD_DATA_N];                                       /** number of model parameters */
  const char*   parameter_name[RDD_DATA_N][RDD_MAX_NUMBER_OF_PARAMETERS];       /** list of names of model parameters */
  const double  parameter_default[RDD_DATA_N][RDD_MAX_NUMBER_OF_PARAMETERS];    /** default values of model parameters */
  const char*   RDD_name[RDD_DATA_N];                                           /** model names */
} AT_rdd_data_struct;


/**
 * Default model parameters and names
 */
static const AT_rdd_data_struct AT_RDD_Data = {
    RDD_DATA_N,
    {  RDD_Test,                     RDD_KatzPoint,                                      RDD_Geiss,                         RDD_KatzSite,                                    RDD_CucinottaPoint,                                 RDD_KatzExtTarget,                        RDD_CucinottaExtTarget},
    {  0,                            2,                                                  1,                                 2,                                               2,                                                  3,                                        3},
    {  {"","",""},                   {"r_min_m", "d_min_Gy",""},                         {"a0_m","",""},                    {"a0_m","d_min_Gy",""},                          {"r_min_m","d_min_Gy",""},                          {"KatzPoint_r_min_m","a0_m","d_min_Gy"},  {"CucinottaPoint_r_min_m","a0_m","d_min_Gy"}},
    {  {0,0,0},                      {1e-10, 1e-10,0},                                   {5e-8,0,0},                        {5e-8,1e-10,0},                                  {5e-11,1e-10,0},                                    {1e-10,1e-8,1e-10},                       {5e-11,1e-8,1e-10}},
    {  "Simple step test function",  "Katz' point target RDD",                           "Geiss' RDD [Geiss et al., 1998]", "Site RDD, as defined in [Edmund et al., 2007]", "Cucinotta, as defined in [Cucinotta et al. 1997]",  "Katz Extended Target",                   "Cucinotta Extended Target"}
};


/**
 * Get index of RDD in AT_RDD_Data for given RDD_number
 * (currently for example RDD with number 2 has index 1)
 *
 * @param[in] RDD_number  RDD number
 * @return            RDD index in AT_RDD_Data table
 */
long AT_RDD_index_from_RDD_number( const long RDD_number );


/**
 * Returns name of the radial dose distribution model from model number
 *
 * @param[in]   RDD_no   radial dose distribution model number
 * @param[out]  RDD_name string containing radial dose distribution model name
 * @return      status code
 */
int AT_RDD_name_from_number( const long RDD_no,
    char* RDD_name);


/**
 * Returns number of the radial dose distribution model from its name
 *
 * @param[in]   RDD_name  string containing radial dose distribution model name
 * @return      RDD_no    radial dose distribution model index
 */
long AT_RDD_number_from_name( const char* RDD_name );


/**
 * Returns number of parameters of the radial dose distribution model from model number
 *
 * @param[in]   RDD_model   radial dose distribution model number
 * return                   number of RDD parameters
 */
int AT_RDD_number_of_parameters( const long RDD_model );


//////////////////////////////////////////////////////// GENERAL FUNCTIONS ////////////////////////////////////////////////////////

/**
 * Returns local dose as a function of distance r_m for a given radial dose distribution model
 *
 * @param[in]   n              number of particles (length of r_m vector)
 * @param[in]   r_m            distance [m] (array of size n)
 * @param[in]   E_MeV_u        particle (ion) energy per nucleon [MeV/u] (single number, no mixed fields)
 * @param[in]   particle_no    particle code number (single number, no mixed fields)
 * @param[in]   material_no    material code number (single number, no mixed fields)
 * @param[in]   rdd_model      radial dose distribution model index
 * @param[in]   rdd_parameter  radial dose distribution model parameters (array of size 4)
 * @param[in]   er_model       electron range / track with model index
 * @param[in]   stopping_power_source_no  TODO
 * @param[out]  D_RDD_Gy       dose [Gy] (array of size n)
 * @return status code
 */
int AT_D_RDD_Gy( const long  n,
    const double  r_m[],
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no,
    double        D_RDD_Gy[]);


/**
 * Returns distance as a function of dose
 *
 * @param[in]   n                   number of particles (length of D_RDD_Gy vector)
 * @param[in]   D_RDD_Gy            dose [Gy] (array of size n)
 * @param[in]   E_MeV_u             particle (ion) energy per nucleon [MeV/u]
 * @param[in]   particle_no         particle code number
 * @param[in]   material_no         material code number
 * @param[in]   rdd_model           Radial Dose Distribution model code number
 * @param[in]   rdd_parameter       Radial Dose Distribution model parameters vector (array of size 4)
 * @param[in]   er_model            delta electron range model code number
 * @param[in]   stopping_power_source_no   TODO
 * @param[out]  r_RDD_m             distance [m] (array of size n)
 * @return status code
 */
int AT_r_RDD_m( const long  n,
    const double  D_RDD_Gy[],
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no,
    double        r_RDD_m[]);


//////////////////////////////////////////////////////// HELPER FUNCTIONS ////////////////////////////////////////////////////////

/**
 * Returns lower integration limit for given RDD (or 0.0 if not used).
 * @param[in] max_electron_range_m           delta electron maximum range [m]
 * @param[in] rdd_model                      Radial Dose Distribution model code number
 * @param[in] rdd_parameter                  Radial Dose Distribution model parameters vector  (array of size 3)
 * @return rmin - lower integration limit for given RDD [m]
 */
double AT_RDD_r_min_m( const double  max_electron_range_m,
    const long    rdd_model,
    const double  rdd_parameter[]);


/**
 * Returns target size (core size) for given RDD (or 0.0 if not used).
 * @param[in] max_electron_range_m           delta electron maximum range [m]
 * @param[in] rdd_model                      Radial Dose Distribution model code number
 * @param[in] rdd_parameter                  Radial Dose Distribution model parameters vector  (array of size 3)
 * @return a0 - target size (core size) for given RDD [m]
 */
double AT_RDD_a0_m( const double  max_electron_range_m,
    const long    rdd_model,
    const double  rdd_parameter[]);


/**
 * Returns precalculated constant for given RDD: normalization factor or some coefficient
 * @param[in] max_electron_range_m           delta electron maximum range [m]
 * @param[in] LET_MeV_cm2_g
 * @param[in] E_MeV_u                        particle (ion) energy per nucleon [MeV/u]
 * @param[in] particle_no                    particle code number
 * @param[in] material_no                    material code number
 * @param[in] rdd_model                      Radial Dose Distribution model code number
 * @param[in] rdd_parameter                  Radial Dose Distribution model parameters vector  (array of size 3)
 * @param[in] er_model                       delta electron range model code number
 * @return precalculated constant for given RDD [Gy]
 */
double AT_RDD_precalculated_constant_Gy( const double  max_electron_range_m,
    const double  LET_MeV_cm2_g,
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model);


/**
 * Returns minimal dose for given RDD: either value given by user is taken
 * (in case of RDDs which gives zero-dose at maximum range (rmax))
 * or value of RDD is calculated at rmax and it is taken as Dmax (if higher than user
 * provided value).
 * @param[in] E_MeV_u                         particle (ion) energy per nucleon [MeV/u]
 * @param[in] particle_no                     particle code number
 * @param[in] material_no                     material code number
 * @param[in] rdd_model                       Radial Dose Distribution model code number
 * @param[in] rdd_parameter                   Radial Dose Distribution model parameters vector  (array of size 3)
 * @param[in] er_model                        delta electron range model code number
 * @param[in] precalculated_constant_Gy
 * @return minimal dose for given RDD [Gy]
 */
double AT_RDD_d_min_Gy( const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const double  precalculated_constant_Gy);


/**
 * Returns maximum value of the dose that it is possible to obtain with given RDD.
 * It is calculated at r = r_min, as most of the RDDs explode to infinity at r = 0
 * @param[in] E_MeV_u                    particle (ion) energy per nucleon [MeV/u]
 * @param[in] particle_no                particle code number
 * @param[in] material_no                material code number
 * @param[in] rdd_model                  Radial Dose Distribution model code number
 * @param[in] rdd_parameter              Radial Dose Distribution model parameters vector  (array of size 3)
 * @param[in] er_model                   delta electron range model code number
 * @param[in] stopping_power_source_no   TODO
 * @return  Maximum dose [Gy]
 */
double AT_RDD_d_max_Gy( const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no);


/**
 * Pre-calculated many useful parameters characterizing RDD.
 * @param[in]  E_MeV_u               particle (ion) energy per nucleon [MeV/u]
 * @param[in]  particle_no           particle code number
 * @param[in]  material_no           material code number
 * @param[in]  rdd_model             Radial Dose Distribution model code number
 * @param[in]  rdd_parameter         Radial Dose Distribution model parameters vector  (array of size 3)
 * @param[in]  er_model              delta electron range model code number
 * @param[in]  stopping_power_source_no     TODO
 * @param[out] f1_parameters <br>
 *     0 - LET_MeV_cm2_g <br>
 *     1 - r_min_m <br>
 *     2 - r_max_m <br>
 *     3 - d_min_Gy <br>
 *     4 - d_max_Gy <br>
 *     5 - normalization constant [Gy] <br>
 *     6 - single_impact_fluence_cm2 <br>
 *     7 - single_impact_dose_Gy
 */
void AT_RDD_f1_parameters_single_field( const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no,
    double        f1_parameters[]);


/**
 * Pre-calculated many useful parameters characterizing RDD.
 * @param[in]  n                     number of particle types in the mixed particle field
 * @param[in]  E_MeV_u               energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no           type of the particles in the mixed particle field (array of size n)
 * @param[in]  material_no           material code number
 * @param[in]  rdd_model             Radial Dose Distribution model code number
 * @param[in]  rdd_parameter         Radial Dose Distribution model parameters vector  (array of size 3)
 * @param[in]  er_model              delta electron range model code number
 * @param[in]  stopping_power_source_no  TODO
 * @param[out] f1_parameters  (array of size 8)
 *     0 - LET_MeV_cm2_g \n
 *     1 - r_min_m \n
 *     2 - r_max_m \n
 *     3 - d_min_Gy \n
 *     4 - d_max_Gy \n
 *     5 - normalization constant [Gy] \n
 *     6 - single_impact_fluence_cm2 \n
 *     7 - single_impact_dose_Gy
 */
void AT_RDD_f1_parameters_mixed_field( const long    n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    stopping_power_source_no,
    double        f1_parameters[]);


#endif // AT_RDD_H_
