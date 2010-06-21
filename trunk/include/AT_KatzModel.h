#ifndef AT_KATZMODEL_H_
#define AT_KATZMODEL_H_

/**
 * @file
 * @brief Katz model algorithm
 */

/*
 *    AT_KatzModel.h
 *    ===========================
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

#include "AT_RDD.h"
#include "AT_GammaResponse.h"


/**
 * TODO
 * @param r_m
 * @param a0_m
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param er_model
 * @param alpha
 * @param Katz_plateau_Gy
 * @param Katz_point_coeff_Gy
 * @param D0_characteristic_dose_Gy
 * @param c_hittedness
 * @param m_number_of_targets
 * @return
 */
double AT_KatzModel_KatzExtTarget_inactivation_probability(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const long    er_model,
    const double  alpha,
    const double  Katz_plateau_Gy,
    const double  Katz_point_coeff_Gy,
    const double  D0_characteristic_dose_Gy,
    const double  c_hittedness,
    const double  m_number_of_targets);


/**
 * TODO
 * @param r_m
 * @param a0_m
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param beta
 * @param C_norm
 * @param Cucinotta_plateau_Gy
 * @param KatzPoint_point_coeff_Gy
 * @param D0_characteristic_dose_Gy
 * @param c_hittedness
 * @param m_number_of_targets
 * @return
 */
double AT_KatzModel_CucinottaExtTarget_inactivation_probability(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy,
    const double  KatzPoint_point_coeff_Gy,
    const double  D0_characteristic_dose_Gy,
    const double  c_hittedness,
    const double  m_number_of_targets);


/**
 * TODO
 * @param n
 * @param r_m
 * @param E_MeV_u
 * @param particle_no
 * @param material_no
 * @param rdd_model
 * @param rdd_parameters
 * @param er_model
 * @param gamma_parameters
 * @param inactivation_probability
 * @return status code
 */
int AT_KatzModel_inactivation_probability(
    const long    n,
    const double  r_m[],
    const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameters[],
    const long    er_model,
    const double  gamma_parameters[5],
    double        inactivation_probability[]);


/**
 * TODO
 */
typedef struct {
  double  a0_m;
  double  KatzPoint_r_min_m;
  double  max_electron_range_m;
  long    er_model;
  double  alpha;
  double  Katz_plateau_Gy;
  double  Katz_point_coeff_Gy;
  double  D0_characteristic_dose_Gy;
  double  c_hittedness;
  double  m_number_of_targets;
} AT_KatzModel_KatzExtTarget_inactivation_probability_parameters;


/**
 * TODO
 * @param t_m
 * @param params
 * @return
 */
double AT_KatzModel_KatzExtTarget_inactivation_cross_section_integrand_m(
    double t_m,
    void* params);


/**
 * TODO
 * @param a0_m
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param er_model
 * @param alpha
 * @param Katz_plateau_Gy
 * @param Katz_point_coeff_Gy
 * @param D0_characteristic_dose_Gy
 * @param c_hittedness
 * @param m_number_of_targets
 * @return
 */
double AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2(
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const long    er_model,
    const double  alpha,
    const double  Katz_plateau_Gy,
    const double  Katz_point_coeff_Gy,
    const double  D0_characteristic_dose_Gy,
    const double  c_hittedness,
    const double  m_number_of_targets);


/**
 * TODO
 */
typedef struct {
  double  a0_m;
  double  KatzPoint_r_min_m;
  double  max_electron_range_m;
  double  beta;
  double  C_norm;
  double  Cucinotta_plateau_Gy;
  double  KatzPoint_coeff_Gy;
  double  D0_characteristic_dose_Gy;
  double  c_hittedness;
  double  m_number_of_targets;
} AT_KatzModel_CucinottaExtTarget_inactivation_probability_parameters;


/**
 * TODO
 * @param t_m
 * @param params
 * @return
 */
double AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_integrand_m(
    double t_m,
    void* params);


/**
 * TODO
 * @param a0_m
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param beta
 * @param C_norm
 * @param Cucinotta_plateau_Gy
 * @param KatzPoint_point_coeff_Gy
 * @param D0_characteristic_dose_Gy
 * @param c_hittedness
 * @param m_number_of_targets
 * @return
 */
double AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2(
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy,
    const double  KatzPoint_point_coeff_Gy,
    const double  D0_characteristic_dose_Gy,
    const double  c_hittedness,
    const double  m_number_of_targets);


/**
 * TODO
 * @param n
 * @param E_MeV_u
 * @param particle_no
 * @param material_no
 * @param rdd_model
 * @param rdd_parameters
 * @param er_model
 * @param gamma_parameters
 * @param inactivation_cross_section_m2
 * @return status code
 */
int AT_KatzModel_inactivation_cross_section_m2(
    const long   n,
    const double E_MeV_u[],
    const long   particle_no,
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double gamma_parameters[5],
    double inactivation_cross_section_m2[]);


/**
 * TODO
 */
double AT_KatzModel_single_field_survival(
    const double fluence_cm2,
	const double E_MeV_u,
    const long   particle_no,
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double D0_characteristic_dose_Gy,
    const double m_number_of_targets,
    const double kappa);


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
  double*  E_MeV_u;
  long*    particle_no;
  long*    material_no;
  long*    rdd_model;
  double*  rdd_parameters;
  long*    er_model;
  double   gamma_parameters[5];
} AT_P_RDD_parameters;


#endif /* AT_KATZMODEL_H_ */
