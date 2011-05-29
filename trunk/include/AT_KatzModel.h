#ifndef AT_KATZMODEL_H_
#define AT_KATZMODEL_H_

/**
 * @brief Katz model algorithm
 */

/*
 *    AT_KatzModel.h
 *    ===========================
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

#include "AT_RDD.h"
#include "AT_GammaResponse.h"


/**
 * TODO
 * @param[in] r_m
 * @param[in] a0_m
 * @param[in] KatzPoint_r_min_m
 * @param[in] max_electron_range_m
 * @param[in] er_model
 * @param[in] alpha
 * @param[in] Katz_plateau_Gy
 * @param[in] Katz_point_coeff_Gy
 * @param[in] D0_characteristic_dose_Gy
 * @param[in] c_hittedness
 * @param[in] m_number_of_targets
 * @return inactivation probability
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
 * @param[in] r_m
 * @param[in] a0_m
 * @param[in] KatzPoint_r_min_m
 * @param[in] max_electron_range_m
 * @param[in] beta
 * @param[in] C_norm
 * @param[in] Cucinotta_plateau_Gy
 * @param[in] KatzPoint_point_coeff_Gy
 * @param[in] D0_characteristic_dose_Gy
 * @param[in] c_hittedness
 * @param[in] m_number_of_targets
 * @return inactivation probability
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
 * @param[in] n                             todo
 * @param[in] r_m                           todo (array of size n)
 * @param[in] E_MeV_u                       todo
 * @param[in] particle_no                   todo
 * @param[in] material_no                   todo
 * @param[in] rdd_model                     todo
 * @param[in] rdd_parameters                todo (array of size 4)
 * @param[in] er_model                      todo
 * @param[in] gamma_parameters              todo (array of size 5)
 * @param[in] stop_power_source             todo
 * @param[out] inactivation_probability     results (array of size n)
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
    const double  gamma_parameters[],
    const long    stop_power_source,
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
 * @param[in] t_m
 * @param[in] params
 * @return TODO
 */
double AT_KatzModel_KatzExtTarget_inactivation_cross_section_integrand_m(
    double t_m,
    void* params);


/**
 * TODO
 * @param[in] a0_m
 * @param[in] KatzPoint_r_min_m
 * @param[in] max_electron_range_m
 * @param[in] er_model
 * @param[in] alpha
 * @param[in] Katz_plateau_Gy
 * @param[in] Katz_point_coeff_Gy
 * @param[in] D0_characteristic_dose_Gy
 * @param[in] c_hittedness
 * @param[in] m_number_of_targets
 * @return TODO
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
 * @param[in] t_m
 * @param[in] params
 * @return TODO
 */
double AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_integrand_m(
    double t_m,
    void* params);


/**
 * TODO
 * @param[in] a0_m
 * @param[in] KatzPoint_r_min_m
 * @param[in] max_electron_range_m
 * @param[in] beta
 * @param[in] C_norm
 * @param[in] Cucinotta_plateau_Gy
 * @param[in] KatzPoint_point_coeff_Gy
 * @param[in] D0_characteristic_dose_Gy
 * @param[in] c_hittedness
 * @param[in] m_number_of_targets
 * @return TODO
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
 * @param[in] n
 * @param[in] E_MeV_u (array of size n)
 * @param[in] particle_no
 * @param[in] material_no
 * @param[in] rdd_model
 * @param[in] rdd_parameters (array of size 4)
 * @param[in] er_model
 * @param[in] gamma_parameters (array of size 5)
 * @param[in] stop_power_source             todo
 * @param[out] inactivation_cross_section_m2 (array of size n)
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
    const double gamma_parameters[],
    const long   stop_power_source,
    double inactivation_cross_section_m2[]);


/**
 * TODO
 * @param[in] fluence_cm2
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @param[in] inactivation_cross_section_m2
 * @param[in] D0_characteristic_dose_Gy
 * @param[in] m_number_of_targets
 * @param[in] sigma0_m2
 * @param[in] stopping_power_source_no
 * @return TODO
 */
double AT_KatzModel_single_field_survival_from_inactivation_cross_section(
    const double fluence_cm2,
	const double E_MeV_u,
    const long   particle_no,
    const long   material_no,
    const double inactivation_cross_section_m2,
    const double D0_characteristic_dose_Gy,
    const double m_number_of_targets,
    const double sigma0_m2,
    const long   stopping_power_source_no);


/**
 * TODO
 * @param[in] fluence_cm2
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @param[in] rdd_model
 * @param[in] rdd_parameters (array of size 4)
 * @param[in] er_model
 * @param[in] D0_characteristic_dose_Gy
 * @param[in] m_number_of_targets
 * @param[in] sigma0_m2
 * @param[in] stopping_power_source_no
 * @param[out] survival
 * @return status code
 */
int AT_KatzModel_single_field_survival(
    const double fluence_cm2,
	const double E_MeV_u,
    const long   particle_no,
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double D0_characteristic_dose_Gy,
    const double m_number_of_targets,
    const double sigma0_m2,
    const long   stopping_power_source_no,
    double*     survival);


/**
 * TODO
 * @param[in] number_of_items
 * @param[in] fluence_cm2 (array of size number_of_items)
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] material_no
 * @param[in] rdd_model
 * @param[in] rdd_parameters (array of size 4)
 * @param[in] er_model
 * @param[in] D0_characteristic_dose_Gy
 * @param[in] m_number_of_targets
 * @param[in] sigma0_m2
 * @param[in] stopping_power_source_no
 * @param[out] survival (array of size number_of_items)
 * @return
 */
int AT_KatzModel_single_field_survival_optimized_for_fluence_vector(
	const long   number_of_items,
    const double fluence_cm2[],
	const double E_MeV_u,
    const long   particle_no,
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double D0_characteristic_dose_Gy,
    const double m_number_of_targets,
    const double sigma0_m2,
    const long   stopping_power_source_no,
    double* survival);

/**
 * TODO
 */
double AT_P_RDD(                    double  r_m,
    void* params);


/**
 * TODO
 */
double AT_sI_int(                   double  r_m,
    void* params);


/**
 * TODO
 */
double AT_D_RDD_Gy_int(             double  r_m,
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
