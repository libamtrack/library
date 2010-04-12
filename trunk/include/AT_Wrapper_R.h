#ifndef AT_WRAPPER_R_H_
#define AT_WRAPPER_R_H_

/**
 * @file
 * @brief Wrapper functions
 */

/*
 *    AT_Wrapper_R.h
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


#include "AmTrack.h"
#include "AT_PhysicsRoutines.h"
#include "AT_ElectronRange.h"
#include "AT_RDD.h"
#include "AT_RDD_ExtendedTarget.h"
#include "AT_GammaResponse.h"
#include "AT_SuccessiveConvolutions.h"


void AT_D_RDD_Gy_R( const int*  n,
    const float*  r_m,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    float*        D_RDD_Gy);


void AT_gamma_response_R( const int*  n,
    const float*  d_Gy,
    const int*    gamma_model,
    const float*  gamma_parameter,
    float*        S);


void AT_LET_MeV_cm2_g_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    float*        LET_MeV_cm2_g);


void AT_max_E_transfer_MeV_R(  const int*  n,
    const float*  E_MeV_u,
    float*        max_E_transfer_MeV);


void AT_max_electron_ranges_m_R(  const int*  number_of_particles,
    const float*  E_MeV_u,
    const int*    material_no,
    const int*    er_model,
    float*        max_electron_range_m);


void AT_r_RDD_m_R  ( const int*  n,
    const float*  D_RDD_Gy,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    float*        r_RDD_m);


void AT_KatzModel_inactivation_probability_R( const int* n,
    const float*   r_m,
    const float*   E_MeV_u,
    const int*     particle_no,
    const int*     material_no,
    const int*     rdd_model,
    const float*   rdd_parameters,
    const int*     er_model,
    const float*   gamma_parameters,
    float*         inactivation_probability);


void AT_KatzModel_inactivation_cross_section_m2_R( const int* n,
    const float*   E_MeV_u,
    const int*     particle_no,
    const int*     material_no,
    const int*     rdd_model,
    const float*   rdd_parameters,
    const int*     er_model,
    const float*   gamma_parameters,
    float*         inactivation_cross_section_m2);


void AT_particle_name_from_particle_no_R(const int* particle_no,
    char** particle_name);


void AT_particle_no_from_particle_name_R(const int* particle_no,
    char** particle_name);


void AT_run_GSM_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    const int*    RDD_model,
    const float*  RDD_parameters,
    const int*    ER_model,
    const int*    gamma_model,
    const float*  gamma_parameters,
    const int*    N_runs,
    const int*    N2,
    const float*  fluence_factor,
    const int*    write_output,
    const int*    nX,
    const float*  voxel_size_m,
    const int*    lethal_events_mode,
    float*        results);


void AT_run_SPIFF_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    const int*    RDD_model,
    const float*  RDD_parameters,
    const int*    ER_model,
    const int*    gamma_model,
    const float*  gamma_parameters,
    int*          N2,
    const float*  fluence_factor,
    const int*    write_output,
    const int*    shrink_tails,
    const float*  shrink_tails_under,
    const int*    adjust_N2,
    const int*    lethal_events_mode,
    float*        results);


void  AT_SC_get_f1_array_size_R(
    const int*    n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    const int*    N2,
    int*          n_bins_f1,
    float*        f1_parameters);


void  AT_SC_get_f1_R(
    const int*    n,
    const float*  E_MeV_u,
    const int*    particle_no,
    const float*  fluence_cm2,
    const int*    material_no,
    const int*    rdd_model,
    const float*  rdd_parameter,
    const int*    er_model,
    const int*    N2,
    const int*    n_bins_f1,
    const float*  f1_parameters,
    float*        norm_fluence,
    float*        dose_contribution_Gy,
    float*        f_parameters,
    float*        f1_d_Gy,
    float*        f1_dd_Gy,
    float*        f1);


void  AT_SC_get_f_array_size_R(
    const float*  u,
    const float*  fluence_factor,
    const int*    N2,
    const int*    n_bins_f1,
    const float*  f1_d_Gy,
    const float*  f1_dd_Gy,
    const float*  f1,
    int*          n_bins_f,
    float*        u_start,
    int*          n_convolutions);


void  AT_SC_get_f_start_R(  const int*    n_bins_f1,
    const int*    N2,
    const float*  f1_d_Gy,
    const float*  f1_dd_Gy,
    const float*  f1,
    const int*    n_bins_f,
    float*        f_d_Gy,
    float*        f_dd_Gy,
    float*        f_start);


void AT_SuccessiveConvolutions_R( const float*  u,
    const int*    n_bins_f,
    int*          N2,
    int*          n_bins_f_used,
    float*        f_d_Gy,
    float*        f_dd_Gy,
    float*        f,
    float*        f0,
    float*        fdd,
    float*        dfdd,
    float*        d,
    const int*    write_output,
    const int*    shrink_tails,
    const float*  shrink_tails_under,
    const int*    adjust_N2);

#endif /* AT_WRAPPER_R_H_ */
