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
    /* radiation field parameters */
    const float*  E_MeV_u,
    const int*  particle_no,
    /* detector parameters */
    const int*  material_no,
    /* radial dose distribution model */
    const int*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const int*  er_model,
    const float*  er_parameter,
    float*  D_RDD_Gy);

void AT_D_RDD_ExtendedTarget_Gy_R( const int*  n,
    const float* r_m,
    const float*  a0_m,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const int*    particle_no,
    /* detector parameters */
    const int*    material_no,
    /* radial dose distribution model */
    const int*    rdd_model,
    const float* rdd_parameter,
    /* electron range model */
    const int*    er_model,
    const float* er_parameter,
    float*       D_RDD_Gy);

void AT_gamma_response_R( const int*  n,
    const float*  d_Gy,
    const int*  gamma_model,
    const float*  gamma_parameter,
    float*  S);

void AT_LET_MeV_cm2_g_R(  const int*  n,
    const float*  E_MeV_u,
    const int*   particle_no,
    const int*   material_no,
    float*       LET_MeV_cm2_g);

void AT_max_E_transfer_MeV_R(  const int*  n,
    const float*  E_MeV_u,
    // results
    float*  max_E_transfer_MeV);

void AT_max_electron_ranges_m_R(  const int*  number_of_particles,
    const float*  E_MeV_u,
    const int*    material_no,
    const int*    er_model,
    // results
    float*  max_electron_range_m);

void AT_r_RDD_m_R  ( const int*  n,
    const float*  D_RDD_Gy,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const int*  particle_no,
    /* detector parameters */
    const int*  material_no,
    /* radial dose distribution model */
    const int*  rdd_model,
    const float*  rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
    /* electron range model */
    const int*  er_model,
    const float*  er_parameter,
    float*  r_RDD_m);

void AT_run_GSM_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*  particle_no,
    const float*  fluence_cm2,
    const int*  material_no,
    const int*  RDD_model,
    const float*  RDD_parameters,
    const int*  ER_model,
    const float*  ER_parameters,
    const int*  gamma_model,
    const float*  gamma_parameters,
    const int*  N_runs,
    const int*   N2,
    const float*  fluence_factor,
    const int*   write_output,
    const int*   nX,
    const float*  voxel_size_m,
    const int*   lethal_events_mode,
    float*  results);

void AT_run_SPIFF_method_R(  const int*  n,
    const float*  E_MeV_u,
    const int*  particle_no,
    const float*  fluence_cm2,
    const int*  material_no,
    const int*  RDD_model,
    const float*  RDD_parameters,
    const int*  ER_model,
    const float*  ER_parameters,
    const int*  gamma_model,
    const float*  gamma_parameters,
    int*  N2, // TODO investigate if this can be changed inside
    const float*  fluence_factor,
    const int*  write_output,
    const int*  shrink_tails,
    const float*  shrink_tails_under,
    const int*  adjust_N2,
    const int*   lethal_events_mode,
    float*  results);

void  AT_SC_get_f1_array_size_R(  /* radiation field parameters */
    const int*  n,
    const float*  E_MeV_u,
    const int*  particle_no,
    /* detector parameters */
    const int*  material_no,
    const int*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const int*  er_model,
    const float*  er_parameter,
    /* algorithm parameters*/
    const int*  N2,
    // from here: return values
    int*  n_bins_f1,
    float*  f1_parameters);

void  AT_SC_get_f1_R(  /* radiation field parameters */
    const int*  n,
    const float*  E_MeV_u,
    const int*  particle_no,
    const float*  fluence_cm2,
    /* detector parameters */
    const int*  material_no,
    const int*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const int*  er_model,
    const float*  er_parameter,
    /* algorithm parameters*/
    const int*  N2,
    const int*  n_bins_f1,
    /* f1 parameters*/
    const float*  f1_parameters,
    // from here: return values
    float*  norm_fluence,
    float*  dose_contribution_Gy,
    float*  f_parameters,
    /*  1 - total fluence_cm2
     *  2 - total_dose_Gy
     *  3 - ave_E_MeV
     *  4 - dw_E_MeV
     *  5 - ave_LET_MeV_cm2_g
     *  6 - dw_LET_MeV_cm2_g
     *  0 - u
     */
    float*  f1_d_Gy,
    float*  f1_dd_Gy,
    float*  f1);

void  AT_SC_get_f_array_size_R(
    const float* u,
    const float* fluence_factor,
    const int* N2,
    const int* n_bins_f1,
    const float* f1_d_Gy,
    const float* f1_dd_Gy,
    const float* f1,
    // from here: return values
    int*  n_bins_f,
    float*  u_start,
    int* n_convolutions);

void  AT_SC_get_f_start_R(  const float*  u_start,
    const int*   n_bins_f1,
    const int*   N2,
    const float*  f1_d_Gy,
    const float*  f1_dd_Gy,
    const float*  f1,
    const int*   n_bins_f,
    // from here: return values
    float*  f_d_Gy,
    float*  f_dd_Gy,
    float*  f_start);

void AT_SuccessiveConvolutions_R( const float*  u,
    const int*  n_bins_f,
    // input + return values
    int*  N2,
    int*  n_bins_f_used,
    float*  f_d_Gy,
    float*  f_dd_Gy,
    float*  f,
    // return values
    float*  f0,
    float*  fdd,
    float*  dfdd,
    float*  d,
    const int*  write_output,
    const int*  shrink_tails,
    const float*  shrink_tails_under,
    const int*  adjust_N2);

#endif /* AT_WRAPPER_R_H_ */
