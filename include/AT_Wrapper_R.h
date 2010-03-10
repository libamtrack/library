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

#include "AT_PhysicsRoutines.h"
#include "AT_ElectronRange.h"
#include "AT_RDD.h"
#include "AT_RDD_ExtendedTarget.h"
#include "AT_GammaResponse.h"

void AT_max_E_transfer_MeV_R(  const int*  n,
    const float*  E_MeV_u,
    // results
    float*  max_E_transfer_MeV);

void AT_max_electron_range_m_R(  const int*  n,
    const float*  E_MeV_u,
    const int*  material_no,
    const int*   er_model,
    // results
    float*  max_electron_range_m);

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

void AT_RDD_ExtendedTarget_Gy_R( const int*  n,
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

#endif /* AT_WRAPPER_R_H_ */
