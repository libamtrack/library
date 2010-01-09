#ifndef AT_SUCCESSIVECONVOLUTIONS_H_
#define AT_SUCCESSIVECONVOLUTIONS_H_

/**
 *    AT_SuccessiveConvolutions.h
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

#include "AT_Constants.h"
#include "AT_RDD.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <malloc.h>

extern int indent_counter;
extern char isp[];
extern FILE * debf;

void  AT_SC_get_f1_array_size(  /* radiation field parameters */
                  long*  n,
                  float*  E_MeV_u,
                  long*  particle_no,
                  /* detector parameters */
                  long*  material_no,
                  long*  rdd_model,
                  float*  rdd_parameter,
                  /* electron range model */
                  long*  er_model,
                  float*  er_parameter,
                  /* algorithm parameters*/
                  long*  N2,
                  // from here: return values
                  long*  n_bins_f1,
                  float*  f1_parameters);


void  AT_SC_get_f1(  /* radiation field parameters */
            long*  n,
            float*  E_MeV_u,
            long*  particle_no,
            float*  fluence_cm2,
            /* detector parameters */
            long*  material_no,
            long*  rdd_model,
            float*  rdd_parameter,
            /* electron range model */
            long*  er_model,
            float*  er_parameter,
            /* algorith parameters*/
            long*  N2,
            long*  n_bins_f1,
            /* f1 parameters*/
            float*  f1_parameters,
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

void AT_SC_get_f_array_size(  float*  u,
    float*  fluence_factor,
    long*  N2,
    long*  n_bins_f1,
    float*  f1_d_Gy,
    float*  f1_dd_Gy,
    float*  f1,
    // from here: return values
    long*  n_bins_f,
    float*  u_start,
    long*  n_convolutions);


void  AT_SC_get_f_start(  float*  u_start,
    long*  n_bins_f1,
    long*  N2,
    float*  f1_d_Gy,
    float*  f1_dd_Gy,
    float*  f1,
    long*  n_bins_f,
    // from here: return values
    float*  f_d_Gy,
    float*  f_dd_Gy,
    float*  f_start);


void AT_SuccessiveConvolutions(  float*  u,
    long*  n_bins_f,
    long*  N2,
    // input + return values
    long*  n_bins_f_used,
    float*  f_d_Gy,
    float*  f_dd_Gy,
    float*  f,
    // return values
    float*  f0,
    float*  fdd,
    float*  dfdd,
    float*  d,
    bool*  write_output,
    bool*  shrink_tails,
    float*  shrink_tails_under,
    bool*  adjust_N2);



#endif // AT_SUCCESSIVECONVOLUTIONS_H_
