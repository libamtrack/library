#ifndef AT_SUCCESSIVECONVOLUTIONS_H_
#define AT_SUCCESSIVECONVOLUTIONS_H_

/**
 * @file
 * @brief Successive Convolution algorithm
 */

/*
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


/**
 * TODO
 */
typedef struct{

  long    array_size;                  /**< size of function arrays F..BI */
  long    N2;
  float   U;

  float   X;
  float   FINAL;

  float   CN;
  long    N1;

  float   CM1;
  float   CM2;
  float   CM3;
  float   CM4;

  float   D1;
  float   D2;
  float   D3;
  float   D4;

  float   F0;
  float*  F;
  long    MIF;
  long    LEF;

  float   H0;
  float*  H;
  long    MIH;
  long    LEH;

  float   E0;
  float*  E;
  float*  DE;
  long    MIE;

  float*  DI;
  float*  A;
  float*  BI;

  bool    write_output;
  FILE*   output_file;

  bool    shrink_tails;
  float   shrink_tails_under;
  bool    adjust_N2;
}     aKList;


/**
 * TODO
 * @param[in]  n                   number of particle types in the mixed particle field
 * @param[in]  E_MeV_u             energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no         type of the particles in the mixed particle field (array of size n)
 * @param[in]  material_no         index number for detector material
 * @param[in]  rdd_model           index number for chosen radial dose distribution
 * @param[in]  rdd_parameters      parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @param[in]  er_model            index number for chosen electron-range model
 * @param[in]  er_parameters       parameters for chosen electron-range model (array of size depending on chosen model)
 * @param[in]  N2                  number of bins per factor of two in local dose array
 * @param[out] n_bins_f1
 * @param[out] f1_parameters
 */
void  AT_SC_get_f1_array_size(  /* radiation field parameters */
    const long   n,
    const float  E_MeV_u[],
    const long   particle_no[],
    /* detector parameters */
    const long   material_no,
    const long   rdd_model,
    const float  rdd_parameter[],
    /* electron range model */
    const long   er_model,
    const float  er_parameter[],
    /* algorithm parameters*/
    const long   N2,
    // from here: return values
    long *       n_bins_f1,
    float        f1_parameters[]);


/**
 * TODO
 * @param[in]  n                     number of particle types in the mixed particle field
 * @param[in]  E_MeV_u               energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no           type of the particles in the mixed particle field (array of size n)
 * @param[in]  fluence_cm2           fluences for the given particles, doses in Gy if negative (array of size n)
 * @param[in]  material_no           index number for detector material
 * @param[in]  rdd_model             index number for chosen radial dose distribution
 * @param[in]  rdd_parameters        parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @param[in]  er_model              index number for chosen electron-range model
 * @param[in]  er_parameters         parameters for chosen electron-range model (array of size depending on chosen model)
 * @param[in]  N2                    number of bins per factor of two in local dose array
 * @param[in]  f1_parameters         TODO
 * @param[out] norm_fluence          TODO
 * @param[out] dose_contribution_Gy  TODO
 * @param[out] f_parameters          TODO
 * @param[out] f1_d_Gy               TODO
 * @param[out] f1_dd_Gy              TODO
 * @param[out] f1                    TODO
 */
void  AT_SC_get_f1(  /* radiation field parameters */
    const long   n,
    const float  E_MeV_u[],
    const long   particle_no[],
    const float  fluence_cm2[],
    /* detector parameters */
    const long   material_no,
    const long   rdd_model,
    const float  rdd_parameter[],
    /* electron range model */
    const long   er_model,
    const float  er_parameter[],
    /* algorithm parameters*/
    const long   N2,
    const long   n_bins_f1,
    /* f1 parameters*/
    const float  f1_parameters[],
    // from here: return values
    float        norm_fluence[],
    float        dose_contribution_Gy[],
    float        f_parameters[],
    /*  1 - total fluence_cm2
     *  2 - total_dose_Gy
     *  3 - ave_E_MeV
     *  4 - dw_E_MeV
     *  5 - ave_LET_MeV_cm2_g
     *  6 - dw_LET_MeV_cm2_g
     *  0 - u
     */
    float  f1_d_Gy[],
    float  f1_dd_Gy[],
    float  f1[]);


/**
 * TODO
 * @param[in]  u
 * @param[in]  fluence_factor
 * @param[in]  N2                  number of bins per factor of two in local dose array
 * @param[in]  n_bins_f1
 * @param[in]  f1_d_Gy
 * @param[in]  f1_dd_Gy
 * @param[in]  f1
 * @param[out] n_bins_f
 * @param[out] u_start
 * @param[out] n_convolutions
 */
void AT_SC_get_f_array_size(  const float   u,
    const float   fluence_factor,
    const long    N2,
    const long    n_bins_f1,
    const float   f1_d_Gy[],
    const float   f1_dd_Gy[],
    const float   f1[],
    // from here: return values
    long*         n_bins_f,
    float*        u_start,
    long*         n_convolutions);


/**
 * TODO
 * @param[in]  u_start             TODO it seems that it is unused in function body
 * @param[in]  n_bins_f1
 * @param[in]  N2                  number of bins per factor of two in local dose array
 * @param[in]  f1_d_Gy
 * @param[in]  f1_dd_Gy
 * @param[in]  f1
 * @param[in]  n_bins_f
 * @param[out] f_d_Gy
 * @param[out] f_dd_Gy
 * @param[out] f_start
 */
void  AT_SC_get_f_start( const float  u_start,
    const long    n_bins_f1,
    const long    N2,
    const float   f1_d_Gy[],
    const float   f1_dd_Gy[],
    const float   f1[],
    const long    n_bins_f,
    // from here: return values
    float   f_d_Gy[],
    float   f_dd_Gy[],
    float   f_start[]);


/**
 * TODO
 * @param theKList
 * @return
 */
aKList  AT_SC_NORMAL(aKList theKList);


/**
 * TODO
 * @param theKList
 * @return
 */
aKList  AT_SC_OUTPUT(aKList theKList);


/**
 * TODO
 * @param theKList
 * @return
 */
aKList  AT_SC_INTERP(aKList theKList);


/**
 * TODO
 * @param theKList
 * @return
 */
aKList  AT_SC_RESET(aKList theKList);


/**
 * TODO
 * @param theKList
 * @return
 */
aKList  AT_SC_ZERO(aKList theKList);


/**
 * TODO
 * @param theKList
 * @return
 */
aKList  AT_SC_SHRINK(aKList theKList);


/**
 * TODO
 * @param theKList
 * @return
 */
aKList AT_SC_FOLD(aKList theKList);


/**
 * TODO
 * @param[in]      u
 * @param[in]      n_bins_f
 * @param[in,out]  N2                    number of bins per factor of two in local dose array
 * @param[in,out]  n_bins_f_used
 * @param[in,out]  f_d_Gy
 * @param[in,out]  f_dd_Gy
 * @param[in,out]  f
 * @param[out]     f0
 * @param[out]     fdd                   frequency:          H * DE        (f * dd)
 * @param[out]     dfdd                  dose contribution:  H * E * DE    (f * d * dd)
 * @param[out]     d                     first moment:                     (<d>)
 * @param[in]      write_output
 * @param[in]      shrink_tails
 * @param[in]      shrink_tails_under
 * @param[in]      adjust_N2
 */
void   AT_SuccessiveConvolutions( const float  u,
    const long   n_bins_f,
    long*        N2,
    long*        n_bins_f_used,
    float        f_d_Gy[],
    float        f_dd_Gy[],
    float        f[],
    float*       f0,
    float        fdd[],
    float        dfdd[],
    float*       d,
    const bool   write_output,
    const bool   shrink_tails,
    const float  shrink_tails_under,
    const bool   adjust_N2);


#endif // AT_SUCCESSIVECONVOLUTIONS_H_
