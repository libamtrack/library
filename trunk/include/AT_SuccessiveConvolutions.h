#ifndef AT_SUCCESSIVECONVOLUTIONS_H_
#define AT_SUCCESSIVECONVOLUTIONS_H_

/**
 * @brief Successive Convolution algorithm
 */

/*
 *    AT_SuccessiveConvolutions.h
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

#include "AT_Constants.h"
#include "AT_RDD.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>

/**
 * Computes the size of the array to hold the f1 (single impact) local dose distribution for a given field, rdd, er
 * @param[in]  n                   number of particle types in the mixed particle field
 * @param[in]  E_MeV_u             energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no         type of the particles in the mixed particle field (array of size n)
 * @param[in]  material_no         index number for detector material
 * @param[in]  rdd_model           index number for chosen radial dose distribution
 * @param[in]  rdd_parameter       parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @param[in]  er_model            index number for chosen electron-range model
 * @param[in]  N2                  number of bins per factor of two in local dose array
 * @return number of bins to hold the f1 distribution
 */
long  AT_n_bins_for_singe_impact_local_dose_distrib(
    const long   n,
    const double E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameter[],
    const long   er_model,
    const long   N2);


/**
 * Computes the f1 (single impact) local dose distribution for a given field, rdd, er
 * Usually step 2 of the CPP-SC method
 * @param[in]  n                     number of particle types in the mixed particle field
 * @param[in]  E_MeV_u               energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no           type of the particles in the mixed particle field (array of size n)
 * @param[in]  fluence_cm2_or_dose_Gy           fluences for the given particles, doses in Gy if negative (array of size n)
 * @param[in]  material_no           index number for detector material
 * @param[in]  rdd_model             index number for chosen radial dose distribution
 * @param[in]  rdd_parameter         parameters for chosen radial dose distribution (array of size depending on chosen model)
 * @param[in]  er_model              index number for chosen electron-range model
 * @param[in]  N2                    number of bins per factor of two in local dose array
 * @param[in]  f1_parameters         array of size n * 8 with n field component characteristics (from AT_SC_get_f1_array_size)
 * @param[in]  n_bins_f1             number of bins holding the f1 distribution (from AT_SC_get_f1_array_size)
 * @param[out] f1_d_Gy               bin midpoints for f1, array of size n_bins_f1
 * @param[out] f1_dd_Gy              bin widths for f1, array of size n_bins_f1
 * @param[out] f1                    f1 values, array of size n_bins_f1
 */
void  AT_single_impact_local_dose_distrib(
    const long   n,
    const double E_MeV_u[],
    const long   particle_no[],
    const double fluence_cm2_or_dose_Gy[],
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameter[],
    const long   er_model,
    const long   N2,
    const long   n_bins_f1,
    const double f1_parameters[],
    double       f1_d_Gy[],
    double       f1_dd_Gy[],
    double       f1[]);


/**
 * Estimates the size of the array to hold the resulting f local dose distribution for a given field, rdd, er
 * Usually step 3 of the CPP-SC method
 * @param[in]  u                   mean number of tracks contributing to the representative point
 * @param[in]  fluence_factor      variable to tweak the total dose from the mixed field (rather than change the single components fluences)
 * @param[in]  N2                  number of bins per factor of two in local dose array
 * @param[in]  n_bins_f1           number of bins holding the f1 distribution (from AT_SC_get_f1_array_size)
 * @param[in]  f1_d_Gy             bin midpoints for f1, array of size n_bins_f1 (from AT_SC_get_f1)
 * @param[in]  f1_dd_Gy            bin width for f1, array of size n_bins_f1 (from AT_SC_get_f1)
 * @param[in]  f1                  f1 values, array of size n_bins_f1 (from AT_SC_get_f1)
 * @param[out] n_bins_f            number of bins holding the resulting f local dose distribution
 * @param[out] u_start             value for u to start convolutions with, between 0.001 and 0.002 where linearization of f into no and one impact is valid
 * @param[out] n_convolutions      number of convolutions necessary to get from u_start to u (u = 2^n_convolutions * u_start)
 */
void AT_n_bins_for_low_fluence_local_dose_distribution(  const double   u,
    const double  fluence_factor,
    const long    N2,
    const long    n_bins_f1,
    const double  f1_d_Gy[],
    const double  f1_dd_Gy[],
    const double  f1[],
    long*         n_bins_f,
    double*       u_start,
    long*         n_convolutions);


/**
 * Computes the (linearized) local dose distribution to start convolutions with, i.e. f at u_start
 * TODO why is it f at u_start ? u_start is not used in function body
 * Usually step 4 of the CPP-SC method
 * @param[in]  n_bins_f1           number of bins holding the f1 distribution (from AT_SC_get_f1_array_size)
 * @param[in]  N2                  number of bins per factor of two in local dose array
 * @param[in]  f1_d_Gy             bin midpoints for f1, array of size n_bins_f1 (from AT_SC_get_f1)
 * @param[in]  f1_dd_Gy            bin width for f1, array of size n_bins_f1 (from AT_SC_get_f1)
 * @param[in]  f1                  f1 values, array of size n_bins_f1 (from AT_SC_get_f1)
 * @param[in]  n_bins_f            number of bins holding the resulting f local dose distribution (from AT_SC_get_f_array_size)
 * @param[out] f_d_Gy              bin midpoints for f, array of size n_bins_f
 * @param[out] f_dd_Gy             bin widths for f, array of size n_bins_f
 * @param[out] f_start             f values to start with, array of size n_bins_f
 */
void  AT_low_fluence_local_dose_distribution(  const long    n_bins_f1,
    const long    N2,
    const double  f1_d_Gy[],
    const double  f1_dd_Gy[],
    const double  f1[],
    const long    n_bins_f,
    double        f_d_Gy[],
    double        f_dd_Gy[],
    double        f_start[]);


/**
 * Routine to perform the convolutions from initial linearized local dose distribution f_start to resulting f
 * as described by Kellerer, 1969. This is a to most extend a reimplementation of Kellerer's original FORTRAN IV code.
 * Usually step 5 of the CPP-SC method
 * @param[in]      u                   value for u to start convolutions with, between 0.001 and 0.002 where linearization of f into no and one impact is valid (from AT_SC_get_f_array_size)
 * @param[in]      n_bins_f            number of bins holding the resulting f local dose distribution (from AT_SC_get_f_array_size)
 * @param[in,out]  N2                  number of bins per factor of two in local dose array f_start, will return new value in case it was adjusted by the routine (higher resolution in case of high fluences)
 * @param[in,out]  n_bins_f_used       in: number of bins used for f_start, \n
 *                                     out: number of bins used for resulting. As tails can be cut and N2 adjusted this is usually not the array size for f_d_Gy, f_dd_Gy, f but smaller (so entries 0..n_bins_f_used-1 are used)
 * @param[in,out]  n_bins_f            number of bins holding the resulting f local dose distribution (from AT_SC_get_f_array_size)
 * @param[in,out]  f_d_Gy              bin midpoints for f, array of size n_bins_f (from AT_SC_get_f_start)
 * @param[in,out]  f_dd_Gy             bin widths for f, array of size n_bins_f (from AT_SC_get_f_start)
 * @param[in,out]  f                   in: f values to start with, array of size n_bins_f (from AT_SC_get_f_start) \n
 *                                     out: resulting f values after convolutions
 * @param[out]     f0                  zero-dose f value (as bins are log this has to be separated)
 * @param[out]     fdd                 frequency:          H * DE        (f * dd), size n_bins_f, used: n_bins_f_used
 * @param[out]     dfdd                dose contribution:  H * E * DE    (f * d * dd), size n_bins_f, used: n_bins_f_used
 * @param[out]     d                   first moment:       (@<d@>)         provides check on convolution quality
 * @param[in]      write_output        if true, a very verbose log file will be written ("SuccessiveConvolutions.log") with results from each convolution
 * @param[in]      shrink_tails        if true, the upper and lower tail of f will be cut after every convolution (bins that contribute less than "shrink_tails_under" to the first moment @<d@>)
 * @param[in]      shrink_tails_under  cut threshold for tails
 * @param[in]      adjust_N2           if true, N2 (i.e. the bin density) can be adjusted. This can be necessary esp. for high doses / fluences where f gets very narrow
 */
void   AT_SuccessiveConvolutions( const double  final_mean_number_of_tracks_contrib,
		const long    n_bins_f,
		long*         N2,
		long*         n_bins_f_used,
		double        f_d_Gy[],
		double        f_dd_Gy[],
		double        f[],
		double*       f0,
		double        fdd[],
		double        dfdd[],
		double*       d,
		const bool    write_output,
		const bool    shrink_tails,
		const double  shrink_tails_under,
		const bool    adjust_N2);

/**
 * Normalized the distribution resulting from last convolution
 */
void AT_Kellerer_normalize(const long array_size,
		const long MIH,
		const long MIE,
		const long LEH,
		const double H0,
		const double E[],
		const double DE[],
		double* CM1,
		double* CM2,
		double H[]);

/**
 * Calculates arrays A and B to be used in quadratic extrapolation of F in AT_SC_FOLD
 */
void AT_Kellerer_interpolation(const long N2,
		const long LEF,
		const long array_size,
		double	F[],
		double	A[],
		double	BI[]);

/**
 * Selects a new coordinate system if F has become to narrow
 */
void AT_Kellerer_reset(long* N2,
		double* U,
		long* array_size,
		long* LEF,
		long* MIE,
		long* MIF,
		const double E0,
		double E[],
		double DE[],
		double F[],
		double A[],
		double BI[],
		double DI[]);


/**
 * Adds the term 2*F0*F(L) to H(L)
 */
void AT_Kellerer_zero(const long MIF,
		const long array_size,
		const long MIE,
		const long LEF,
		const double F0,
		const double F[],
		const double DE[],
		long* MIH,
		long* LEH,
		double H[]);

/**
 * Cuts tails of distribution that contribute less that shrink_tails_under to @<f@>
 */
void AT_Kellerer_shrink(const long array_size,
		const long MIE,
		const double shrink_tails_under,
		const double DE[],
		long* MIH,
		long* LEH,
		double H[]);

/**
 * Does actual convolution, makes use of the symmetry
 */
void AT_Kellerer_folding(		const long N2,
		const long array_size,
		const long LEF,
		const long MIE,
		const long MIF,
		const double DE[],
		const double DI[],
		long* MIH,
		long* LEH,
		const double F0,
		double* H0,
		double F[],
		double H[],
		double A[],
		double BI[]);



#endif // AT_SUCCESSIVECONVOLUTIONS_H_
