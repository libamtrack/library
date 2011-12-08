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
#include "AT_Histograms.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
	#include <sys/malloc.h>
#else
	#include <malloc.h>
#endif

/**
 * Computes the size of the array to hold the f1 (single impact) local dose distribution for a given field, rdd, er
 * @param[in]  n                   number of particle types in the mixed particle field
 * @param[in]  E_MeV_u             energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no         type of the particles in the mixed particle field (array of size n)
 * @param[in]  material_no         index number for detector material
 * @param[in]  rdd_model           index number for chosen radial dose distribution
 * @param[in]  rdd_parameter       parameters for chosen radial dose distribution (array of size 4)
 * @param[in]  er_model            index number for chosen electron-range model
 * @param[in]  N2                  number of bins per factor of two in local dose array
 * @param[in]  stopping_power_source_no  TODO
 * @return number of bins to hold the f1 distribution
 */
long AT_n_bins_for_single_impact_local_dose_distrib( const long   n,
    const double E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    const long   rdd_model,
    const double rdd_parameter[],
    const long   er_model,
    const long   N2,
    const long   stopping_power_source_no);


/**
 * Computes the f1 (single impact) local dose distribution for a given field,
 * radial dose distribution and electron-range model. The routine can handle
 * mixed fields by adding the single contributions according to their respective
 * fluence.
 * This routine can handle all monotonously falling RDDs which should be the case
 * for every realistic approach.
 * Usually step 2 of the CPP-SC method
 * @param[in]  n                     number of particle types in the mixed particle field
 * @param[in]  E_MeV_u               energy of particles in the mixed particle field (array of size n)
 * @param[in]  particle_no           type of the particles in the mixed particle field (array of size n)
 * @param[in]  fluence_cm2_or_dose_Gy           fluences for the given particles, doses in Gy if negative (array of size n)
 * @param[in]  material_no           index number for detector material
 * @param[in]  rdd_model             index number for chosen radial dose distribution
 * @param[in]  rdd_parameter         parameters for chosen radial dose distribution (array of size 4)
 * @param[in]  er_model              index number for chosen electron-range model
 * @param[in]  N2                    number of bins per factor of two in local dose array
 * @param[in]  f1_parameters         n field component characteristics (array of size 8) // TODO in R package should also work for 8*n, but is not working
 * @param[in]  n_bins_f1             number of bins holding the f1 distribution
 * @param[in]  stopping_power_source_no     TODO
 * @param[out] f1_d_Gy               bin midpoints for f1 (array of size n_bins_f1)
 * @param[out] f1_dd_Gy              bin widths for f1 (array of size n_bins_f1)
 * @param[out] f1                    f1 values (array of size n_bins_f1)
 */
void AT_single_impact_local_dose_distrib( const long   n,
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
    const long   stopping_power_source_no,
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
 * @param[in]  f1_d_Gy             bin midpoints for f1 (array of size n_bins_f1)
 * @param[in]  f1_dd_Gy            bin width for f1 (array of size n_bins_f1)
 * @param[in]  f1                  f1 values (array of size n_bins_f1)
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
 * @param[in]  f1_d_Gy             bin midpoints for f1 (array of size n_bins_f1)
 * @param[in]  f1_dd_Gy            bin width for f1 (array of size n_bins_f1)
 * @param[in]  f1                  f1 values (array of size n_bins_f1)
 * @param[in]  n_bins_f            number of bins holding the resulting f local dose distribution (from AT_SC_get_f_array_size)
 * @param[out] f_d_Gy              bin midpoints for f (array of size n_bins_f)
 * @param[out] f_dd_Gy             bin widths for f (array of size n_bins_f)
 * @param[out] f_start             f values to start with (array of size n_bins_f)
 */
void AT_low_fluence_local_dose_distribution(  const long    n_bins_f1,
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
 *
 * @param[in]      final_mean_number_of_tracks_contrib                   value for u to start convolutions with, between 0.001 and 0.002 where linearization of f into no and one impact is valid (from AT_SC_get_f_array_size)
 * @param[in]      n_bins              Size of arrays convolutions are performed on
 * @param[in,out]  N2                  number of bins per factor of two in local dose array f_start, will return new value in case it was adjusted by the routine (higher resolution in case of high fluences)
 * @param[in,out]  n_bins_f_used       in number of bins used for f_start, out number of bins used for resulting. As tails can be cut and N2 adjusted this is usually not the array size for f_d_Gy, f_dd_Gy, f but smaller (so entries 0 to n_bins_f_used-1 are used)
 * @param[in,out]  f_d_Gy              bin midpoints for f (array of size n_bins)
 * @param[in,out]  f_dd_Gy             bin widths for f (array of size n_bins)
 * @param[in,out]  f                   in low fluence approx values to start with, out resulting values after convolutions (array of size n_bins)
 * @param[out]     f0                  zero-dose f value (as bins are log this has to be separated)
 * @param[out]     fdd                 frequency f x f_dd_Gy precomputed for comfort (array of size n_bins)
 * @param[out]     dfdd                dose contribution f x f_dd_Gy precomputed for comfort (array of size n_bins)
 * @param[out]     d                   first moment of distribution f - should coincide with given dose and provides check on convolution quality
 * @param[in]      write_output        if true, a very verbose log file will be written ("SuccessiveConvolutions.log") with results from each convolution
 * @param[in]      shrink_tails        if true, the upper and lower tail of f will be cut after every convolution (bins that contribute less than "shrink_tails_under" to the first moment @<d@>)
 * @param[in]      shrink_tails_under  cut threshold for tails
 * @param[in]      adjust_N2           if true, N2 (i.e. the bin density) can be adjusted. This can be necessary esp. for high doses or fluences where f gets very narrow
 */
void AT_SuccessiveConvolutions( const double  final_mean_number_of_tracks_contrib,
		const long    n_bins,
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
 *
 * @param[in]      values_first_bin    TODO
 * @param[in]      value_midpoints    TODO (array of size frequency_n_bins)
 * @param[in]      value_bin_widths    TODO (array of size frequency_n_bins)
 * @param[in]      frequency_n_bins    TODO
 * @param[in]      frequency_zero_bin    TODO
 * @param[in]      frequency_first_bin    TODO
 * @param[out]     frequency               TODO (array of size frequency_n_bins)
 */
void AT_Kellerer_normalize( const long values_first_bin,
		const double value_midpoints[],
		const double value_bin_widths[],
		const long frequency_n_bins,
		const double frequency_zero_bin,
		const long frequency_first_bin,
		double frequency[]);

/**
 * Calculates arrays A and B to be used in quadratic extrapolation of F in AT_SC_FOLD
 *
 * @param[in]       N2             TODO
 * @param[in]       LEF             TODO
 * @param[in]       array_size             TODO
 * @param[in]       F             TODO (array of size array_size)
 * @param[out]       A             TODO (array of size array_size)
 * @param[out]       BI             TODO (array of size array_size)
 */
void AT_Kellerer_interpolation( const long N2,
		const long LEF,
		const long array_size,
		double	F[],
		double	A[],
		double	BI[]);


/**
 * Selects a new coordinate system if F has become to narrow
 *
 * @param[in]       N2             TODO
 * @param[in]       array_size             TODO
 * @param[in,out]       LEF             TODO
 * @param[in,out]       MIE             TODO
 * @param[in,out]       MIF             TODO
 * @param[in]       E0             TODO
 * @param[in,out]       E             TODO (array of size array_size)
 * @param[in,out]       DE             TODO (array of size array_size)
 * @param[in,out]       F             TODO (array of size array_size)
 * @param[in,out]       DI             TODO (array of size array_size)
 */
void AT_Kellerer_reset( long* N2,
		const long array_size,
		long* LEF,
		long* MIE,
		long* MIF,
		const double E0,
		double E[],
		double DE[],
		double F[],
		double DI[]);


/**
 * Adds the term 2*F0*F(L) to H(L)
 *
 * @param[in]       MIF             TODO
 * @param[in]       array_size             TODO
 * @param[in]       MIE             TODO
 * @param[in]       LEF             TODO
 * @param[in]       F0             TODO
 * @param[in]       F             TODO (array of size array_size)
 * @param[in]       DE             TODO (array of size array_size)
 * @param[in,out]       MIH             TODO
 * @param[in,out]       LEH             TODO
 * @param[in,out]       H             TODO (array of size array_size)
 */
void AT_Kellerer_zero( const long MIF,
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
 *
 * @param[in]       array_size             TODO
 * @param[in]       MIE             TODO
 * @param[in]       shrink_tails_under             TODO
 * @param[in]       DE             TODO (array of size array_size)
 * @param[in,out]       MIH             TODO
 * @param[in,out]       LEH             TODO
 * @param[in,out]       H             TODO (array of size array_size)
 */
void AT_Kellerer_shrink( const long array_size,
		const long MIE,
		const double shrink_tails_under,
		const double DE[],
		long* MIH,
		long* LEH,
		double H[]);


/**
 * Does actual convolution, makes use of the symmetry
 *
 * @param[in]       n_bins             TODO
 * @param[in]       bins_per_factor_2             TODO
 * @param[in]       delta_i             TODO (array of size array_size)
 * @param[in]       values_first_bin             TODO
 * @param[in]       values_bin_widths             TODO (array of size array_size)
 * @param[in]       frequency_n_bins_last             TODO
 * @param[in]       frequency_first_bin_last             TODO
 * @param[in]       frequency_zero_bin_last             TODO
 * @param[in,out]   frequency_last             TODO (array of size array_size)
 * @param[in,out]   frequency_n_bins             TODO
 * @param[in,out]   frequency_first_bin             TODO
 * @param[in,out]   frequency_zero_bin             TODO
 * @param[in,out]   frequency             TODO (array of size array_size)
 */
void AT_Kellerer_folding( const long n_bins,
		const long bins_per_factor_2,
		const double delta_i[],
		const long values_first_bin,
		const double values_bin_widths[],
		const long frequency_n_bins_last,
		const long frequency_first_bin_last,
		const double frequency_zero_bin_last,
		double frequency_last[],
		long* frequency_n_bins,
		long* frequency_first_bin,
		double* frequency_zero_bin,
		double frequency[]);


#endif // AT_SUCCESSIVECONVOLUTIONS_H_
