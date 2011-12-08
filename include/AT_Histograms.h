#ifndef AT_HISTOGRAMS_H_
#define AT_HISTOGRAMS_H_

/**
 * @brief Histograms handling
 */


/*
 *    AT_Constants.h
 *    ==============
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

#include <math.h>
#include <stdlib.h>
#include <assert.h>

/**
 * @enum AT_histo_type
 * Histogram type enumerator
 */
enum AT_histo_type{
  AT_histo_linear       = 0, /**< Histogram with linear bins and arithmetic midpoints, step is the bin limit difference **/
  AT_histo_log          = 1  /**< Histogram with logarithmic bins and geometrical midpoints, step is the bin limit factor **/
};


///////////////////////////////// Left limit routines ////////////////////////////////////
/**
 * Returns left limit of bin with given index number for linear histograms
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  step                  step between bin limits
 * @param[in]  bin_no                index number of bin (zero-based)
 * @param[out] left_limit            left limit of bin
 * @return     status code
 */
int AT_histo_linear_left_limit( const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no,
                                double* left_limit);


/**
 * Returns left limit of bin with given index number for logarithmic histograms
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  step                  step between bin limits
 * @param[in]  bin_no                index number of bin (zero-based)
 * @param[out] left_limit            left limit of bin
 * @return     status code
 */
int AT_histo_logarithmic_left_limit(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no,
                                double* left_limit);


/**
 * Returns left limit of bin with given index number
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  step                  step between bin limits
 * @param[in]  bin_no                index number of bin (zero-based)
 * @param[in]  histo_type            type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] left_limit            left limit of bin
 * @return                      status code
 */
int AT_histo_left_limit( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const long bin_no,
    double* left_limit);


/**
 * Returns vector with all left limits of a histogram
 * @param[in] number_of_bins            number of bins in histogram
 * @param[in] lowest_left_limit         left limit of first bin
 * @param[in] step                      step between bin limits
 * @param[in] histo_type                type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] left_limits              array of bin width
 * length number_of_bins + 1
 * @return status code
 */
int AT_histo_left_limits( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double left_limits[]);


///////////////////////////////// Bin widths routines ////////////////////////////////////
/**
 * Returns width of bin with given index number for linear histograms
 * @param[in] number_of_bins        number of bins in histogram
 * @param[in] lowest_left_limit     left limit of first bin
 * @param[in] step                  step between bin limits
 * @param[in] bin_no                index number of bin (zero-based)
 * @return                      width of bin
 */
double AT_histo_linear_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no);


/**
 * Returns width of bin with given index number for logarithmic histograms
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  step                  step between bin limits
 * @param[in]  bin_no                index number of bin (zero-based)
 * @param[out] bin_width             width of bin
 * @return status code
 */
int AT_histo_logarithmic_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no,
                                double* bin_width);


/**
 * Returns width of bin with given index number
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  step                  step between bin limits
 * @param[in]  bin_no                index number of bin (zero-based)
 * @param[in] histo_type             type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] bin_width             width of bin
 * @return status code
 */
int AT_histo_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long histo_type,
                                const long bin_no,
                                double* bin_width);


/**
 * Returns vector with all bin widths of a histogram
 * @param[in] number_of_bins            number of bins in histogram
 * @param[in] lowest_left_limit         left limit of first bin
 * @param[in] step                      step between bin limits
 * @param[in] histo_type                type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] bin_widths              array of bin width
 * length number_of_bins
 * @return status code
 */
int AT_histo_bin_widths( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double bin_widths[]);


///////////////////////////////// Midpoint routines ////////////////////////////////////
/**
 * Returns midpoint of bin with given index number for linear histograms
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  step                  step between bin limits
 * @param[in]  bin_no                index number of bin (zero-based)
 * @param[out] midpoint              midpoint of bin
 * @return status code
 */
int AT_histo_linear_midpoint(   const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no,
                                double* midpoint);


/**
 * Returns midpoint of bin with given index number for logarithmic histograms
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  step                  step between bin limits
 * @param[in]  bin_no                index number of bin (zero-based)
 * @param[out] midpoint              midpoint of bin
 * @return status code
 */
int AT_histo_logarithmic_midpoint(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no,
                                double* midpoint);


/**
 * Returns midpoint of bin with given index number
 * @param[in]  number_of_bins       number of bins in histogram
 * @param[in]  lowest_left_limit    left limit of first bin
 * @param[in]  step                 step between bin limits
 * @param[in]  bin_no               index number of bin (zero-based)
 * @param[in] histo_type            type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] midpoint             midpoint of bin
 * @return status code
 */
int AT_histo_midpoint(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long histo_type,
                                const long bin_no,
                                double* midpoint);


/**
 * Returns vector with all midpoints of a histogram
 * @param[in] number_of_bins            number of bins in histogram
 * @param[in] lowest_left_limit         left limit of first bin
 * @param[in] step                      step between bin limits
 * @param[in] histo_type                type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] midpoints                array of midpoints
 * length number_of_bins
 * @return status code
 */
int AT_histo_midpoints( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double midpoints[]);

///////////////////////////////// Step routines ////////////////////////////////////
/**
 * Returns step for given range for linear histograms
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  highest_left_limit    right (exclusive) limit of last bin
 * @param[out] step                  step between bins
 * @return status code
 */
int AT_histo_linear_step(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double highest_left_limit,
                                double* step);

/**
 * Returns step for given range for logarithmic histograms
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  highest_left_limit    right (exclusive) limit of last bin
 * @param[out] step                  step between bins
 * @return status code
 */
int AT_histo_logarithmic_step(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double highest_left_limit,
                                double* step);

/**
 * Returns step for given range
 * @param[in]  number_of_bins        number of bins in histogram
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  highest_left_limit    right (exclusive) limit of last bin
 * @param[in] histo_type             type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] step                  step between bins
 * @return status code
 */
int AT_histo_step(      const long number_of_bins,
                        const double lowest_left_limit,
                        const double highest_left_limit,
                        const long histo_type,
                        double* step);

///////////////////////////////// Number of bin routines ////////////////////////////////////
/**
 * Returns number of bins for given range for linear histograms
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  highest_left_limit    right (exclusive) limit of last bin
 * @param[in]  step                  step between bins
 * @param[out] number_of_bins        number of bins in histogram
 * @return status code
 */
int AT_histo_linear_n_bins(     const double lowest_left_limit,
                                const double highest_left_limit,
                                const double step,
                                long* number_of_bins);

/**
 * Returns number of bins for given range for logarithmic histograms
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  highest_left_limit    right (exclusive) limit of last bin
 * @param[in]  step                  step between bins
 * @param[out] number_of_bins        number of bins in histogram
 * @return status code
 */
int AT_histo_logarithmic_n_bins(     const double lowest_left_limit,
                                const double highest_left_limit,
                                const double step,
                                long* number_of_bins);

/**
 * Returns number of bins for given range for logarithmic histograms
 * @param[in]  lowest_left_limit     left limit of first bin
 * @param[in]  highest_left_limit    right (exclusive) limit of last bin
 * @param[in]  step                  step between bins
 * @param[in] histo_type             type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] number_of_bins        number of bins in histogram
 * @return status code
 */
int AT_histo_n_bins(     const double lowest_left_limit,
    const double highest_left_limit,
    const double step,
    const long histo_type,
    long* number_of_bins);

///////////////////////////////// Access routines ////////////////////////////////////
/**
 * Returns bin index number for a given value for linear histograms
 * @param[in] number_of_bins        number of bins in histogram
 * @param[in] lowest_left_limit     left limit of first bin
 * @param[in] step                  step between bin limits
 * @param[in] value                 value
 * @return                      bin index number (zero-based, positive), -1 if value below first, -2 if value above last bin
 */
long AT_histo_linear_bin_no(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const double value);


/**
 * Returns bin index number for a given value for logarithmic histograms
 * @param[in] number_of_bins        number of bins in histogram
 * @param[in] lowest_left_limit     left limit of first bin
 * @param[in] step                  step between bin limits
 * @param[in] value                 value
 * @return                      bin index number (zero-based, positive), -1 if value below first, -2 if value above last bin
 */
long AT_histo_logarithmic_bin_no(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const double value);


/**
 * Returns bin index number for a given value for logarithmic histograms
 * @param[in] number_of_bins      number of bins in histogram
 * @param[in] lowest_left_limit   left limit of first bin
 * @param[in] step                step between bin limits
 * @param[in] histo_type          type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[in] value               value
 * @return                        bin index number (zero-based, positive), -1 if value below first, -2 if value above last bin
 */
long AT_histo_bin_no(      const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const double value);


/**
 * Increases bin frequency content by 'weight' for a given value
 * @param[in] number_of_bins       number of bins in histogram
 * @param[in] lowest_left_limit    left limit of first bin
 * @param[in] step                 step between bin limits
 * @param[in] histo_type           type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[in] value                value
 * @param[in] weight               weight by which the bin frequency content for 'value' is increase (usually 1)
 * @param[in,out] frequency        vector of frequencies for the histogram (array of number_of_bins)
 */
void AT_histo_add_single(      const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const double value,
    const double weight,
    double frequency[]);

/**
 * Add an array of values with individual weight to histogram
 * TODO: optimize routine, now histotype is checked unnecessarily often
 * @param[in] number_of_bins       number of bins in histogram
 * @param[in] lowest_left_limit    left limit of first bin
 * @param[in] step                 step between bin limits
 * @param[in] histo_type           type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[in] n_values             number of values given
 * @param[in] value                values (array of size number_of_values)
 * @param[in] weight               weights by which the bin frequency content for 'value' is increase (array of size number_of_values)
 * @param[in,out] frequency        vector of frequencies for the histogram (array of size number_of_bins)
 */
void AT_histo_add_multi(      const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const long n_values,
    const double value[],
    const double weight[],
    double frequency[]);

/**
 * TODO
 * @param[in] number_of_bins
 * @param[in] lowest_left_limit
 * @param[in] step
 * @param[in] histo_type
 * @param[in] frequency (array of size number_of_bins)
 * @param[out] sum
 */
void AT_histo_sum(	const long number_of_bins,
		const double lowest_left_limit,
		const double step,
		const long histo_type,
		const double frequency[],
		double* sum);


/**
 * TODO
 * @param[in] number_of_bins
 * @param[in] lowest_left_limit
 * @param[in] step
 * @param[in] histo_type
 * @param[out] frequency  (array of size number_of_bins)
 */
void AT_histo_normalize(	const long number_of_bins,
		const double lowest_left_limit,
		const double step,
		const long histo_type,
		double frequency[]);

/* TRANSIENT ROUTINES FOR TRANSFORMING OLD-STYLE KELLERER HISTOGRAMS INTO NEW STYLE */
/**
 * TODO
 * @param[in] N2
 * @return step
 */
double AT_N2_to_step( double N2 );

/**
 * TODO
 * @param[in] step
 * @return N2
 */
double AT_step_to_N2( double step );


/* OLD ROUTINES, KEPT FOR COMPATIBILITY */
/**
 * Returns bin width
 *
 * @param[in] number_of_bins               number of bin in histogram
 * @param[in] bin_centers                  bin centers (array of size number_of_bins)
 * @return bin_width
 */
double AT_histoOld_log_bin_width(	const long number_of_bins,
								const double bin_centers[]);


/**
 * Returns lower bin limit for single bin
 *
 * @param[in] number_of_bins               number of bin in histogram
 * @param[in] bin_centers                  bin centers (array of size number_of_bins)
 * @param[in] bin_no                       index of bin (zero based)
 * @return lower_left_limit
 */
double AT_histoOld_lower_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * Returns lower bin limit for single bin
 *
 * @param[in] number_of_bins               number of bin in histogram
 * @param[in] bin_centers                  bin centers (array of size number_of_bins)
 * @param[in] bin_no                       index of bin (zero based)
 * @return upper_left_limit
 */
double AT_histoOld_upper_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * Returns bin width for single bin
 *
 * @param[in] number_of_bins               number of bin in histogram
 * @param[in] bin_centers                  bin centers (array of size number_of_bins)
 * @param[in] bin_no                       index of bin (zero based)
 * @return bin width
 */
double AT_histoOld_get_bin_width(	const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * Returns bin widths
 *
 * @param[in] number_of_bins               number of bin in histogram
 * @param[in] bin_centers                  bin centers (array of size number_of_bins)
 * @param[out] bin_widths                  resulting bin widths (array of size number_of_bins)
 */
void AT_histoOld_get_bin_widths(	const long number_of_bins,
								const double bin_centers[],
								double bin_widths[]);


/**
 * Returns bin index number for a value
 *
 * @param[in] number_of_bins               number of bin in histogram
 * @param[in] bin_centers                  bin centers (array of size number_of_bins)
 * @param[in] value                        value to put into histogram
 * @return bin_no
 */
long AT_histoOld_bin_no(	const long number_of_bins,
						const double bin_centers[],
						const double value);

#endif /* AT_CONSTANTS_H_ */
