#ifndef AT_HISTOGRAMS_H_
#define AT_HISTOGRAMS_H_

/**
 * @file
 * @brief This files handels histograms used by libamtrack
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

#ifndef NAN
#define NAN (0.0/0.0)
#endif

/**
 * @enum AT_histo_type
 * Histogram type enumerator
 */
enum AT_histo_type{
  AT_histo_linear       = 0, /** Histogram with linear bins and arithmetic midpoints, step is the bin limit difference **/
  AT_histo_log          = 1  /** Histogram with logarithmic bins and geometrical midpoints, step ist the bin limit factor **/
};


///////////////////////////////// Left limit routines ////////////////////////////////////
/**
 * Returns left limit of bin with given index number for linear histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @return                      left limit of bin
 */
double AT_histo_linear_left_limit(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no);

/**
 * Returns left limit of bin with given index number for logarithmic histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @return                      left limit of bin
 */
double AT_histo_logarithmic_left_limit(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no);


/**
 * Returns left limit of bin with given index number
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @param histo_type            type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @return                      left limit of bin
 */
double AT_histo_left_limit( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const long bin_no);


/**
 * Returns vector with all left limits of a histogram
 * @param[in] number_of_bins            number of bins in histogram
 * @param[in] lowest_left_limit         left limit of first bin
 * @param[in] step                      step between bin limits
 * @param[in] histo_type                type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] left_limits              array of bin width
 * @length number_of_bins + 1
 */
void AT_histo_left_limits( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double left_limits[]);


///////////////////////////////// Bin widths routines ////////////////////////////////////
/**
 * Returns width of bin with given index number for linear histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @return                      width of bin
 */
double AT_histo_linear_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no);


/**
 * Returns width of bin with given index number for logarithmic histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @return                      width of bin
 */
double AT_histo_logarithmic_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no);


/**
 * Returns width of bin with given index number
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @param histo_type            type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @return                      width of bin
 */
double AT_histo_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long histo_type,
                                const long bin_no);


/**
 * Returns vector with all bin widths of a histogram
 * @param[in] number_of_bins            number of bins in histogram
 * @param[in] lowest_left_limit         left limit of first bin
 * @param[in] step                      step between bin limits
 * @param[in] histo_type                type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] bin_widths              array of bin width
 * @length number_of_bins
 */
void AT_histo_bin_widths( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double bin_widths[]);


///////////////////////////////// Midpoint routines ////////////////////////////////////
/**
 * Returns midpoint of bin with given index number for linear histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @return                      midpoint of bin
 */
double AT_histo_linear_midpoint(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no);


/**
 * Returns midpoint of bin with given index number for logarithmic histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @return                      midpoint of bin
 */
double AT_histo_logarithmic_midpoint(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no);


/**
 * Returns midpoint of bin with given index number
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param bin_no                index number of bin (zero-based)
 * @param histo_type            type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @return                      midpoint of bin
 */
double AT_histo_midpoint(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long histo_type,
                                const long bin_no);


/**
 * Returns vector with all midpoints of a histogram
 * @param[in] number_of_bins            number of bins in histogram
 * @param[in] lowest_left_limit         left limit of first bin
 * @param[in] step                      step between bin limits
 * @param[in] histo_type                type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param[out] bin_widths               array of midpoints
 * @length number_of_bins
 */
void AT_histo_midpoints( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double midpoints[]);


///////////////////////////////// Access routines ////////////////////////////////////
/**
 * Returns bin index number for a given value for linear histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param value                 value
 * @return                      bin index number (zero-based, positive), -1 if value below first, -2 if value above last bin
 */
long AT_histo_linear_bin_no(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const double value);


/**
 * Returns bin index number for a given value for logarithmic histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param value                 value
 * @return                      bin index number (zero-based, positive), -1 if value below first, -2 if value above last bin
 */
long AT_histo_logarithmic_bin_no(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const double value);


/**
 * Returns bin index number for a given value for logarithmic histograms
 * @param number_of_bins        number of bins in histogram
 * @param left_limit            left limit of first bin
 * @param step                  step between bin limits
 * @param[in] histo_type        type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param value                 value
 * @return                      bin index number (zero-based, positive), -1 if value below first, -2 if value above last bin
 */
long AT_histo_bin_no(      const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const double value);


/**
 * Increases bin frequency content by 'weight' for a given value
 * @param[in] number_of_bins    number of bins in histogram
 * @param[in] left_limit        left limit of first bin
 * @param[in] step              step between bin limits
 * @param[in] histo_type        type of histogram (linear or logarithmic)
 * @see AT_histo_type
 * @param value[in]             value
 * @param weight[in]            weight by which the bin frequency content for 'value' is increase (usually 1)
 * @param frequency[in/out]     vector of frequencies for the histogram
 * @length number_of_bins
 */
void AT_histo_add(      const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const double value,
    const double weight,
    double frequency[]);


/* OLD ROUTINES, KEPT FOR COMPATIBILITY */
/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @return
 */
double AT_histoOld_log_bin_width(	const long number_of_bins,
								const double bin_centers[]);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param bin_no
 * @return
 */
double AT_histoOld_lower_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param bin_no
 * @return
 */
double AT_histoOld_upper_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param bin_no
 * @return
 */
double AT_histoOld_get_bin_width(	const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param bin_widths
 */
void AT_histoOld_get_bin_widths(	const long number_of_bins,
								const double bin_centers[],
								double bin_widths[]);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param value
 * @return
 */
long AT_histoOld_bin_no(	const long number_of_bins,
						const double bin_centers[],
						const double value);

#endif /* AT_CONSTANTS_H_ */
