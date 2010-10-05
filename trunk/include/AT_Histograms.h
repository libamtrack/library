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

/**
 * @enum AT_histo_type
 * Histogram type enumerator
 */
enum AT_histo_type{
  AT_histo_linear       = 0, /** Histogram with linear bins and arithmetic mid-points **/
  AT_hiso_log           = 1  /** Histogram with logarithmic bins and geometrical mid-points **/
};


/**
 * TODO
 * @param number_of_bins
 * @param left_limits
 * @param bin_widths
 */
void AT_histo_bin_width(        const long number_of_bins,
                                const double left_limits[],
                                double bin_widths[]);


/**
 * TODO
 * @param number_of_bins
 * @param left_limits
 * @param histo_type
 * @param midpoints
 */
void AT_histo_midpoints(        const long number_of_bins,
                                const double left_limits[],
                                const long histo_type,
                                double midpoints[]);


/**
 * TODO
 * @param number_of_bins
 * @param left_limits
 * @param midpoints
 */
void AT_histo_arithmetic_midpoints( const long number_of_bins,
    const double left_limits[],
    double midpoints[]);


/**
 * TODO
 * @param number_of_bins
 * @param left_limits
 * @param midpoints
 */
void AT_histo_geometric_midpoints( const long number_of_bins,
    const double left_limits[],
    double midpoints[]);


/**
 * TODO
 * @param number_of_bins
 * @param midpoints
 * @param histo_type
 * @param left_limits
 */
void AT_histo_left_limits( const long number_of_bins,
    const double midpoints[],
    const long histo_type,
    double left_limits[]);


/**
 * TODO
 * @param number_of_bins
 * @param midpoints
 * @param left_limits
 */
void AT_histo_arithmetic_left_limits( const long number_of_bins,
    const double midpoints[],
    double left_limits[]);


/**
 * TODO
 * @param number_of_bins
 * @param midpoints
 * @param left_limits
 */
void AT_histo_geometric_left_limits( const long number_of_bins,
    const double midpoints[],
    double left_limits[]);

/* OLD ROUTINES, KEPT FOR COMPATIBILITY */
/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @return
 */
double AT_histo_log_bin_width(	const long number_of_bins,
								const double bin_centers[]);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param bin_no
 * @return
 */
double AT_histo_lower_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param bin_no
 * @return
 */
double AT_histo_upper_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param bin_no
 * @return
 */
double AT_histo_get_bin_width(	const long number_of_bins,
								const double bin_centers[],
								const long bin_no);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param bin_widths
 */
void AT_histo_get_bin_widths(	const long number_of_bins,
								const double bin_centers[],
								double bin_widths[]);


/**
 * TODO
 * @param number_of_bins
 * @param bin_centers
 * @param value
 * @return
 */
long AT_histo_bin_no(	const long number_of_bins,
						const double bin_centers[],
						const double value);

#endif /* AT_CONSTANTS_H_ */
