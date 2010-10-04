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

double AT_histo_log_bin_width(	const long number_of_bins,
								const double bin_centers[]);

double AT_histo_lower_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no);

double AT_histo_upper_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no);

double AT_histo_get_bin_width(	const long number_of_bins,
								const double bin_centers[],
								const long bin_no);

void AT_histo_get_bin_widths(	const long number_of_bins,
								const double bin_centers[],
								double bin_widths[]);

long AT_histo_bin_no(	const long number_of_bins,
						const double bin_centers[],
						const double value);

#endif /* AT_CONSTANTS_H_ */
