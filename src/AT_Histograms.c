/**
 * @file
 * @brief This files handles histograms used by libamtrack
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
#include "AT_Histograms.h"

double AT_histo_log_bin_width(	const long number_of_bins,
								const double bin_centers[])
{
	return(log(bin_centers[number_of_bins-1]) - log(bin_centers[0]))/(number_of_bins-1.);
}

double AT_histo_lower_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no)
{
	return(exp(log(bin_centers[bin_no]) - 0.5 * AT_histo_log_bin_width(	number_of_bins,
																		bin_centers)));
}

double AT_histo_upper_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no)
{
	return(exp(log(bin_centers[bin_no]) + 0.5 * AT_histo_log_bin_width(	number_of_bins,
																		bin_centers)));
}

double AT_histo_get_bin_width(	const long number_of_bins,
								const double bin_centers[],
								const long bin_no)
{
	double lower_limit		=	AT_histo_lower_bin_limit(	number_of_bins,
															bin_centers,
															bin_no);
	double upper_limit		=	AT_histo_upper_bin_limit(	number_of_bins,
															bin_centers,
															bin_no);
	return(upper_limit - lower_limit);
}

void AT_histo_get_bin_widths(	const long number_of_bins,
								const double bin_centers[],
								double bin_widths[])
{
	long i;
	for (i = 0; i < number_of_bins; i++){
		bin_widths[i] = AT_histo_get_bin_width(	number_of_bins,
												bin_centers,
												i);
	}
}

long AT_histo_bin_no(	const long number_of_bins,
						const double bin_centers[],
						const double value)
{
	double lower_limit		=	AT_histo_lower_bin_limit(	number_of_bins,
															bin_centers,
															0);
	double log_bin_width	=	AT_histo_log_bin_width(	number_of_bins,
														bin_centers);

	return(floor(  (log(value) - log(lower_limit)) / log_bin_width ));
}
