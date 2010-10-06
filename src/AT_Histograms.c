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


///////////////////////////////// Left limit routines ////////////////////////////////////
int AT_histo_linear_left_limit( const long number_of_bins,
								const double lowest_left_limit,
								const double step,
								const long bin_no,
								double * left_limit){
	if((bin_no >= 0) && (bin_no <= number_of_bins)){
		*left_limit = lowest_left_limit + bin_no * step;
		return EXIT_SUCCESS;
	}
	return EXIT_FAILURE;
}


int AT_histo_logarithmic_left_limit(    const long number_of_bins,
										const double lowest_left_limit,
										const double step,
										const long bin_no,
										double * left_limit){
	if((bin_no >= 0) && (bin_no <= number_of_bins)){
		*left_limit = lowest_left_limit * (bin_no + 1.0) * step;
		return EXIT_SUCCESS;
	}
	return EXIT_FAILURE;
}


int AT_histo_left_limit( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const long bin_no,
    double * left_limit)
{
  int status_code;
  if (histo_type == AT_histo_linear){
    status_code = AT_histo_linear_left_limit( number_of_bins,
                                   lowest_left_limit,
                                   step,
                                   bin_no,
                                   left_limit);
  }else{
	status_code = AT_histo_logarithmic_left_limit( number_of_bins,
                                   lowest_left_limit,
                                   step,
                                   bin_no,
                                   left_limit);
  }
  return EXIT_FAILURE;
}


int AT_histo_left_limits( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double left_limits[])
{
  long bin_no;
  int global_status_code = EXIT_SUCCESS;
  for (bin_no = 0; bin_no < number_of_bins + 1; bin_no++){
	  int status_code = AT_histo_left_limit( number_of_bins,
                                          lowest_left_limit,
                                          step,
                                          histo_type,
                                          bin_no,
                                          &left_limits[bin_no]);
	  if( status_code == EXIT_FAILURE)
		  global_status_code = EXIT_FAILURE;
  }
  return global_status_code;
}

///////////////////////////////// Bin widths routines ////////////////////////////////////
double AT_histo_linear_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no)
{
  return(step);
}


int AT_histo_logarithmic_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no,
                                double * bin_width)
{
  double left_limit, right_limit;
  int status_code = AT_histo_logarithmic_left_limit(       number_of_bins,
                                          lowest_left_limit,
                                          step,
                                          bin_no + 1,
                                          &left_limit);
  if( status_code == EXIT_FAILURE )
	  return status_code;
  status_code = AT_histo_logarithmic_left_limit(       number_of_bins,
                                          lowest_left_limit,
                                          step,
                                          bin_no,
                                          &right_limit);
  if( status_code == EXIT_FAILURE )
	  return status_code;
  *bin_width = right_limit - left_limit;
  return EXIT_SUCCESS;

}


int AT_histo_bin_width(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long histo_type,
                                const long bin_no,
                                double * bin_width){
  int status_code = EXIT_FAILURE;
  if (histo_type == AT_histo_linear){
    AT_histo_linear_bin_width( number_of_bins,
                                   lowest_left_limit,
                                   step,
                                   bin_no);
    status_code = EXIT_SUCCESS;
  }else{
    status_code = AT_histo_logarithmic_bin_width( number_of_bins,
                                   lowest_left_limit,
                                   step,
                                   bin_no,
                                   bin_width);
  }
  return status_code;
}


int AT_histo_bin_widths( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double bin_widths[])
{
  long bin_no;
  int global_status_code = EXIT_SUCCESS;
  for (bin_no = 0; bin_no < number_of_bins; bin_no++){
     int status_code = AT_histo_bin_width( number_of_bins,
                                          lowest_left_limit,
                                          step,
                                          histo_type,
                                          bin_no,
                                          &bin_widths[bin_no]);
	  if( status_code == EXIT_FAILURE)
		  global_status_code = EXIT_FAILURE;
  }
  return global_status_code;
}


///////////////////////////////// Midpoint routines ////////////////////////////////////
int AT_histo_linear_midpoint(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no,
                                double * midpoint)
{
  double left_limit;
  int status_code = AT_histo_linear_left_limit( number_of_bins,
                                        lowest_left_limit,
                                        step,
                                        bin_no,
                                        &left_limit);
  if( status_code == EXIT_FAILURE)
	  return status_code;

  *midpoint = left_limit + 0.5 * step;
  return EXIT_SUCCESS;
}


int AT_histo_logarithmic_midpoint(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const long bin_no,
                                double * midpoint)
{
  double left_limit, right_limit;

  int status_code = AT_histo_logarithmic_left_limit(       number_of_bins,
          lowest_left_limit,
          step,
          bin_no + 1,
          &right_limit);

  if( status_code == EXIT_FAILURE )
	  return status_code;

  status_code = AT_histo_logarithmic_left_limit(       number_of_bins,
                                            lowest_left_limit,
                                            step,
                                            bin_no,
                                            &left_limit);

  if( status_code == EXIT_FAILURE )
	  return status_code;

  *midpoint = sqrt( left_limit * right_limit );

  return EXIT_SUCCESS;
}


int AT_histo_midpoint(      const long number_of_bins,
						const double lowest_left_limit,
						const double step,
						const long histo_type,
						const long bin_no,
						double * midpoint){
	int status_code = EXIT_FAILURE;
	if (histo_type == AT_histo_linear){
		status_code = AT_histo_linear_midpoint( number_of_bins,
				lowest_left_limit,
				step,
				bin_no,
				midpoint);
	}else{
		status_code = AT_histo_logarithmic_midpoint( number_of_bins,
				lowest_left_limit,
				step,
				bin_no,
				midpoint);
	}
	return status_code;
}


int AT_histo_midpoints( const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    double midpoints[])
{
  long bin_no;
  int global_status_code = EXIT_SUCCESS;
  for (bin_no = 0; bin_no < number_of_bins; bin_no++){
	  int status_code = AT_histo_midpoint( number_of_bins,
                                          lowest_left_limit,
                                          step,
                                          histo_type,
                                          bin_no,
                                          &midpoints[bin_no]);
	  if( status_code == EXIT_FAILURE)
		  global_status_code = EXIT_FAILURE;
  }
  return global_status_code;
}

///////////////////////////////// Access routines ////////////////////////////////////
long AT_histo_linear_bin_no(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const double value)
{
  double difference = value - lowest_left_limit;
  if (difference < 0){
    return (-1);
  }
  long bin_no = floor(difference / step);

  if (bin_no > number_of_bins - 1){
    return(-2);
  }

  return(bin_no);
}


long AT_histo_logarithmic_bin_no(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double step,
                                const double value)
{
  double factor = value / lowest_left_limit;
  if (factor < 1.0){
    return (-1);
  }
  long bin_no = floor(log10(factor) - log10(step));

  if (bin_no > number_of_bins - 1){
    return(-2);
  }

  return(bin_no);
}


long AT_histo_bin_no(      const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const double value){
  if (histo_type == AT_histo_linear){
    return(AT_histo_linear_bin_no( number_of_bins,
                                   lowest_left_limit,
                                   step,
                                   value));
  }else{
    return(AT_histo_logarithmic_bin_no( number_of_bins,
                                   lowest_left_limit,
                                   step,
                                   value));
  }
}


void AT_histo_add(      const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const double value,
    const double weight,
    double frequency[])
{
  long bin_no;
  if (histo_type == AT_histo_linear){
    bin_no = AT_histo_linear_bin_no( number_of_bins,
                                   lowest_left_limit,
                                   step,
                                   value);
  }else{
    bin_no = AT_histo_logarithmic_bin_no( number_of_bins,
                                   lowest_left_limit,
                                   step,
                                   value);
  }
  if(bin_no > 0){
    frequency[bin_no] += weight;
  }
}


/* OLD ROUTINES, KEPT FOR COMPATIBILITY */
double AT_histoOld_log_bin_width(	const long number_of_bins,
								const double bin_centers[])
{
	return(log(bin_centers[number_of_bins-1]) - log(bin_centers[0]))/(number_of_bins-1.);
}


double AT_histoOld_lower_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no)
{
	return(exp(log(bin_centers[bin_no]) - 0.5 * AT_histoOld_log_bin_width(	number_of_bins,
																		bin_centers)));
}


double AT_histoOld_upper_bin_limit(const long number_of_bins,
								const double bin_centers[],
								const long bin_no)
{
	return(exp(log(bin_centers[bin_no]) + 0.5 * AT_histoOld_log_bin_width(	number_of_bins,
																		bin_centers)));
}


double AT_histoOld_get_bin_width(	const long number_of_bins,
								const double bin_centers[],
								const long bin_no)
{
	double lower_limit		=	AT_histoOld_lower_bin_limit(	number_of_bins,
															bin_centers,
															bin_no);
	double upper_limit		=	AT_histoOld_upper_bin_limit(	number_of_bins,
															bin_centers,
															bin_no);
	return(upper_limit - lower_limit);
}


void AT_histoOld_get_bin_widths(	const long number_of_bins,
								const double bin_centers[],
								double bin_widths[])
{
	long i;
	for (i = 0; i < number_of_bins; i++){
		bin_widths[i] = AT_histoOld_get_bin_width(	number_of_bins,
												bin_centers,
												i);
	}
}


long AT_histoOld_bin_no(	const long number_of_bins,
						const double bin_centers[],
						const double value)
{
	double lower_limit		=	AT_histoOld_lower_bin_limit(	number_of_bins,
															bin_centers,
															0);
	double log_bin_width	=	AT_histoOld_log_bin_width(	number_of_bins,
														bin_centers);

	return(floor(  (log(value) - log(lower_limit)) / log_bin_width ));
}
