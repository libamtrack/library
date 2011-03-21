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
		*left_limit = lowest_left_limit * pow(step, bin_no);
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
  for (bin_no = 0; bin_no < number_of_bins; bin_no++){
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
                                          bin_no,
                                          &left_limit);
  if( status_code == EXIT_FAILURE )
	  return status_code;
  status_code = AT_histo_logarithmic_left_limit(       number_of_bins,
                                          lowest_left_limit,
                                          step,
                                          bin_no + 1,
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

///////////////////////////////// Step routines ////////////////////////////////////
int AT_histo_linear_step(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double highest_left_limit,
                                double * step)
{
  if ((number_of_bins <= 0)||(highest_left_limit <= lowest_left_limit)){
    return EXIT_FAILURE;
  }

  * step = (highest_left_limit - lowest_left_limit) / number_of_bins;
  return EXIT_SUCCESS;
}


int AT_histo_logarithmic_step(      const long number_of_bins,
                                const double lowest_left_limit,
                                const double highest_left_limit,
                                double * step)
{
  if ((number_of_bins <= 0)||(highest_left_limit <= lowest_left_limit)
     ||(lowest_left_limit <= 0)||(highest_left_limit <= 0)){
    return EXIT_FAILURE;
  }

  * step = pow(highest_left_limit / lowest_left_limit, 1.0 / (double)number_of_bins);
  return EXIT_SUCCESS;
}


int AT_histo_step(      const long number_of_bins,
                        const double lowest_left_limit,
                        const double highest_left_limit,
                        const long histo_type,
                        double * step)
{
  int status_code = EXIT_FAILURE;
  if (histo_type == AT_histo_linear){
          status_code = AT_histo_linear_step( number_of_bins,
              lowest_left_limit,
              highest_left_limit,
              step);
  }else{
          status_code = AT_histo_logarithmic_step( number_of_bins,
              lowest_left_limit,
              highest_left_limit,
              step);
  }
  return status_code;
}

///////////////////////////////// Number of bin routines ////////////////////////////////////

int AT_histo_linear_n_bins(     const double lowest_left_limit,
                                const double highest_left_limit,
                                const double step,
                                long * number_of_bins)
{
	if ((step <= 0)||(highest_left_limit <= lowest_left_limit)){
		return EXIT_FAILURE;
	}

	*number_of_bins = (highest_left_limit - lowest_left_limit) / step;
	return EXIT_SUCCESS;
}


int AT_histo_logarithmic_n_bins(     const double lowest_left_limit,
                                const double highest_left_limit,
                                const double step,
                                long * number_of_bins)
{
	  if ((step <= 1.0)||(highest_left_limit <= lowest_left_limit)
	     ||(lowest_left_limit <= 0)||(highest_left_limit <= 0)){
	    return EXIT_FAILURE;
	  }

	  *number_of_bins = floor(log(highest_left_limit / lowest_left_limit) / log(step)) + 1;
	  return EXIT_SUCCESS;
}


int AT_histo_n_bins(     const double lowest_left_limit,
    const double highest_left_limit,
    const double step,
    const long histo_type,
    long * number_of_bins)
{
	  int status_code = EXIT_FAILURE;
	  if (histo_type == AT_histo_linear){
	          status_code = AT_histo_linear_n_bins( lowest_left_limit,
	              highest_left_limit,
	              step,
	              number_of_bins);
	  }else{
	          status_code = AT_histo_logarithmic_n_bins( lowest_left_limit,
	              highest_left_limit,
	              step,
	              number_of_bins);
	  }
	  return status_code;
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

  assert( step > 0.0 );
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
  long bin_no = floor(log10(factor) / log10(step));

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


void AT_histo_add_single(      const long number_of_bins,
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
  if(bin_no >= 0){
    frequency[bin_no] += weight;
  }
}

void AT_histo_add_multi(      const long number_of_bins,
    const double lowest_left_limit,
    const double step,
    const long histo_type,
    const long n_values,
    const double value[],
    const double weight[],
    double frequency[])
{
	long i;
	for(i = 0; i< n_values; i++){
		AT_histo_add_single(number_of_bins,
				lowest_left_limit,
				step,
				histo_type,
				value[i],
				weight[i],
				frequency);
	}
}

//TODO rewrite in such way that it has return type "double"
void AT_histo_sum(	const long number_of_bins,
		const double lowest_left_limit,
		const double step,
		const long histo_type,
		const double frequency[],
		double* sum)
{
	long 	i;
	double midpoint = 0.0;
	double bin_width = 0.0;
	*sum	= 0.0;
	for(i = 0; i < number_of_bins; i++){
		AT_histo_midpoint( number_of_bins,
				lowest_left_limit,
				step,
				histo_type,
				i,
				&midpoint);
		AT_histo_bin_width( number_of_bins,
				lowest_left_limit,
				step,
				histo_type,
				i,
				&bin_width);
		*sum	+= frequency[i] * midpoint * bin_width;
	}
}


void AT_histo_normalize(	const long number_of_bins,
		const double lowest_left_limit,
		const double step,
		const long histo_type,
		double frequency[])
{
	/* get histogram sum */
	double 	sum;
	AT_histo_sum(	number_of_bins,
			lowest_left_limit,
			step,
			histo_type,
			frequency,
			&sum);

	/* divide frequencies by sum */
	assert(sum > 0.0);
	long 	i;
	for(i = 0; i < number_of_bins; i++){
		frequency[i] /= sum;
	}
}

/* TRANSIENT ROUTINES FOR TRANSFORMING OLD-STYLE KELLERER HISTOGRAMS INTO NEW STYLE */
double AT_N2_to_step(double N2){
	return(pow(2.0, 1.0/(double)(N2)));
}

double AT_step_to_N2(double step){
	return(log(2.0)/log(step));
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
