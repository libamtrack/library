/**
 * @brief Gamma Response models
 */


/*
 *    AT_GammaResponse.c
 *    ==================
 *
 *    Created on: 28.07.2009
 *    Creator: greilich
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

#include "AT_GammaResponse.h"


long AT_Gamma_index_from_material_number( const long Gamma_no ){
  long  index          =  -1;
  long  number_of_GRs  =  1;
  find_elements_int(  &Gamma_no,
      number_of_GRs,
      AT_GR_Data.GR_no,
      AT_GR_Data.n,
      &index);   // TODO replace call to find_elements_int by call to simpler function which will find the index just for one argument
  return index;
}


void AT_Gamma_name_from_number( const long Gamma_no,
    char* Gamma_name){

  long  index = AT_Gamma_index_from_material_number( Gamma_no );

  if( index != -1){
    strcpy(Gamma_name, AT_GR_Data.GR_name[index]);
  } else {
    strcpy(Gamma_name,"*** invalid choice ***");
  }
}


long AT_Gamma_number_of_parameters( const long Gamma_no){
    long  index = AT_Gamma_index_from_material_number( Gamma_no );
    if( index == -1){
#ifndef NDEBUG
      printf("gamma model no %ld not found\n", Gamma_no);
#endif
      return 0;
    }
    return AT_GR_Data.n_parameters[index];
}


void AT_gamma_response( const long  number_of_doses,
    const double   d_Gy[],
    const long     gamma_model,
    const double   gamma_parameter[],
    const bool	   lethal_events_mode,
    double         S[]){

  assert( number_of_doses > 0);
  assert( d_Gy != NULL);
  assert( gamma_parameter != NULL);

  long  i,j;
  /*
   * (0) Test model, response m*x + c
   * parameters:  m - gamma_parameter[0]
   *              c - gamma_parameter[1]
   */
  if(gamma_model == GR_Test){
    double  m  = gamma_parameter[0];
    double  c  = gamma_parameter[1];
    for (i = 0; i < number_of_doses; i++){
      if(lethal_events_mode){
    	  S[i]    =  0;								// No lethal_event_mode that would make sense
      }else{
    	  S[i]    =  m * d_Gy[i] + c;
      }
    }
  }
  /*
   *  (1) General multi-hit, multi-target model
   *
   *  N.B.: In this case the array length IS NOT (*n_parameter) but 4*(*n_parameter) !!
   *
   *  4 * i parameters:  k    - relative contribution from component i (= Smax_i), preferably adding to 1 for all components (not mandatory though!)
   *             D1    - characteristic dose (one hit per target in average) of component i
   *             c    - hittedness of component i
   *             m    - number of targets for component i
   *  if 4*ith parameter (k_i == 0) -> end of parameter list
   **/
  if(gamma_model == GR_GeneralTarget){
	  long  n_gamma_parameter = 0;
	  while  (gamma_parameter[n_gamma_parameter] != 0){
		  n_gamma_parameter  += 4;
	  }

	  for (i = 0; i < number_of_doses; i++){
		  S[i]  = 0.0;
	  }

	  for (j = 0; j < n_gamma_parameter / 4; j++){            // loop over all components
		  double k      =  gamma_parameter[j * 4];
		  double D1     =  gamma_parameter[j * 4 + 1];
		  double c      =  gamma_parameter[j * 4 + 2];
		  double m      =  gamma_parameter[j * 4 + 3];

		  assert( D1 > 0);

		  double tmp    =  0.0;

		  for (i = 0; i < number_of_doses; i++){
			  if(c == 1){                    // in case of c = 1 use simplified, faster formula
				  tmp      =  1.0 - exp(-1.0 * d_Gy[i] / D1);
			  }else{                      // for c > 1 use incomplete gamma function
				  assert( c > 0 );
				  assert( d_Gy[i] >= 0 );
				  tmp      =  gsl_sf_gamma_inc_P(c, d_Gy[i] / D1);
			  }

			  if(m == 1){                    // in case of m = 1 do not use pow for speed reasons
				  tmp      *=  k;
			  }else{
				  tmp      = k * pow(tmp, m);
			  }

			  S[i]    +=  tmp;
		  }
	  }
	  if(lethal_events_mode){
		  for(i = 0; i < number_of_doses; i++){
			  assert( S[i] < 1.0);
			  S[i]	= -1.0 * log(1.0 - S[i]);
		  }
	  }
  }
  /*
   *  (2) RL ACCUMULATED COUNTS MODEL
   *
   *  parameters:  Smax   - max. RATE signal (in contrast to model == 1, where Smax is the maximum response signal)
   *         chi    - dynamics, ration between Smax and start count rate c0
   *         D1    - characteristic dose at which rate signal reaches saturation
   *
   *         are transformed into
   *
   *         c0    - RATE signal at D = 0
   *         B    - linear slope of RATE signal below D = D1
   *         D1    - see above
   */

  if(gamma_model == GR_Radioluminescence){
	  double Smax   =  gamma_parameter[0];
	  double D1     =  gamma_parameter[1];
	  double chi    =  gamma_parameter[2];

	  assert( chi > 0);
	  assert( D1 > 0);

	  // transform parameters
	  double c0     =  Smax / chi;
	  double B      =  (Smax - c0) / D1;

	  for (i = 0; i < number_of_doses; i++){
		  if(d_Gy[i]  <= D1){
			  S[i]    =  c0 * d_Gy[i] + 0.5 * B * gsl_pow_2(d_Gy[i]);
		  }else{
			  S[i]    =  c0 * D1 + 0.5 * B * gsl_pow_2(D1) + Smax * (d_Gy[i] - D1);
		  }
	  }
	  if(lethal_events_mode){
		  for(i = 0; i < number_of_doses; i++){
			  assert( S[i] < 1.0);
			  S[i]	= -1.0 * log(1.0 - S[i]);
		  }
	  }
  }

  /*
   *  (3) EXPONENTIAL-SATURATION MODEL, obsolete
   *
   *  parameters:  Smax   - max. response signal
   *         D1    - characteristic dose at which rate signal reaches saturation
   */
  if(gamma_model == GR_ExpSaturation){
    double  Smax  =  gamma_parameter[0];
    double  D0    =  gamma_parameter[1];

    assert( D0 > 0);

    for (i = 0; i < number_of_doses; i++){
      if(lethal_events_mode){
    	  S[i]	  =  d_Gy[i] / D0;
      }else{
     	  S[i]    =  Smax * (1.0 - exp(-1.0 * d_Gy[i] / D0));
     }
    }
  }

  /*
   *  (4) LINEAR-QUADRATIC MODEL
   *
   *  parameters:  alpha   - 1st parameter in LQ equation
   *               beta    - 2nd parameter in LQ equation
   *               D0      - 3rd parameter - transition-dose
   */
  if(gamma_model == GR_LinQuad){

    double  alpha =  gamma_parameter[0];
    double  beta  =  gamma_parameter[1];
    double  D0    =  gamma_parameter[2];

    assert( alpha > 0. );
    assert( beta >= 0. );
    assert( D0 > 0. );

    for (i = 0; i < number_of_doses; i++){
		if( d_Gy[i] < D0 ){
			S[i]    =       alpha * d_Gy[i] + beta * gsl_pow_2(d_Gy[i]);
		} else {
			S[i]    =       alpha * D0 + beta * gsl_pow_2(D0) + ( alpha + 2.0 * beta * D0) * (d_Gy[i] - D0);
		}
    	if(!lethal_events_mode){  // non-lethal events mode
    		S[i]    =       exp( -S[i] );
    	}
    }
  }

  /*
   *  (6) Multi-component model as given in Geiss, 1997 (PhD thesis)
   *
   *  N.B.: This function is supposed to be a mixed of one- and two-hit characteristics. Due to missing linear terms in the second
   *        component, however, it cannot be expressed using the general hit-target model (implemented as (1))
   *
   *  5 parameters:  c  - constant, arbitrary for efficiency prediction
   *                 k1 - contribution from linear component (0.8 for Geiss, 1997)
   *                 a1 - dose prefactor for linear component (3e-4 for Geiss, 1997)
   *                 k2 - contribution from quadratic component (0.2 for Geiss, 1997)
   *                 a2 - dose prefactor for quadratic component (1e-6 for Geiss, 1997)
   **/
  if(gamma_model == GR_Geiss){
	  double co	= gamma_parameter[0];
	  double k1	= gamma_parameter[1];
	  double a1	= gamma_parameter[2];
	  double k2	= gamma_parameter[3];
	  double a2	= gamma_parameter[4];

	  for (i = 0; i < number_of_doses; i++){
		  if(lethal_events_mode){
			  // TODO: Use more efficient formula
			  assert( 1.0 - k1 * (1.0 - exp(-a1 * d_Gy[i])) - k2 * (1.0 - exp(-a2 * d_Gy[i] * d_Gy[i])) > 0. );
			  S[i]    = -log(1.0 - k1 * (1.0 - exp(-a1 * d_Gy[i])) - k2 * (1.0 - exp(-a2 * gsl_pow_2(d_Gy[i]))));
		  } else {
			  S[i]    = co * (k1 * (1.0 - exp(-a1 * d_Gy[i])) + k2 * (1.0 - exp(-a2 * gsl_pow_2(d_Gy[i]))));
		  }
	  }
  }

}


double AT_get_gamma_response_for_average_dose(  const long  number_of_bins,
    const double   dose_Gy_bin_position[],
    const double   dose_Gy_bin_width[],
    const double   dose_bin_frequency[],
    const long     gamma_model,
    const double   gamma_parameter[],
    const bool     lethal_events_mode)
{
  long i;
  double total_dose =  0.0;
  for(i = 0; i < number_of_bins; i++){
    total_dose   +=  dose_Gy_bin_position[i] * dose_Gy_bin_width[i] * dose_bin_frequency[i];
  }

  double gamma_response = 0.0;
  i  = 1;
  AT_gamma_response(  i,
      &total_dose,
      gamma_model,
      gamma_parameter,
      lethal_events_mode,
      &gamma_response);

  return gamma_response;
}


void AT_get_response_distribution_from_dose_distribution(  const long  number_of_bins,
    const double   dose_Gy_bin_position[],
    const double   dose_bin_frequency[],
    const long     gamma_model,
    const double   gamma_parameter[],
    const bool     lethal_events_mode,
    double         response_bin_frequency[])
{
  AT_gamma_response(  number_of_bins,
      dose_Gy_bin_position,
      gamma_model,
      gamma_parameter,
      lethal_events_mode,
      response_bin_frequency);
}


double AT_get_ion_response_from_response_distribution(  const long  number_of_bins,
	const double   dose_Gy_bin_width[],
    const double   dose_bin_frequency[],
    const double   ion_response_bin_frequency[])
{
  long i;
  double ion_response =  0.0;
  for(i = 0; i < number_of_bins; i++){
    ion_response    +=  ion_response_bin_frequency[i] * dose_Gy_bin_width[i] * dose_bin_frequency[i];
  }
  return ion_response;
}


double AT_get_ion_response_from_dose_distribution(  const long  number_of_bins,
	    const double   dose_Gy_bin_position[],
	    const double   dose_Gy_bin_width[],
	    const double   dose_bin_frequency[],
	    const long     gamma_model,
	    const double   gamma_parameter[],
	    const bool     lethal_events_mode)
{
	assert( number_of_bins > 0);
	double* ion_response_bin_frequency = (double*)calloc(sizeof(double),number_of_bins);

	AT_get_response_distribution_from_dose_distribution(  number_of_bins,
	    dose_Gy_bin_position,
	    dose_bin_frequency,
	    gamma_model,
	    gamma_parameter,
	    lethal_events_mode,
	    ion_response_bin_frequency);

	double ion_response = AT_get_ion_response_from_response_distribution( number_of_bins,
			  dose_Gy_bin_width,
			  dose_bin_frequency,
			  ion_response_bin_frequency);

	free( ion_response_bin_frequency);

	return ion_response;
}


double AT_get_ion_efficiency_from_dose_distribution(  const long  number_of_bins,
	    const double   dose_Gy_bin_position[],
	    const double   dose_Gy_bin_width[],
	    const double   dose_bin_frequency[],
	    const long     gamma_model,
	    const double   gamma_parameter[],
	    const bool     lethal_events_mode)
{
  double ion_response = AT_get_ion_response_from_dose_distribution(  number_of_bins,
		    dose_Gy_bin_position,
		    dose_Gy_bin_width,
		    dose_bin_frequency,
		    gamma_model,
		    gamma_parameter,
		    lethal_events_mode);

  double gamma_response = AT_get_gamma_response_for_average_dose( number_of_bins,
		  dose_Gy_bin_position,
		  dose_Gy_bin_width,
		  dose_bin_frequency,
		  gamma_model,
		  gamma_parameter,
		  lethal_events_mode);

  assert( ion_response >= 0);
  assert( gamma_response > 0);

  return ion_response / gamma_response;
}


double AT_get_ion_efficiency_from_response_distribution(  const long  number_of_bins,
	    const double   dose_Gy_bin_position[],
	    const double   dose_Gy_bin_width[],
	    const double   dose_bin_frequency[],
	    const double   ion_response_bin_frequency[],
	    const long     gamma_model,
	    const double   gamma_parameter[],
	    const bool     lethal_events_mode)
{
  double ion_response = AT_get_ion_response_from_response_distribution( number_of_bins,
		  dose_Gy_bin_width,
		  dose_bin_frequency,
		  ion_response_bin_frequency);

  double gamma_response = AT_get_gamma_response_for_average_dose( number_of_bins,
		  dose_Gy_bin_position,
		  dose_Gy_bin_width,
		  dose_bin_frequency,
		  gamma_model,
		  gamma_parameter,
		  lethal_events_mode);

  assert( ion_response >= 0);
  assert( gamma_response > 0);

  return ion_response / gamma_response;
}


void AT_get_gamma_response(  const long  number_of_bins,
    const double   d_Gy[],
    const double   dd_Gy[],
    const double   f[],
    const double   f0,
    const long     gamma_model,
    const double   gamma_parameter[],
    const bool    lethal_events_mode,
    double   S[],
    double*  S_HCP,
    double*  S_gamma,
    double*  efficiency)
{
  long i;

  AT_gamma_response(  number_of_bins,
      d_Gy,
      gamma_model,
      gamma_parameter,
      lethal_events_mode,
      S);

  *S_HCP         =  0.0;
  double D_gamma =  0.0;

  for(i = 0; i < number_of_bins; i++){
    D_gamma   +=  d_Gy[i] * dd_Gy[i] * f[i];
    *S_HCP    +=  S[i] * dd_Gy[i] * f[i];
  }

  i  = 1;
  AT_gamma_response(  i,
      &D_gamma,
      gamma_model,
      gamma_parameter,
      lethal_events_mode,
      S_gamma);

  assert( *S_HCP >= 0);
  assert( *S_gamma > 0);

  *efficiency    =  *S_HCP / *S_gamma;
}
