/**
 * @file
 * @brief Gamma Response models
 */


/*
*    AT_GammaResponse.c
*    ==================
*
*    Created on: 28.07.2009
*    Author: greilich
*
*    Copyright 2006, 2009 Steffen Greilich / the libamtrack team
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


void getGammaName( const long Gamma_no,
    char* Gamma_name){

  // find look-up index for GR number in GR data table
  long  match;
  const long n_tmp = 1;

  find_elements_int(  &Gamma_no,
      n_tmp,
      AT_GR_Data.GR_no,
      AT_GR_Data.n,
      &match);

  if( match != -1){
    strcpy(Gamma_name, AT_GR_Data.GR_name[match]);
  } else {
    strcpy(Gamma_name,"*** invalid choice ***");
  }
}


void AT_gamma_response( const long  number_of_doses,
    const float   d_Gy[],
    const long    gamma_model,
    const float   gamma_parameter[],
    // return
    float         S[]){

  long  i,j;
  /*
   * (0) Testmodel, response m*x + c
   * parameters:  m - gamma_parameter[0]
   *         c - gamma_parameter[1]
   */
  if(gamma_model == GR_Test){
    float  m  = gamma_parameter[0];
    float  c  = gamma_parameter[1];
    for (i = 0; i < number_of_doses; i++){
      S[i]    =  m * d_Gy[i] + c;
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
      S[i]  = 0.0f;
    }
    for (j = 0; j < n_gamma_parameter / 4; j++){            // loop over all components
      float k      =  gamma_parameter[j * 4];
      float D1     =  gamma_parameter[j * 4 + 1];
      float c      =  gamma_parameter[j * 4 + 2];
      float m      =  gamma_parameter[j * 4 + 3];

      float tmp    =  0.0f;

      for (i = 0; i < number_of_doses; i++){
        if(c == 1){                    // in case of c = 1 use simplified, faster formula
          tmp      =  1.0f - exp(-1.0f * d_Gy[i] / D1);
        }else{                      // for c > 1 use incomplete gamma function
          tmp      =  gammp(c, d_Gy[i] / D1);
        }

        if(m == 1){                    // in case of m = 1 do not use pow for speed reasons
          tmp      *=  k;
        }else{
          tmp      = k * pow(tmp, m);
        }

        S[i]    +=  tmp;
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
    float Smax   =  gamma_parameter[0];
    float D1     =  gamma_parameter[1];
    float chi    =  gamma_parameter[2];
    // transform parameters
    float c0     =  Smax / chi;
    float B      =  (Smax - c0) / D1;

    for (i = 0; i < number_of_doses; i++){
      if(d_Gy[i]  <= D1){
        S[i]    =  c0 * d_Gy[i] + 0.5f * B * gsl_pow_2(d_Gy[i]);
      }else{
        S[i]    =  c0 * D1 + 0.5f * B * gsl_pow_2(D1) + Smax * (d_Gy[i] - D1);
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
    // exp-saturation model
    float  Smax  =  gamma_parameter[0];
    float  D0    =  gamma_parameter[1];

    for (i = 0; i < number_of_doses; i++){
      S[i]    =  Smax * (1.0f - (float)exp(-1.0f * d_Gy[i] / D0));
    }
  }

  /*
   *  (4) LINEAR-QUADRATIC MODEL
   *
   *  parameters:  alpha   - 1st parameter in equation SF = exp( - alpha *D^2 - beta *D)
   *               beta    - 2nd parameter in equation SF = exp( - alpha *D^2 - beta *D)
   *               D0      - 3rd parameter - transition-dose
   */
  if(gamma_model == GR_LinQuad){

    float  alpha =  gamma_parameter[0];
    float  beta  =  gamma_parameter[1];
    float  D0    =  gamma_parameter[2];

    if( alpha < 0.0f )
      alpha = 0.0f;
    if( beta < 0.0f )
      beta = 0.0f;

    for (i = 0; i < number_of_doses; i++){
      if( d_Gy[i] < D0 ){
        S[i]    =  expf( -alpha * d_Gy[i] - beta * gsl_pow_2(d_Gy[i]));
      } else {
        S[i]    =  expf( -alpha * D0 - beta * gsl_pow_2(D0) - ( alpha + 2.0f * beta * D0) * (d_Gy[i] - D0) );
      }
    }
  }

  /*
   *  (4) LINEAR-QUADRATIC MODEL
   *
   *      parameters:     alpha   - 1st parameter in equation SF = exp( - alpha *D^2 - beta *D)
   *                      beta    - 2nd parameter in equation SF = exp( - alpha *D^2 - beta *D)
   *                      D0      - 3rd parameter - transition-dose
   */

  if(gamma_model == GR_LinQuad_Log){

    float   alpha   =       gamma_parameter[0];
    float   beta    =       gamma_parameter[1];
    float   D0      =       gamma_parameter[2];

    if( alpha < 0.0f )
      alpha = 0.0f;
    if( beta < 0.0f )
      beta = 0.0f;

    for (i = 0; i < number_of_doses; i++){
      if( d_Gy[i] < D0 ){
        S[i]    =       alpha * d_Gy[i] + beta * d_Gy[i] * d_Gy[i];
      } else {
        S[i]    =       alpha * D0 + beta * gsl_pow_2(D0) + ( alpha + 2.0f * beta * D0) * (d_Gy[i] - D0);
      }
    }
  }
}


void AT_get_gamma_response(  const long  number_of_bins,
    const float   d_Gy[],
    const float   dd_Gy[],
    const float   f[],
    const float   f0,
    const long    gamma_model,
    const float   gamma_parameter[],
    const bool    lethal_events_mode,
    // return
    float   S[],
    float*  S_HCP,
    float*  S_gamma,
    float*  efficiency)
{
  long i;

  AT_gamma_response(  number_of_bins,
      d_Gy,
      gamma_model,
      gamma_parameter,
      // return
      S);

  *S_HCP        =  0.0f;
  float D_gamma =  0.0f;

  for(i = 0; i < number_of_bins; i++){
    D_gamma   +=  d_Gy[i] * dd_Gy[i] * f[i];
    *S_HCP    +=  S[i] * dd_Gy[i] * f[i];
  }

  if( lethal_events_mode ){
    *S_HCP = expf( -(*S_HCP) );
  }

  i  = 1;
  AT_gamma_response(  i,
      &D_gamma,
      gamma_model,
      gamma_parameter,
      // return
      S_gamma);

  if( lethal_events_mode ){
    *S_gamma = expf( -(*S_gamma) );
  }

  *efficiency    =  *S_HCP / *S_gamma;
}
