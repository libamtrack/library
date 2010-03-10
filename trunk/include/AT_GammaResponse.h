#ifndef AT_GAMMARESPONSE_H_
#define AT_GAMMARESPONSE_H_

/**
 * @file
 * @brief Gamma response models
 */

/*
*    AT_GammaResponse.h
*    ==================
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "AT_Constants.h"
#include "AT_NumericalRoutines.h"

/**
 * Gamma Response code numbers
 */
enum GammaResponseModels{
  GR_Test                  = 1,      /**< no parameters */
      GR_GeneralTarget     = 2,      /**< TODO */
      GR_Radioluminescence = 3,      /**< 0 - Smax, 1 - D0, 2 - dyn */
      GR_ExpSaturation     = 4,      /**< 0 - Smax, 1 - D0 */
      GR_LinQuad           = 5,      /**< 0 - alpha, 1 - beta, 2 - D0 */
      GR_LinQuad_Log       = 6       /**< 0 - alpha, 1 - beta, 2 - D0 */
};


#define GR_DATA_N    6


typedef struct {
  long    n;
  long    GR_no[GR_DATA_N];
  long    n_parameters[GR_DATA_N];
  char*   parameter_name[GR_DATA_N][4];
  float   parameter_default[GR_DATA_N][4];
  char*   GR_name[GR_DATA_N];
} gr_data;


static const gr_data AT_GR_Data = {
    GR_DATA_N,
    {  GR_Test,          GR_GeneralTarget,          GR_Radioluminescence,        GR_ExpSaturation, GR_LinQuad, GR_LinQuad_Log},
    {  0, 4, 3, 2, 3, 3},
    {  {"","","",""},{"S_max", "D0_Gy", "c", "m"},{"S_max","D0_Gy","dyn",""},{"S_max","D0_Gy","",""},{"alpha","beta","D0_Gy",""},{"alpha","beta","D0_Gy",""}},
    {  {0,0,0,0}, {1, 10, 1, 1}, {1,10,5,0}, {1,10,0,0}, {1, 1, 10, 0}, {1, 1, 10, 0}},
    {  "simple test gamma response",  "generalized multi-target/multi-hit gamma response",  "radioluminescence gamma response",    "exp.-sat. gamma response (obsolete, use gen. target/hit instead)", "linear-quadratic gamma response","lethal events number response"}
};


/**
 * Returns name of the gamma response model from index
 *
 * PUBLIC
 *
 * @param[in]  Gamma_no    gamma response model index
 * @param[out] Gamma_name  string containing gamma response model name
 */
void getGammaName(  const long* Gamma_no,
    char* Gamma_name);


/**
 * Returns name of the response model from index
 *
 * PUBLIC
 *
 * TODO this method is not implemented yet !
 *
 * @param[in]  Method_no    response model index
 * @param[out] Method_name  string containing response model name
 */
void getMethodName( const long* Method_no,
    char* Method_name);


/**
 * Returns the detector / cell gamma response for a vector of given doses
 * according to the chosen gamma response model
 *
 * PUBLIC
 *
 * @param[in]  n                number of doses given in vector d_Gy
 * @param[in]  d_Gy             doses in Gy (vector of length n)
 * @param[in]  gamma_model      gamma response model index
 * @param[in]  gamma_parameter  vector holding necessary parameters for the chose gamma response model
 * @param[out] S                gamma responses (vector of length n)
 */
void AT_gamma_response( const long*  n,
    const float*  d_Gy,
    const long*  gamma_model,
    const float*  gamma_parameter,
    // return
    float*  S);

/**
 * Returns the detector / cell gamma response for a local dose distribution
 * according to the chosen gamma response model, used by reponse model
 * routines in AmTrack.c
 *
 * PRIVATE
 *
 * @param[in]  n                   number of bin in given local dose distribution
 * @param[in]  d_Gy                local dose bin position in Gy (vector of length n)
 * @param[in]  dd_Gy               local dose bin width in Gy (vector of length n)
 * @param[in]  d_Gy                local dose frequency (vector of length n)
 * @param[in]  d_Gy                frequency of zero local dose (pointer to float)
 * @param[in]  gamma_model         gamma response model index
 * @param[in]  gamma_parameter     vector holding necessary parameters for the chose gamma response model
 * @param[in]  lethal_events_mode  if true, allows to do calculations for cell survival
 * @see  AmTrack.c/AT_IGK
 * @param[in]  S                   gamma responses for given bins (vector of length n)
 * @param[in]  S_HCP               HCP response for given local dose distribution (expectation value of S distribution)
 * @param[in]  S_gamma             gamma response for given local dose distribution (gamma response of expectation value of d distribution)
 * @param[in]  efficiency          RE = S_HCP/S_gamma for given local dose distribution
 */
void AT_get_gamma_response(  const long*  n,
    const float*  d_Gy,
    const float*  dd_Gy,
    const float*  f,
    const float*  f0,
    const long*  gamma_model,
    const float*  gamma_parameter,
    const bool* lethal_events_mode,
    // return
    float*  S,
    float*  S_HCP,
    float*  S_gamma,
    float*  efficiency);

#endif // AT_GAMMARESPONSE_H_
