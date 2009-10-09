#ifndef AT_GAMMARESPONSE_H_
#define AT_GAMMARESPONSE_H_

/**
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
#include <stdbool.h>

#include "AT_Utils.h"

extern int indent_counter;
extern char isp[];
extern FILE * debf;


void AT_gamma_response(  long*  n,
              float*  d_Gy,
              long*  gamma_model,
              float*  gamma_parameter,
              // return
              float*  S);

void AT_get_gamma_response(  long*  n,
                float*  d_Gy,
                float*  dd_Gy,
                float*  f,
                float*  f0,
                long*  gamma_model,
                float*  gamma_parameter,
                bool* lethal_events_mode,
                // return
                float*  S,
                float*  S_HCP,
                float*  S_gamma,
                float*  efficiency);

#endif // AT_GAMMARESPONSE_H_
