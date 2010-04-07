/**
 * @file
 * @brief ...
 */

/*
 *    AT_KatzModel.c
 *    ===========================
 *
 *    Created on: 01.03.2010
 *    Author: kongruencja
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

#include "AT_KatzModel.h"


double AT_D_RDD_Gy_int( double  r_m,
    void*   params){
  double D_Gy;
  long  n_tmp              = 1;
  AT_P_RDD_parameters* par = (AT_P_RDD_parameters*)params;
  double fr_m               = r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  n_tmp,
      &fr_m,
      *(par->E_MeV_u),
      *(par->particle_no),
      *(par->material_no),
      *(par->rdd_model),
      par->rdd_parameters,
      *(par->er_model),
      par->er_parameters,
      &D_Gy);

  return (2.0 * M_PI * r_m * D_Gy);
}


double AT_sI_int( double  r_m,
    void*   params){
  double  P = AT_P_RDD(r_m, params);
  return (r_m * P);
}


double AT_P_RDD( double  r_m,
    void*   params)
{
  double  D_Gy;
  long   n_tmp              = 1;
  AT_P_RDD_parameters* par  = (AT_P_RDD_parameters*)params;
  double  fr_m               = r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  n_tmp,
      &fr_m,
      *(par->E_MeV_u),
      *(par->particle_no),
      *(par->material_no),
      *(par->rdd_model),
      par->rdd_parameters,
      *(par->er_model),
      par->er_parameters,
      &D_Gy);

  long gamma_model = GR_GeneralTarget;
  double P;
  AT_gamma_response(  n_tmp,
      &D_Gy,
      gamma_model,
      par->gamma_parameters,
      // return
      &P);
  return P;
}
