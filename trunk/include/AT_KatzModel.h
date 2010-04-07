#ifndef AT_KATZMODEL_H_
#define AT_KATZMODEL_H_

/**
 * @file
 * @brief Katz model algorithm
 */

/*
 *    AT_KatzModel.h
 *    ===========================
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

#include "AT_RDD.h"
#include "AT_GammaResponse.h"


/**
 * TODO
 */
double         AT_P_RDD(                    double  r_m,
    void* params);


/**
 * TODO
 */
double         AT_sI_int(                   double  r_m,
    void* params);


/**
 * TODO
 */
double         AT_D_RDD_Gy_int(             double  r_m,
    void* params);


/**
 * TODO
 */
typedef struct {
  /* radiation field parameters */
  double*  E_MeV_u;
  long*    particle_no;
  /* detector parameters */
  long*    material_no;
  /* radial dose distribution model */
  long*    rdd_model;
  double*  rdd_parameters;
  /* electron range model */
  long*    er_model;
  double*  er_parameters;
  /* gamma response parameter*/
  double   gamma_parameters[5];
} AT_P_RDD_parameters;


#endif /* AT_KATZMODEL_H_ */
