#ifndef AT_RDD_EXTENDEDTARGET_H_
#define AT_RDD_EXTENDEDTARGET_H_

/**
 * @brief
 * File contains implementation of RDD
 * averaged over finite size, circular target of radius a0.
 * Here also inverse RDD (radius as function of dose) are implemented.
 */

/*
 *    AT_ExtendedTarget.h
 *    ===========================
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

#include "AT_RDD_Simple.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


/**
 * Geometry function Phi \n
 *
 * Length of arc segment is given as  2 Phi(r,a0,t) r
 *
 * Phi(r,a0,t) = 2 arctan( sqrt( (a0^2 - (r-t)^2) / ((t + r)^2  - a0^2) ) )   for t > | r - a0 |
 *
 * Phi(r,a0,t) = pi                                                          for t <= | r - a0 |
 *
 * @param[in]        r_m      distance from the center of the ion to the center of the target [m]
 * @param[in]        a0_m     radius of the target [m]
 * @param[in]        t_m      distance at which segment length is calculated [m]
 * @return   Phi(r,a0,t)
 */
double geometryFunctionPhi(         const double r_m,
    const double a0_m,
    const double t_m);


/* --------------------------------------------------- KATZ EXT TARGET RDD ---------------------------------------------------*/


/**
 * TODO
 */
typedef struct {
  double  r_m;
  double  a0_m;
  double  max_electron_range_m;
  long    er_model;
  double  alpha;
  double  KatzPoint_r_min_m;
  double  Katz_point_coeff_Gy;
} AT_RDD_ExtendedTarget_KatzPoint_parameters;


/**
 * TODO
 * @param t_m
 * @param params
 * @return
 */
double AT_RDD_ExtendedTarget_KatzPoint_integrand_Gy(
    double t_m,
    void* params);


/**
 * TODO
 * @param r_m
 * @param a0_m
 * @param er_model
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param alpha
 * @param Katz_point_coeff_Gy
 * @return
 */
double AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration(
    const double  r_m,
    const double  a0_m,
    const long    er_model,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  alpha,
    const double  Katz_point_coeff_Gy);


/**
 * TODO
 * @param r_m
 * @param a0_m
 * @param er_model
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param alpha
 * @param Katz_plateau_Gy
 * @param Katz_point_coeff_Gy
 * @return
 */
double AT_RDD_ExtendedTarget_KatzPoint_Gy(
    const double  r_m,
    const double  a0_m,
    const long    er_model,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  alpha,
    const double  Katz_plateau_Gy,
    const double  Katz_point_coeff_Gy);


/**
 * TODO
 */
typedef struct {
  double  D_Gy;
  double  a0_m;
  double  max_electron_range_m;
  long    er_model;
  double  alpha;
  double  Katz_plateau_Gy;
  double  Katz_point_r_min_m;
  double  Katz_point_coeff_Gy;
} AT_inverse_RDD_ExtendedTarget_KatzPoint_parameters;


/**
 * TODO
 * @param r_m
 * @param params
 * @return
 */
double AT_inverse_RDD_ExtendedTarget_KatzPoint_solver_function_Gy( const double r_m , void * params );


/**
 * TODO
 * @param D_Gy
 * @param r_min_m
 * @param max_electron_range_m
 * @param a0_m
 * @param er_model
 * @param alpha
 * @param Katz_plateau_Gy
 * @param Katz_point_coeff_Gy
 * @return
 */
double AT_inverse_RDD_ExtendedTarget_KatzPoint_m( const double D_Gy,
    const double r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double Katz_plateau_Gy,
    const double Katz_point_coeff_Gy);

/* --------------------------------------------------- CUCINOTTA EXT TARGET RDD ---------------------------------------------------*/

/**
 * TODO
 */
typedef struct {
  double  r_m;
  double  a0_m;
  double  KatzPoint_r_min_m;
  double  max_electron_range_m;
  double  beta;
  double  Katz_point_coeff_Gy;
  double  C_norm;
} AT_RDD_ExtendedTarget_CucinottPoint_parameters;


/**
 * TODO
 * @param t_m
 * @param params
 * @return
 */
double AT_RDD_ExtendedTarget_CucinottaPoint_integrand_Gy(
    double t_m,
    void* params);


/**
 * TODO
 * @param r_m
 * @param a0_m
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param beta
 * @param Katz_point_coeff_Gy
 * @param C_norm
 * @return
 */
double AT_RDD_ExtendedTarget_CucinottaPoint_Gy_by_integration(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  Katz_point_coeff_Gy,
    const double  C_norm);


/**
 * TODO
 * @param r_m
 * @param a0_m
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param beta
 * @param Katz_point_coeff_Gy
 * @param C_norm
 * @param Cucinotta_plateau_Gy
 * @return
 */
double AT_RDD_ExtendedTarget_CucinottaPoint_Gy(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  Katz_point_coeff_Gy,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy);


/**
 * TODO
 */
typedef struct {
  double  D_Gy;
  double  a0_m;
  double  KatzPoint_r_min_m;
  double  max_electron_range_m;
  double  beta;
  double  Katz_point_coeff_Gy;
  double  C_norm;
  double  Cucinotta_plateau_Gy;
} AT_inverse_RDD_ExtendedTarget_CucinottaPoint_parameters;


/**
 * TODO
 * @param r_m
 * @param params
 * @return
 */
double AT_inverse_RDD_ExtendedTarget_CucinottaPoint_solver_function_Gy( const double r_m , void * params );


/**
 * TODO
 * @param D_Gy
 * @param a0_m
 * @param KatzPoint_r_min_m
 * @param max_electron_range_m
 * @param beta
 * @param Katz_point_coeff_Gy
 * @param C_norm
 * @param Cucinotta_plateau_Gy
 * @return
 */
double AT_inverse_RDD_ExtendedTarget_CucinottaPoint_m( const double D_Gy,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  Katz_point_coeff_Gy,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy);


#endif /* AT_RDD_EXTENDEDTARGET_H_ */
