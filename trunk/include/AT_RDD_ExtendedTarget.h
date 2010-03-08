#ifndef AT_RDD_EXTENDEDTARGET_H_
#define AT_RDD_EXTENDEDTARGET_H_

/**
 * @file
 * @brief Extended Target RDDs
 */

/*
 *    AT_ExtendedTarget.h
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

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

///////////////////////////////////////////// SITE /////////////////////////////////////////////

/**
 * Returns RDD as a function of distance r_m for target with radius a0_m
 * approximated version, called RDD_Site
 *
 * @param[in]   n
 * @param[in]   r_m            distance [m]
 * @param[in]   a0_m           target radius
 * @param[in]   E_MeV_u
 * @param[in]   particle_no
 * @param[in]   material_no
 * @param[in]   rdd_model
 * @param[in]   rdd_parameter
 * @param[in]   er_model
 * @param[in]   er_parameter
 * @param[out]  D_RDD_Gy       dose [Gy]
 */
void AT_RDD_Site_Gy( const long*  n,
    const float*  r_m,
    const float*  a0_m,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const long*   particle_no,
    /* detector parameters */
    const long*   material_no,
    /* radial dose distribution model */
    const long*   rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const long*   er_model,
    const float*  er_parameter,
    float*        D_RDD_Gy);

///////////////////////////////////////////// EXTENDED TARGET /////////////////////////////////////////////

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
double          geometryFunctionPhi(         const double r_m,
    const double a0_m,
    const double t_m);


/**
 * Returns RDD as a function of distance r_m for target with radius a0_m
 *
 * Extended target RDD is defined as:
 *
 * Dext( r , a0 ) = 1 / (pi a0^2) \int_low^up D(t) 2 Phi(t,a0,r) t dt
 *
 * where
 * low = r-a0   for r > a0
 * low = 0      for r <= a0
 *
 * up = r + a0
 *
 * For low values of r (r << a0) we have Dext(r,a0) ~= Dext(0,a0) =  1 / (pi a0^2) \int_0^a0 D(t) 2 pi t dt
 *
 * For large values of r (r >> a0) we have Dext(r,a0) ~= D(r)
 *
 * @param[in]   n
 * @param[in]   r_m            distance [m]
 * @param[in]   a0_m           target radius
 * @param[in]   E_MeV_u
 * @param[in]   particle_no
 * @param[in]   material_no
 * @param[in]   rdd_model
 * @param[in]   rdd_parameter
 * @param[in]   er_model
 * @param[in]   er_parameter
 * @param[out]  D_RDD_Gy       dose [Gy]
 */
void AT_RDD_ExtendedTarget_Gy( const long  n,
    const float* r_m,
    const float  a0_m,
    /* radiation field parameters */
    const float  E_MeV_u,
    const long   particle_no,
    /* detector parameters */
    const long   material_no,
    /* radial dose distribution model */
    const long   rdd_model,
    const float* rdd_parameter,
    /* electron range model */
    const long   er_model,
    const float* er_parameter,
    float*       D_RDD_Gy);


/**
 * Katz extended target integrand
 * D(t) 2 Phi(t,a0,r) t
 *
 * @param[in]   t_m
 * @param[in]   params  structure AT_RDD_ExtendedTarget_parameters
 */
double         AT_RDD_Katz_ext_integrand_Gy(double t_m,
    void * params);

/**
 * Calculates directly
 *
 * Dext( r , a0 ) = 1 / (pi a0^2) \int_low^up D(t) 2 Phi(t,a0,r) t dt
 *
 */
double AT_RDD_ExtendedTarget_integrate_Gy(  const double r_m,
    const double a0_m,
    const double r_min_m,
    const double r_max_m,
    const float  E_MeV_u,
    const long   particle_no,
    /* detector parameters */
    const long   material_no,
    /* radial dose distribution model */
    const long   rdd_model,
    const float* rdd_parameter,
    /* electron range model */
    const long   er_model,
    const float* er_parameter);


//TODO implement inverse extended target RDD

typedef struct {
  float   r_m;
  float   a0_m;
  /* radiation field parameters */
  float   E_MeV_u;          /**< energy per nucleon */
  long    particle_no;
  /* detector parameters */
  long    material_no;
  /* radial dose distribution model */
  long    rdd_model;
  float*  rdd_parameter;
  /* electron range model */
  long    er_model;
  float*   er_parameter;
} AT_RDD_ExtendedTarget_parameters;

#endif /* AT_RDD_EXTENDEDTARGET_H_ */
