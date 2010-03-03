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
 * Geometry function
 *
 * @param[in]        r0_m     distance from the center of the ion to the center of the target [m]
 * @param[in]        a0_m     radius of the target [m]
 * @param[in]        r_m      distance at which segment length is calculated [m]
 * @return   Phi(r0,a0,r)
 */
double          geometryFunctionPhi(         const double r0_m,
    const double a0_m,
    const double r_m);


/**
 * Returns RDD as a function of distance r_m for target with radius a0_m
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
void AT_RDD_ExtendedTarget_Gy( const long*  n,
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


/**
 * TODO
 */
double   AT_RDD_Katz_ext_kernel_Gy(   const double t_m,
    const double r_m,
    const double a0_m,
    const double alpha,
    const double r_min_m,
    const double r_max_m,
    const double Katz_point_coeff_Gy);

/**
 * TODO
 */
double         AT_RDD_Katz_ext_integrand_Gy(double t_m,
    void * params);

/**
 * TODO
 */
double         AT_RDD_Katz_ext_Gy(          const double r_m,
    const double a0_m,
    const double alpha,
    const double r_min_m,
    const double r_max_m,
    const double Katz_point_coeff_Gy);


//TODO implement inverse extended target RDD


#endif /* AT_RDD_EXTENDEDTARGET_H_ */
