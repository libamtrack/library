#ifndef AT_RDD_SIMPLE_H_
#define AT_RDD_SIMPLE_H_

/**
 * @brief Radial Dose Distribution models, mostly point-like.
 * Radial Dose Distribution describes dose deposition around
 * particle track: dose D in point as a function of distance r.
 * TODO put better description
 * This file contains single-particle case implementation,
 * for fast multi-particle case see AT_RDD.h file.
 * Here also inverse RDD (radius as function of dose) are implemented.
 */

/*
 *    AT_RDD_Simple.h
 *    ========
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

#include <string.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h>

#include "AT_Constants.h"
#include "AT_DataLET.h"
#include "AT_PhysicsRoutines.h"

/* --------------------------------------------------- TEST RDD ---------------------------------------------------*/

/**
 * Test RDD, flat for given range of r.
 * @param r_m
 * @param r_min_m                  minimum radius cut-off distance [m]
 * @param r_max_m
 * @param norm_constant_Gy
 * @return
 */
 double AT_RDD_Test_Gy( const double r_m,
    const double r_min_m,
    const double r_max_m,
    const double norm_constant_Gy);

/**
 * TODO
 * @param D_Gy
 * @param r_max_m                  delta electron maximum range rmax [m]
 * @return
 */
 double AT_inverse_RDD_Test_m( const double D_Gy,
    const double r_max_m);

/* --------------------------------------------------- GEISS RDD ---------------------------------------------------*/

/**
 * TODO
 * @param r_m
 * @param r_min_m                  minimum radius cut-off distance [m]
 * @param max_electron_range_m
 * @param a0_m
 * @param norm_constant_Gy
 * @return
 */
 double AT_RDD_Geiss_Gy( const double r_m,
    const double r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const double norm_constant_Gy);


/**
 * TODO
 * @param D_Gy
 * @param d_min_Gy
 * @param d_max_Gy
 * @param a0_m
 * @param norm_constant_Gy
 * @return
 */
 double AT_inverse_RDD_Geiss_m( const double D_Gy,
    const double d_min_Gy,
    const double d_max_Gy,
    const double a0_m,
    const double norm_constant_Gy);


/* --------------------------------------------------- KATZ RDD ---------------------------------------------------*/


/**
 * Katz C constant given by equation:
 *
 * C = 2 pi N e^4 / ( m c^2 (4 pi eps_0)^2)
 *
 * In other form:
 *
 * C = N * C1 , where C1 = e^4 / ( 8 m pi (c * eps_0)^2)
 *
 * For water: C = 1.36662e-12 [J/m] = 8.53 [MeV/m]
 */
static const double AT_Katz_C1_J_m2 = \
    GSL_CONST_MKSA_ELECTRON_CHARGE * GSL_CONST_MKSA_ELECTRON_CHARGE * GSL_CONST_MKSA_ELECTRON_CHARGE * GSL_CONST_MKSA_ELECTRON_CHARGE * \
    M_1_PI * 0.125 / \
    (GSL_CONST_MKSA_MASS_ELECTRON * GSL_CONST_MKSA_SPEED_OF_LIGHT * GSL_CONST_MKSA_SPEED_OF_LIGHT * \
        GSL_CONST_MKSA_VACUUM_PERMITTIVITY * GSL_CONST_MKSA_VACUUM_PERMITTIVITY);


/**
 * Calculates coefficient
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] Z_eff                    effective ion charge Zeff
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] density_kg_m3   material density rho [kg/m^3]
 * @param[in] electron_density_m3      electron density of given material [1/m^3]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @return coeff [Gy]                  calculated coefficient
 */
 double AT_RDD_Katz_coeff_Gy( const double Z_eff,
    const double beta,
    const double density_kg_m3,
    const double electron_density_m3,
    const double max_electron_range_m);


/**
 * TODO
 * @param E_MeV_u
 * @param particle_no
 * @param material_no
 * @param er_model
 * @return
 */
 double AT_RDD_Katz_coeff_Gy_general(     const double  E_MeV_u,
    const long    particle_no,
    const long    material_no,
    const long    er_model);


/**
 * Calculates "old" Katz RDD (derived from linear (on wmax) ER model):
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho *  1/r * (1/r - 1/rmax)
 *
 * or:
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/rmax^2 * rmax/r * (rmax/r - 1.)
 *
 * thus:
 *
 * D(r) = coeff * rmax/r * (rmax/r - 1.)
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r_m                      distance r [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] radial dose distribution at distance r
 */
 double AT_RDD_Katz_LinearER_Dpoint_Gy(        const double r_m,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy);


/**
 * Calculates "new" Katz RDD (derived from power-law (on wmax) ER model):
 *
 * D(r) = (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1/r^2 * 1/alpha * (1 - r/rmax)^(1/alpha)
 *
 * thus:
 *
 * D(r) = coeff * 1/alpha * (rmax/r)^2 * (1 - r/rmax)^(1/alpha)
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r_m                      distance r [m]
 * @param[in] alpha                    parameter of ER model
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] radial dose distribution at distance r
 */
 double AT_RDD_Katz_PowerLawER_Dpoint_Gy(        const double r_m,
    const double alpha,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy);


/**
 * Single particle implementation of Katz model RDD
 * @param r_m                      distance r [m]
 * @param r_min_m                  minimum radius cut-off distance [m]
 * @param max_electron_range_m     delta electron maximum range rmax [m]
 * @param er_model                 delta electron range model code number
 * @param alpha                    parameter of ER model
 * @param Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return
 */
 double AT_RDD_KatzPoint_Gy( const double r_m,
    const double r_min_m,
    const double max_electron_range_m,
    const long   er_model,
    const double alpha,
    const double Katz_point_coeff_Gy);


/**
 * TODO
 */
typedef struct {
 double  D_Gy;
 long    er_model;
 double  r_min_m;
 double  max_electron_range_m;
 double  alpha;
 double  Katz_point_coeff_Gy;
} AT_inverse_RDD_KatzPoint_PowerLawER_parameters;


/**
 * TODO
 * @param D_Gy
 * @param r_min_m                  minimum radius cut-off distance [m]
 * @param max_electron_range_m     delta electron maximum range rmax [m]
 * @param Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return
 */
double AT_inverse_RDD_KatzPoint_LinearER_m( const double D_Gy,
   const double r_min_m,
   const double max_electron_range_m,
   const double Katz_point_coeff_Gy);


/**
 * TODO
 * @param r_m
 * @param params
 * @return
 */
double AT_inverse_RDD_KatzPoint_PowerLawER_solver_function_Gy( const double r_m , void * params );


/**
 * Single particle implementation of inverse Katz model RDD
 * @param D_Gy
 * @param r_min_m                  minimum radius cut-off distance [m]
 * @param max_electron_range_m     delta electron maximum range rmax [m]
 * @param er_model                 delta electron range model code number
 * @param alpha                    parameter of ER model
 * @param Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return
 */
double AT_inverse_RDD_KatzPoint_m( const double D_Gy,
   const double r_min_m,
   const double max_electron_range_m,
   const long   er_model,
   const double alpha,
   const double Katz_point_coeff_Gy);


/* --------------------------------------------------- CUCINOTTA RDD ---------------------------------------------------*/

/**
 * Calculates short range modification function fS(r)
 * for Cucinotta RDD
 *
 * fS(r) = 1.0/( r0/r + 0.6 + 1.7 beta + 1.1 beta^2)           [here r0 = 10^(-9) [m]]
 *
 * @param[in] r_m                      distance [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @return fS(r)
 */
 double AT_RDD_Cucinotta_f_shortRange( const double r_m,
    const double beta);


/**
 * Calculates long range modification function fL(r)
 * for Cucinotta RDD
 *
 * fL(r) = exp( -(r/(0.37rmax))^2 )
 *
 * @param[in] r_m                      distance [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @return fL(r)
 */
 double AT_RDD_Cucinotta_f_longRange( const double r_m,
    const double max_electron_range_m);


/**
 * Calculates radial component D_delta
 * for Cucinotta RDD
 *
 * Ddelta(r) = C z^2 / beta^2 1/rho fS(r) fL(r) /r^2
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * Ddelta(r) = coeff * fS(r) * fL(r) * rmax^2/r^2
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r_m                      distance [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return Ddelta(r) [Gy]
 */
 double AT_RDD_Cucinotta_Ddelta_Gy( const double r_m,
    const double max_electron_range_m,
    const double beta,
    const double Katz_point_coeff_Gy);


/**
 * In Cucinotta RDD we have: d = (beta/2) * (hbar * c / wr) and wr = 13eV
 * Part of this can be saved as a constant: 0.5 * hbar * c / wr , where wr = 13eV
 */
static const double AT_RDD_Cucinotta_C_dm_wr  =  0.5 * GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR * GSL_CONST_MKSA_SPEED_OF_LIGHT / (13.0 * GSL_CONST_MKSA_ELECTRON_VOLT);


/**
 * Calculates excitation component D_exc
 * for Cucinotta RDD
 *
 * Dexc(r) = C exp( - r / 2d ) / r^2            [where d = (beta/2) * (hbar * c / wr) and wr = 13eV ]
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * Dexc(r) = Cnorm * coeff * exp( - r / 2d ) * (rmax/r)^2
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * Cnorm      =  Cnorm = (LET / (2 pi rho * ((pi rmax^2 - pi rmin^2)) ) - Ddelta_average(rmin,rmax) ) / Dexc_average(rmin,rmax)
 *
 * @param[in] r_m                      distance [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] C_norm                   TODO
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return Dexc(r) [Gy]
 */
double AT_RDD_Cucinotta_Dexc_Gy( const double r_m,
    const double max_electron_range_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy);


/**
 * Calculates Cucinotta point RDD
 *
 * D(r)    = Dexc(r) + Ddelta(r)
 *
 * Ddelta(r) = C z^2 / beta^2 1/rho fS(r) fL(r) /r^2
 *
 * Dexc(r) = C exp( - r / 2d ) / r^2
 *
 * @param[in] r_m                      distance [m]
 * @param[in] r_min_m                  minimum cut-off [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] C_norm                   TODO
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy]
 */
 double AT_RDD_CucinottaPoint_Gy( const double r_m,
    const double r_min_m,
    const double max_electron_range_m,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy);


/**
 * TODO
 */
typedef struct {
  double  D_Gy;
  long    er_model;
  double  r_min_m;
  double  max_electron_range_m;
  double  beta;
  double  C_norm;
  double  Katz_point_coeff_Gy;
} AT_inverse_RDD_Cucinotta_parameters;


/**
 * TODO
 * @param r_m
 * @param params
 * @return
 */
double AT_inverse_RDD_Cucinotta_solver_function_Gy( const double r_m , void * params );


/**
 * TODO
 * @param D_Gy
 * @param r_min_m                  minimum cut-off [m]
 * @param max_electron_range_m     delta electron maximum range rmax [m]
 * @param er_model
 * @param beta                     relative ion speed beta = v/c
 * @param C_norm
 * @param Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return
 */
double AT_inverse_RDD_Cucinotta_m( const double D_Gy,
    const double r_min_m,
    const double max_electron_range_m,
    const long   er_model,
    const double beta,
    const double C_norm,
    const double Katz_point_coeff_Gy);


#endif /* AT_RDD_SIMPLE_H_ */
