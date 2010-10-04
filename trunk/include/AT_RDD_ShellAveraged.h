#ifndef AT_RDD_SHELLAVERAGED_H_
#define AT_RDD_SHELLAVERAGED_H_

/**
 * @file
 * @brief
 * This file contains implementation of so called
 * Site RDD which is can be created from any PointLike
 * RDD by setting some fixed value of RDD in the core
 * in such way that whole RDD is LET-normalized.
 * TODO describe better...
 * Moreover some useful function as average
 * dose delivered in shell spanned between radius r1 and r2
 * are implemented here.
 * Here also inverse RDD (radius as function of dose) are implemented.
 */

/*
 *    AT_RDD_ShellAveraged.h
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
#include "AT_RDD_Simple.h"

/* --------------------------------------------------- SHELL AVERAGE DOSE ---------------------------------------------------*/

/**
 * Calculates average dose for "old" Katz RDD (derived from linear (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr
 *
 * thus:
 *
 * Dav(r1,r2) = coeff / (pi r2^2 - pi r1^2)  *  \int_r1^r2 rmax/r * (rmax/r - 1.) 2 pi r dr
 * Dav(r1,r2) = 2 * coeff * rmax^2/ (r2^2 - r1^2)  *  \int_r1^r2 (1/r - 1/rmax) dr
 * Dav(r1,r2) = 2 * coeff / ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 (1/r - 1/rmax) dr
 *
 * Let us calculate integral:
 *
 * \int_r1^r2 (1/r - 1/rmax) dr = log(r) - r/rmax |_r1^r2 = log(r2/r1) - (r2 - r1)/rmax
 *
 * in other words:
 *
 * Dav(r1,r2) = 2 * coeff * ( log(r2/r1) - (r2 - r1)/rmax ) / ((r2/rmax)^2 - (r1/rmax)^2)
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
inline double   AT_RDD_Katz_LinearER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy);


/**
 * Calculates average dose kernel for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * kernel = 2/(x2^2 - x1^2) * \int_x1^x2 1/x^2 * 1/alpha * (1 - x)^(1/alpha) x dx =
 *        = 2/(x2^2 - x1^2) * \int_x1^x2 1/x * 1/alpha * (1 - x)^(1/alpha) dx
 *
 * now we use the information that:
 *
 * \int 1/x * 1/alpha * (1 - x)^(1/alpha) dx = (1-x)^(1/alpha) ((x-1)/x)^(-1/alpha) _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x)+constant
 *
 * thus:
 *
 * kernel =  2/(x2^2 - x1^2) * (F2 - F1)
 *
 * where:
 *
 * F1 = (1-x1)^(1/alpha) ((x1-1)/x1)^(-1/alpha) _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x1)
 * F2 = (1-x2)^(1/alpha) ((x2-1)/x2)^(-1/alpha) _2F_1(-1/alpha,-1/alpha;(alpha-1)/alpha;1/x2)
 *
 * here _2F_1 is the special hypergeometric function
 *
 * @param[in] x1                     inner radius x1 (lower integration limit)
 * @param[in] x2                     outer radius x2 (upper integration limit)
 * @param[in] alpha                  parameter of ER model
 * @return kernel                    calculated kernel
 */
double   AT_RDD_Katz_PowerLawER_DaverageKernel(  const double x1,
    const double x2,
    const double alpha);


/**
 * Calculates approximate value of average dose kernel for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * kernel = 2/(x2^2 - x1^2) * \int_x1^x2 1/x^2 * 1/alpha * (1 - x)^(1/alpha) x dx =
 *        = 2/(x2^2 - x1^2) * \int_x1^x2 1/x * 1/alpha * (1 - x)^(1/alpha) dx
 *
 * now we use the series expansion:
 *
 * 1/x * 1/alpha * (1 - x)^(1/alpha) = 1 / (x alpha) - 1 / alpha^2  + (1/alpha  - 1) x / (2 alpha^2) + O(x^2)
 *
 * to calculate the integral:
 *
 * \int 1/x * 1/alpha * (1 - x)^(1/alpha) dx \approx x / alpha^2 ( (x / 4alpha) * (1/alpha - 1) - 1 ) + log(x) + C
 *
 * thus:
 *
 * kernel =  2/(x2^2 - x1^2) * (F2 - F1)
 *
 * where:
 *
 * F1 = x1 / alpha^2 ( (x1 / 4alpha) * (1/alpha - 1) - 1 ) + log(x1)
 * F2 = x2 / alpha^2 ( (x2 / 4alpha) * (1/alpha - 1) - 1 ) + log(x2)
 *
 * @param[in] x1                     inner radius x1 (lower integration limit)
 * @param[in] x2                     outer radius x2 (upper integration limit)
 * @param[in] alpha                  parameter of ER model
 * @return kernel                    calculated kernel
 */
double   AT_RDD_Katz_PowerLawER_DaverageKernel_approx(  const double x1,
    const double x2,
    const double alpha);


/**
 * Calculates average dose for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr
 *
 * thus:
 *
 * D(r) = coeff * kernel(r)
 *
 * Dav(r1,r2) = coeff * 2 / (r2^2 - r1^2) * \int_r1^r2 kernel(r) r dr
 *
 * substituting x1 = r1/rmax , x2 = r2/rmax we will have:
 *
 * Dav(r1,r2) = coeff * 2 / (x2^2 - x1^2) * \int_x1^x2 kernel(x) x dx
 *
 * in other words:
 *
 * Dav(r1,r2) = coeff * kernel_av( x1, x2 )
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] alpha                    parameter of ER model
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double   AT_RDD_Katz_PowerLawER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double alpha,
    const double Katz_point_coeff_Gy);


/**
 * Maximum electron range and beta
 */
typedef struct {
  double max_electron_range_m;
  double beta;
} AT_RDD_Cucinotta_Ddelta_average_integrand_m_parameters;


/**
 * Integrand in Cucinotta Ddelta calculation:
 *
 *  fS(r) * fL(r) /r;
 * @param[in] r_m                      distance [m]
 * @param[in] params                   vector of parameters (table of doubles)
 * @return integrand  fS(r) * fL(r) /r
 */
double AT_RDD_Cucinotta_Ddelta_average_integrand_m(  double r_m,
    void * params);


/**
 * Calculates average dose for Cucinotta delta RDD.
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 Ddelta(r) 2 pi r dr
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * Ddelta(r) = coeff * fS(r) * fL(r) * rmax^2/r^2
 *
 * Thus
 *
 * Dav(r1,r2) = 2 * coeff/ ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 fS(r) * fL(r) * 1/r dr
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double  AT_RDD_Cucinotta_Ddelta_average_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double beta,
    const double Katz_point_coeff_Gy);


/**
 * Calculates average dose for Cucinotta excitation RDD with Cnorm = 1
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 Dexc(r) 2 pi r dr
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * Dexc(r) = Cnorm * coeff * exp( - r / 2d ) * (rmax/r)^2            [where d = (beta/2) * (hbar * c / wr) and wr = 13eV ]
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 * Cnorm      =  1
 *
 * Thus
 *
 * Dav(r1,r2) = 2 coeff/ ((r2/rmax)^2 - (r1/rmax)^2) * \int_r1^r2 exp( - r / 2d ) * 1/r dr
 *
 * where:
 *
 * coeff      =  (C / 2 pi) * (Zeff/beta)^2 * 1/rho * 1 /rmax^2
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double   AT_RDD_Cucinotta_Dexc_average_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double beta,
    const double Katz_point_coeff_Gy);


/**
 * Calculates normalization constant
 * for Cucinotta RDD
 *
 * We should have:
 *
 * LET = 2 pi rho \int_rmin^rmax D(r) r dr
 *
 * Thus:
 *
 * LET = 2 pi rho \int_rmin^rmax Ddelta(r) r dr + 2 pi rho \int_rmin^rmax Dexc(r) r dr
 *
 * and
 *
 * LET / (2 pi rho) = (pi rmax^2 - pi rmin^2) Ddelta_average(rmin,rmax) + (pi rmax^2 - pi rmin^2) * Cnorm * Dexc_average(rmin,rmax)
 *
 * so
 *
 * LET / (2 pi rho * ((pi rmax^2 - pi rmin^2)) )  = Ddelta_average(rmin,rmax) + Cnorm * Dexc_average(rmin,rmax)
 *
 * finally:
 *
 * Cnorm = (LET / (2 pi rho * ((pi rmax^2 - pi rmin^2)) ) - Ddelta_average(rmin,rmax) ) / Dexc_average(rmin,rmax)
 *
 * @param[in] r_min_m                  minimum radius cut-off distance [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return C norm
 */
double   AT_RDD_Cucinotta_Cnorm( const double r_min_m,
    const double max_electron_range_m,
    const double beta,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double Katz_point_coeff_Gy);


/**
 * TODO
 * @param r1_m
 * @param r2_m
 * @param a0_m
 * @param max_electron_range_m
 * @param norm_Gy
 * @return
 */
double   AT_RDD_Geiss_average_Gy(  const double r1_m,
    const double r2_m,
    const double a0_m,
    const double max_electron_range_m,
    const double norm_Gy);

/* --------------------------------------------------- dEdx IN OUTER SHELL ---------------------------------------------------*/


/**
 * Calculates energy delivered to shell between radius a_0 and r_max
 * for "old" Katz RDD (derived from linear (on wmax) ER model).
 *
 * dEdx = rho \int_a0^rmax  D(r) 2 pi r dr =
 *      = rho * (pi rmax^2 - pi a0^2) D_av(a0,rmax)
 *
 * because:
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double   AT_RDD_Katz_LinearER_dEdx_J_m(  const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double Katz_point_coeff_Gy);


/**
 * Calculates energy delivered to shell between radius a_0 and r_max
 * for "new" Katz RDD (derived from power-law (on wmax) ER model).
 *
 * dEdx = rho \int_a0^rmax D(r) 2 pi r dr =
 *      = rho * (pi rmax^2 - pi a0^2) D_av(a0,rmax)
 *
 * because:
 *
 * Dav(r1,r2) = 1/ (pi r2^2 - pi r1^2) * \int_r1^r2 D(r) 2 pi r dr
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] alpha                    parameter of ER model
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double   AT_RDD_Katz_PowerLawER_dEdx_J_m(  const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double alpha,
    const double Katz_point_coeff_Gy);

/* --------------------------------------------------- SITE RDD ---------------------------------------------------*/

// TODO check why in SPIFF KatzSite and Edmund are not reproducing dose

/**
 * Calculates Site RDD, which is LET-normalized
 * for "old" Katz RDD (derived from linear (on wmax) ER model).
 *
 * Dsite(r) = 1 / (rho pi a0^2) * (LET - dEdx)   for r < a0
 * Dsite(r) = D(r)                               for r >= a0
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] dEdx_J_m                 dEdx [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double   AT_RDD_Katz_LinearER_DSite_Gy( const double r_m,
    const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy);


/**
 * Calculates Site RDD, which is LET-normalized
 * for "new" Katz RDD (derived from power-law (on wmax) ER model).
 *
 * Dsite(r) = 1 / (rho pi a0^2) * (LET - dEdx)   for r < a0
 * Dsite(r) = D(r)                               for r >= a0
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density rho [kg/m^3]
 * @param[in] alpha                    parameter of ER model
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] dEdx_J_m                 dEdx [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double   AT_RDD_Katz_PowerLawER_DSite_Gy( const double r_m,
    const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double alpha,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy);

/**
 * TODO
 * @param r_m
 * @param r_min_m                  minimum radius cut-off distance [m]
 * @param max_electron_range_m     delta electron maximum range rmax [m]
 * @param a0_m                     inner radius a0 (lower integration limit) [m]
 * @param er_model                 delta electron range model code number
 * @param alpha                    parameter of ER model
 * @param density_kg_m3            material density rho [kg/m^3]
 * @param LET_J_m                  LET [J/m]
 * @param dEdx_J_m                 dEdx [J/m]
 * @param Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return
 */
inline double  AT_RDD_KatzSite_Gy( const double r_m,
    const double r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double density_kg_m3,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy);


/**
 * TODO
 * @param D_Gy
 * @param r_min_m                  minimum radius cut-off distance [m]
 * @param max_electron_range_m     delta electron maximum range rmax [m]
 * @param a0_m                     radius target [m]
 * @param er_model                 delta electron range model code number
 * @param alpha                    parameter of ER model
 * @param Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return
 */
double  AT_inverse_RDD_KatzSite_m( const double D_Gy,
    const double r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double d_max_Gy,
    const double Katz_point_coeff_Gy);


#endif /* AT_RDD_SHELLAVERAGED_H_ */
