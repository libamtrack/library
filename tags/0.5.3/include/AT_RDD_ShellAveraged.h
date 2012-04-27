#ifndef AT_RDD_SHELLAVERAGED_H_
#define AT_RDD_SHELLAVERAGED_H_

/**
 * @brief RDD Shell averaged
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
 * @f[ Dav(r1,r2) = 1/(\pi r_2^2 - \pi r_1^2) * \int_{r_1}^{r_2} D(r) 2 \pi r dr @f]
 *
 * thus:
 *
 * @f[ Dav(r1,r2) = coeff / (\pi r_2^2 - \pi r_1^2)  *  \int_{r_1}^{r_2} rmax/r * (rmax/r - 1.) 2 \pi r dr @f]
 * @f[ Dav(r1,r2) = 2 * coeff * rmax^2/ (r_2^2 - r_1^2)  *  \int_{r_1}^{r_2} (1/r - 1/rmax) dr @f]
 * @f[ Dav(r1,r2) = 2 * coeff / ((r_2/rmax)^2 - (r_1/rmax)^2) * \int_{r_1}^{r_2} (1/r - 1/rmax) dr @f]
 *
 * Let us calculate integral:
 *
 * @f[ \int_{r_1}^{r_2} (1/r - 1/rmax) dr = log(r) - r/rmax |_{r_1}^{r_2} = log(r_2/r_1) - (r_2 - r_1)/rmax @f]
 *
 * in other words:
 *
 * @f[ Dav(r_1,r_2) = 2 * coeff * ( log(r_2/r_1) - (r_2 - r_1)/rmax ) / ((r_2/rmax)^2 - (r_1/rmax)^2) @f]
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
 double AT_RDD_Katz_LinearER_Daverage_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double Katz_point_coeff_Gy);


/**
 * Calculates average dose kernel for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius @f$r_1@f$ and @f$r_2@f$
 *
 * @f[kernel = 2/(x_2^2 - x_1^2) * \int_{x_1}^{x_2} 1/x^2 * 1/\alpha * (1 - x)^{1/\alpha} x dx = @f]
 * @f[       = 2/(x_2^2 - x_1^2) * \int_{x_1}^{x_2} 1/x * 1/\alpha * (1 - x)^{1/\alpha} dx @f]
 *
 * now we use the information that:
 *
 * @f[ \int 1/x * 1/\alpha * (1 - x)^{1/\alpha} dx = (1-x)^{1/\alpha} ((x-1)/x)^{-1/\alpha} _2F_1(-1/\alpha,-1/\alpha;(\alpha-1)/\alpha;1/x)+constant @f]
 *
 * thus:
 *
 * @f[ kernel =  2/(x_2^2 - x_1^2) * (F2 - F1) @f]
 *
 * where:
 *
 * @f[ F1 = (1-x_1)^{1/\alpha} ((x_1-1)/x_1)^{-1/\alpha} _2F_1(-1/\alpha,-1/\alpha;(\alpha-1)/\alpha;1/x_1) @f]
 * @f[ F2 = (1-x_2)^{1/\alpha} ((x_2-1)/x_2)^{-1/\alpha} _2F_1(-1/\alpha,-1/\alpha;(\alpha-1)/\alpha;1/x_2) @f]
 *
 * here @f$_2F_1@f$ is the special hypergeometric function
 *
 * @param[in] x1                     inner radius x1 (lower integration limit)
 * @param[in] x2                     outer radius x2 (upper integration limit)
 * @param[in] alpha                  parameter of ER model
 * @return                           calculated kernel
 */
double AT_RDD_Katz_PowerLawER_DaverageKernel(  const double x1,
    const double x2,
    const double alpha);


/**
 * Calculates approximate value of average dose kernel for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * @f[ kernel = 2/(x_2^2 - x_1^2) * \int_{x_1}^{x_2} 1/x^2 * 1/\alpha * (1 - x)^(1/\alpha) x dx = @f]
 * @f[        = 2/(x_2^2 - x_1^2) * \int_{x_1}^{x_2} 1/x * 1/\alpha * (1 - x)^(1/\alpha) dx @f]
 *
 * now we use the series expansion:
 *
 * @f[ 1/x * 1/\alpha * (1 - x)^(1/\alpha) = 1 / (x \alpha) - 1 / \alpha^2  + (1/\alpha  - 1) x / (2 \alpha^2) + O(x^2) @f]
 *
 * to calculate the integral:
 *
 * @f[ \int 1/x * 1/\alpha * (1 - x)^(1/\alpha) dx \approx x / \alpha^2 ( (x / 4\alpha) * (1/\alpha - 1) - 1 ) + log(x) + C @f]
 *
 * thus:
 *
 * @f[ kernel =  2/(x_2^2 - x_1^2) * (F2 - F1) @f]
 *
 * where:
 *
 * @f[ F1 = x_1 / \alpha^2 ( (x_1 / 4\alpha) * (1/\alpha - 1) - 1 ) + log(x_1) @f]
 * @f[ F2 = x_2 / \alpha^2 ( (x_2 / 4\alpha) * (1/\alpha - 1) - 1 ) + log(x_2) @f]
 *
 * @param[in] x1                     inner radius x1 (lower integration limit)
 * @param[in] x2                     outer radius x2 (upper integration limit)
 * @param[in] alpha                  parameter of ER model
 * @return kernel                    calculated kernel
 */
double AT_RDD_Katz_PowerLawER_DaverageKernel_approx(  const double x1,
    const double x2,
    const double alpha);


/**
 * Calculates average dose for "new" Katz RDD (derived from power-law (on wmax) ER model).
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * @f[ Dav(r1,r2) = 1/ (\pi r2^2 - \pi r1^2) * \int_{r_1}^{r_2} D(r) 2 \pi r dr @f]
 *
 * thus:
 *
 * @f[ D(r) = coeff * kernel(r) @f]
 *
 * @f[ Dav(r1,r2) = coeff * 2 / (r2^2 - r1^2) * \int_{r_1}^{r_2} kernel(r) r dr @f]
 *
 * substituting x1 = r1/rmax , x2 = r2/rmax we will have:
 *
 * @f[ Dav(r1,r2) = coeff * 2 / (x2^2 - x1^2) * \int_x1^x2 kernel(x) x dx @f]
 *
 * in other words:
 *
 * @f[ Dav(r1,r2) = coeff * kernel_av( x1, x2 ) @f]
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] alpha                    parameter of ER model
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double AT_RDD_Katz_PowerLawER_Daverage_Gy(  const double r1_m,
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
 * @param[in] params                   vector of parameters
 * @return integrand  fS(r) * fL(r) /r
 */
double AT_RDD_Cucinotta_Ddelta_average_integrand_m(  double r_m,
    void* params);


/**
 * Calculates average dose for Cucinotta delta RDD.
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * @f[ Dav(r1,r2) = 1/ (\pi r_2^2 - \pi r_1^2) * \int_{r_1}^{r_2} Ddelta(r) 2 \pi r dr @f]
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * @f[ Ddelta(r) = coeff * fS(r) * fL(r) * rmax^2/r^2 @f]
 *
 * Thus
 *
 * @f[ Dav(r_1,r_2) = 2 * coeff/ ((r_2/rmax)^2 - (r_1/rmax)^2) * \int_{r_1}^{r_2} fS(r) * fL(r) * 1/r dr @f]
 *
 * where:
 *
 * @f[ coeff      =  (C / 2 \pi) * (Zeff/beta)^2 * 1/\rho * 1 /rmax^2 @f]
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double AT_RDD_Cucinotta_Ddelta_average_Gy(  const double r1_m,
    const double r2_m,
    const double max_electron_range_m,
    const double beta,
    const double Katz_point_coeff_Gy);


/**
 * Calculates average dose for Cucinotta excitation RDD with Cnorm = 1
 * Here averaging is done over a shell between radius r_1 and r_2
 *
 * @f[ Dav(r_1,r_2) = 1/ (\pi r_2^2 - \pi r_1^2) * \int_{r_1}^{r_2} Dexc(r) 2 \pi r dr @f]
 *
 * We calculate using pre-calculated constant in following manner:
 *
 * @f[ Dexc(r) = Cnorm * coeff * exp( - r / 2d ) * (rmax/r)^2 @f]
 *
 * where:
 *
 * @f[ d = (beta/2) * (hbar * c / wr) @f]
 * @f[ wr = 13eV @f]
 * @f[ coeff      =  (C / 2 \pi) * (Zeff/beta)^2 * 1/\rho * 1 /rmax^2 @f]
 * @f[ Cnorm      =  1 @f]
 *
 * Thus
 *
 * @f[ Dav(r_1,r_2) = 2 coeff/ ((r_2/rmax)^2 - (r_1/rmax)^2) * \int_{r_1}^{r_2} exp( - r / 2d ) * 1/r dr @f]
 *
 * where:
 *
 * @f[ coeff      =  (C / 2 \pi) * (Zeff/beta)^2 * 1/\rho * 1 /rmax^2 @f]
 *
 * @param[in] r1_m                     inner radius r1 (lower integration limit) [m]
 * @param[in] r2_m                     outer radius r2 (upper integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed beta = v/c
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return D(r) [Gy] average radial dose distribution between r1 and r2
 */
double AT_RDD_Cucinotta_Dexc_average_Gy(  const double r1_m,
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
 * @f[ LET = 2 \pi \rho \int_{rmin}^{rmax} D(r) r dr @f]
 *
 * Thus:
 *
 * @f[ LET = 2 \pi \rho \int_{rmin}^{rmax} Ddelta(r) r dr + 2 \pi \rho \int_{rmin}^{rmax} Dexc(r) r dr @f]
 *
 * and
 *
 * @f[ LET / (2 \pi \rho) = (\pi rmax^2 - \pi rmin^2) Ddelta_{average}(rmin,rmax) + (\pi rmax^2 - \pi rmin^2) * Cnorm * Dexc_{average}(rmin,rmax) @f]
 *
 * so
 *
 * @f[ LET / (2 \pi \rho * ((\pi rmax^2 - \pi rmin^2)) )  = Ddelta_{average}(rmin,rmax) + Cnorm * Dexc_{average}(rmin,rmax) @f]
 *
 * finally:
 *
 * @f[ Cnorm = (LET / (2 \pi \rho * ((\pi rmax^2 - \pi rmin^2)) ) - Ddelta_{average}(rmin,rmax) ) / Dexc_{average}(rmin,rmax) @f]
 *
 * @param[in] r_min_m                  minimum radius cut-off distance [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] beta                     relative ion speed @f$\beta@f$ = v/c
 * @param[in] material_density_kg_m3   material density @f$\rho@f$ [kg/m^3]
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return C norm
 */
double AT_RDD_Cucinotta_Cnorm( const double r_min_m,
    const double max_electron_range_m,
    const double beta,
    const double material_density_kg_m3,
    const double LET_J_m,
    const double Katz_point_coeff_Gy);


/**
 * TODO
 * @param[in] r1_m
 * @param[in] r2_m
 * @param[in] a0_m
 * @param[in] max_electron_range_m
 * @param[in] norm_Gy
 * @return
 */
double AT_RDD_Geiss_average_Gy(  const double r1_m,
    const double r2_m,
    const double a0_m,
    const double max_electron_range_m,
    const double norm_Gy);

/* --------------------------------------------------- dEdx IN OUTER SHELL ---------------------------------------------------*/


/**
 * Calculates energy delivered to shell between radius a_0 and r_max
 * for "old" Katz RDD (derived from linear (on wmax) ER model).
 *
 * @f[ dEdx = \rho \int_{a_0}^{rmax}  D(r) 2 \pi r dr = @f]
 * @f[      = \rho * (\pi rmax^2 - \pi a_0^2) D_{av}(a_0,rmax) @f]
 *
 * because:
 *
 * @f[ Dav(r1,r2) = 1/ (\pi r_2^2 - \pi r_1^2) * \int_{r_1}^{r_2} D(r) 2 \pi r dr @f]
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density @f$\rho@f$ [kg/m^3]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double AT_RDD_Katz_LinearER_dEdx_J_m(  const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double Katz_point_coeff_Gy);


/**
 * Calculates energy delivered to shell between radius a_0 and r_max
 * for "new" Katz RDD (derived from power-law (on wmax) ER model).
 *
 * @f[ dEdx = \rho \int_{a_0}^{rmax} D(r) 2 \pi r dr = @f]
 * @f[      = \rho * (\pi rmax^2 - \pi a_0^2) D_{av}(a_0,rmax) @f]
 *
 * because:
 *
 * @f[ Dav(r1,r2) = 1/ (\pi r_2^2 - \pi r_1^2) * \int_{r_1}^{r_2} D(r) 2 \pi r dr @f]
 *
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density @f$\rho@f$ [kg/m^3]
 * @param[in] alpha                    parameter of ER model
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double AT_RDD_Katz_PowerLawER_dEdx_J_m(  const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double alpha,
    const double Katz_point_coeff_Gy);

/* --------------------------------------------------- SITE RDD ---------------------------------------------------*/

// TODO check why in CPPSC KatzSite and Edmund are not reproducing dose

/**
 * Calculates Site RDD, which is LET-normalized
 * for "old" Katz RDD (derived from linear (on wmax) ER model).
 *
 * @f[ Dsite(r) = 1 / (\rho \pi a_0^2) * (LET - dEdx)   for r < a_0 @f]
 * @f[ Dsite(r) = D(r)                                  for r >= a_0 @f]
 *
 * @param[in] r_m                      TODO
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density @f$\rho@f$ [kg/m^3]
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] dEdx_J_m                 dEdx [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double AT_RDD_Katz_LinearER_DSite_Gy( const double r_m,
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
 * @f[ Dsite(r) = 1 / (\rho \pi a_0^2) * (LET - dEdx)   for r < a_0 @f]
 * @f[ Dsite(r) = D(r)                                  for r >= a_0 @f]
 *
 * @param[in] r_m                      TODO
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] material_density_kg_m3   material density @f$\rho@f$ [kg/m^3]
 * @param[in] alpha                    parameter of ER model
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] dEdx_J_m                 dEdx [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return dEdx [J/m] energy delivered to shell between radius a_0 and r_max
 */
double AT_RDD_Katz_PowerLawER_DSite_Gy( const double r_m,
    const double a0_m,
    const double max_electron_range_m,
    const double material_density_kg_m3,
    const double alpha,
    const double LET_J_m,
    const double dEdx_J_m,
    const double Katz_point_coeff_Gy);

/**
 * TODO
 * @param[in] r_m
 * @param[in] r_min_m                  minimum radius cut-off distance [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] a0_m                     inner radius a0 (lower integration limit) [m]
 * @param[in] er_model                 delta electron range model code number
 * @param[in] alpha                    parameter of ER model
 * @param[in] density_kg_m3            material density @f$\rho@f$ [kg/m^3]
 * @param[in] LET_J_m                  LET [J/m]
 * @param[in] dEdx_J_m                 dEdx [J/m]
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return
 */
 double AT_RDD_KatzSite_Gy( const double r_m,
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
 * @param[in] D_Gy
 * @param[in] r_min_m                  minimum radius cut-off distance [m]
 * @param[in] max_electron_range_m     delta electron maximum range rmax [m]
 * @param[in] a0_m                     radius target [m]
 * @param[in] er_model                 delta electron range model code number
 * @param[in] alpha                    parameter of ER model
 * @param[in] d_max_Gy                 TODO
 * @param[in] Katz_point_coeff_Gy      precalculated coefficient [Gy]
 * @return
 */
double AT_inverse_RDD_KatzSite_m( const double D_Gy,
    const double r_min_m,
    const double max_electron_range_m,
    const double a0_m,
    const long   er_model,
    const double alpha,
    const double d_max_Gy,
    const double Katz_point_coeff_Gy);


#endif /* AT_RDD_SHELLAVERAGED_H_ */
