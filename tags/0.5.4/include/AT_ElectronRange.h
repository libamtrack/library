#ifndef AT_ELECTRONRANGE_H_
#define AT_ELECTRONRANGE_H_

/**
 * @brief Electron range models
 */

/*
 *    AT_ElectronRange.h
 *    ==============
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

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "gsl/gsl_pow_int.h"
#include "gsl/gsl_sf_log.h"

#include "AT_NumericalRoutines.h"
#include "AT_DataMaterial.h"
#include "AT_DataParticle.h"
#include "AT_PhysicsRoutines.h"
#include "AT_Error.h"

/**
 * Electron range models code numbers
 */
enum AT_ERModels{
  ER_Test                  = 1,     /**< dummy electron range models */
      ER_ButtsKatz         = 2,     /**< Butts&Katz(?) electron range model, R = k * w, valid for ?? < w < 2keV , TODO ref needed [Katz et al., 1972] */
      ER_Waligorski        = 3,     /**< Waligorski(?) electron range model, R = k * w ^ alpha, valid for ?? < w < ? TODO ref needed*/
      ER_Geiss             = 4,     /**< Geiss(?) electron range model, R = k * E ^ alpha, valid for ?? < w < ? TODO ref needed [Geiss, 1997]*/
      ER_Scholz            = 5,     /**< Scholz(?) electron range model, R = k * E ^ alpha, valid for ?? < w < ? TODO ref needed [Scholz, 2001]*/
      ER_Edmund            = 6,     /**< Edmund(?) electron range model, R = k * w ^ alpha, valid for ?? < w < ? TODO ref needed*/
      ER_Tabata            = 7      /**< Tabata electron range model, valid for 0.3keV < w < 30MeV TODO ref needed [Tabata, 1972]*/
};

#define ER_DATA_N    7


/**
 * TODO
 */
typedef struct {
  int     n;
  int     ER_no[ER_DATA_N];
  char*   ER_name[ER_DATA_N];
} AT_ER_data_struct;


/**
 * TODO
 */
static const AT_ER_data_struct AT_ER_Data = {
    ER_DATA_N,
    {  ER_Test,                 ER_ButtsKatz,                       ER_Waligorski,                             ER_Geiss,                         ER_Scholz,                       ER_Edmund,                          ER_Tabata },
    {  "simple test ER model",  "Butts & Katz' ER model (linear)",  "Waligorski's ER model (power-law wmax)",  "Geiss' ER model (power-law E)", "Scholz' ER model (power-law E)", "Edmund' ER model (power-law wmax)","Tabata  ER model"}
};


/**
 * Returns name of the electron model from index
 *
 * @param[in]   ER_no    electron-range-model index
 * @param[out]  ER_name  string containing the electron-range model name
 * @return      Status code
 */
int getERName( const int ER_no,
    char* ER_name);


/**
 * 1e-5 * wmax_keV
 * @param[in] wmax_keV
 * @return
 */
double AT_ER_ButtsKatz_range_g_cm2( double wmax_keV );


/**
 * 6e-6 * pow( wmax_keV, alpha )
 * @param[in] wmax_keV
 * @return
 */
double AT_ER_Waligorski_range_g_cm2( double wmax_keV );


/**
 * 6.13*1e-6  * pow( wmax_keV, alpha )
 * @param[in] wmax_keV
 * @return
 */
double AT_ER_Edmund_range_g_cm2( double wmax_keV );


/**
 * 4e-5 * pow(E_MeV_u, 1.5)
 * @param[in] E_MeV_u
 * @return
 */
double AT_ER_Geiss_range_g_cm2( double E_MeV_u );


/**
 * 5e-5 * pow(E_MeV_u, 1.7)
 * @param[in] E_MeV_u
 * @return
 */
double AT_ER_Scholz_range_g_cm2( double E_MeV_u );


/**
 *  tau = 2.0 * gsl_pow_2(beta) / (1.0 - gsl_pow_2(beta))
 * (a1_g_cm2)*(((gsl_sf_log(1.0 + a2 * tau))/a2) - ((a3*tau)/(1.0 + a4*pow(tau,a5))) )
 *
 * @param[in] beta
 * @param[in] a1_g_cm2
 * @param[in] a2
 * @param[in] a3
 * @param[in] a4
 * @param[in] a5
 * @return range
 */
double AT_ER_Tabata_range_g_cm2( double beta,
    double a1_g_cm2,
    double a2,
    double a3,
    double a4,
    double a5);


/**
 * Alpha exponent in power law ER models. It is defined as follows:
 * alpha = 1.667   when energy of ejected delta electron wmax is higher than 1keV
 * alpha = 1.079   when energy of ejected delta electron wmax is less or equal than 1keV
 *
 * @param[in] E_MeV_u  kinetic energy for particle [in MeV/u]
 * @return
 */
double AT_ER_PowerLaw_alpha( const double E_MeV_u );


/**
 * Tabata model constants
 * @param[in] average_A
 * @param[in] average_Z
 * @param[out] a1_g_cm2
 * @param[out] a2
 * @param[out] a3
 * @param[out] a4
 * @param[out] a5
 */
void AT_ER_Tabata_constants( const double average_A,
    const double average_Z,
    double* a1_g_cm2,
    double* a2,
    double* a3,
    double* a4,
    double* a5);


/**
 * Returns the maximum electron range (track radius) in m
 * for a given parametrization
 *
 * @param[in]  number_of_particles          number of particles in the incident field
 * @param[in]  E_MeV_u                      kinetic energy for particles in the given field (array of size number_of_particles)
 * @param[in]  material_no                  material index
 * @param[in]  er_model                     electron-range model index
 * @param[out] max_electron_range_m         electron range (track radius) in m  (array of size number_of_particles)
 */
void AT_max_electron_ranges_m( const long number_of_particles ,
    const double  E_MeV_u[],
    const int     material_no,
    const int     er_model,
    double        max_electron_range_m[]);


/**
 * Returns the maximum electron range (track diameter) in m for given energy
 *
 * @param[in]  E_MeV_u                      kinetic energy for particles in the given field
 * @param[in]  material_no                  index for detector material
 * @param[in]  er_model                     index for electron-range model chosen
 * @return                                  electron range (track diameter) in m
 */
double AT_max_electron_range_m(  const double E_MeV_u,
    const int    material_no,
    const int    er_model);


#endif /* AT_ELECTRONRANGE_H_ */
