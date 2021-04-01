#ifndef AMTRACK_AT_KATZMODEL_H
#define AMTRACK_AT_KATZMODEL_H

/**
 * @brief Katz model algorithm implementation
 */

/*
 *    AT_KatzModel.h
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

#include "AT_KatzModel_Implementation.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/**
 * Katz model flavours
 */
enum KatzModelFlavour {
    Katz_Linear = 1,      /**< based on Butts-Katz linear ER model */
    Katz_PowerLaw = 2,      /**< based on Zhang RDD and power law ER model */
    Katz_Cucinotta = 3       /**< based on Cucinotta RDD and Tabata ER model (no approximated version!) */
};


/**
 * Calculates inactivation cross-section (sigma) for given energy.
 * Sigma is evaluated by double integration:
 *  - first of radial dose distribution (over site with radius a0)
 *  - second of inactivation probability (over all radii)
 *  Sigma can be calculated using 3 possible combination of radial dose and electron range models
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] m
 * @param[in] D0_Gy
 * @param[in] a0_um
 * @param[in] katz_model_flavour
 * @param[in] stop_power_source
 * @return inactivation cross-section (sigma) [um2]
 */
double AT_KatzModel_sigma_um2_single(
        const double E_MeV_u,
        const long particle_no,
        const double m,
        const double D0_Gy,
        const double a0_um,
        const long katz_model_flavour,
        const long stop_power_source);

/**
 * Calculates inactivation cross-section (sigma) for given energy.
 * Vectorised version of AT_KatzModel_sigma_um2_single
 * @param[in] n
 * @param[in] E_MeV_u  (array of size n)
 * @param[in] particle_no
 * @param[in] m
 * @param[in] D0_Gy
 * @param[in] a0_um
 * @param[in] katz_model_flavour
 * @param[in] stop_power_source
 * @param[out] sigma_um2 (array of size n)
 * @return status code
 */
int AT_KatzModel_sigma_um2(
        const long n,
        const double E_MeV_u[],
        const long particle_no,
        const double m,
        const double D0_Gy,
        const double a0_um,
        const long katz_model_flavour,
        const long stop_power_source,
        double sigma_um2[]);

/**
 * Calculates approximated inactivation cross-section (sigma) for given energy.
 * approximation via "trkwid" function
 * Sigma can be calculated using 2 possible combination of radial dose and electron range models
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] m
 * @param[in] sigma0_um2
 * @param[in] kappa
 * @param[in] katz_model_flavour
 * @return inactivation cross-section (sigma) [um2]
 */
double AT_KatzModel_sigma_approx_um2_single(
        const double E_MeV_u,
        const long particle_no,
        const double m,
        const double sigma0_um2,
        const double kappa,
        const long katz_model_flavour);

/**
 * Calculates approximated inactivation cross-section (sigma) for given energy.
 * Vectorised version of AT_KatzModel_sigma_approx_um2_single
 * @param[in] n
 * @param[in] E_MeV_u (array of size n)
 * @param[in] particle_no
 * @param[in] m
 * @param[in] sigma0_um2
 * @param[in] kappa
 * @param[in] katz_model_flavour
 * @param[out] sigma_um2 (array of size n)
 * @return status code
 */
int AT_KatzModel_sigma_approx_um2(
        const long n,
        const double E_MeV_u[],
        const long particle_no,
        const double m,
        const double sigma0_um2,
        const double kappa,
        const long katz_model_flavour,
        double sigma_um2[]);

/**
 * Calculates cell survival for given dose, energy and particle type.
 * @param[in] dose_Gy
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] m
 * @param[in] D0_Gy
 * @param[in] sigma0_um2
 * @param[in] kappa
 * @param[in] a0_um
 * @param[in] katz_model_flavour
 * @param[in] approximate
 * @param[in] stopping_power_source_no
 * @return cell survival
 */
double AT_KatzModel_survival_single(
        const double dose_Gy,
        const double E_MeV_u,
        const long particle_no,
        const double m,
        const double D0_Gy,
        const double sigma0_um2,
        const double kappa,
        const double a0_um,
        const long katz_model_flavour,
        const bool approximate,
        const long stopping_power_source_no);


/**
 * Calculates cell survival for given dose, energy and particle type.
 * Vectorised version of AT_KatzModel_survival_single
 * @param[in] n
 * @param[in] dose_Gy (array of size n)
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] m
 * @param[in] D0_Gy
 * @param[in] sigma0_um2
 * @param[in] kappa
 * @param[in] a0_um
 * @param[in] katz_model_flavour
 * @param[in] approximate
 * @param[in] stopping_power_source_no
 * @param[out] survival (array of size n)
 * @return status code
 */
int AT_KatzModel_survival(
        const long n,
        const double dose_Gy[],
        const double E_MeV_u,
        const long particle_no,
        const double m,
        const double D0_Gy,
        const double sigma0_um2,
        const double kappa,
        const double a0_um,
        const long katz_model_flavour,
        const bool approximate,
        const long stopping_power_source_no,
        double survival[]);


/**
 * Calculates RBE for given energy and particle type.
 * @param[in] E_MeV_u
 * @param[in] particle_no
 * @param[in] m
 * @param[in] D0_Gy
 * @param[in] sigma0_um2
 * @param[in] kappa
 * @param[in] a0_um
 * @param[in] katz_model_flavour
 * @param[in] approximate
 * @param[in] stopping_power_source_no
 * @param[in] level
 * @return rbe
 */
double AT_KatzModel_RBE_single(
        const double E_MeV_u,
        const long particle_no,
        const double m,
        const double D0_Gy,
        const double sigma0_um2,
        const double kappa,
        const double a0_um,
        const long katz_model_flavour,
        const bool approximate,
        const long stopping_power_source_no,
        const double level);

/**
 * Calculates RBE for given energy and particle type.
 * Vectorised version of AT_KatzModel_RBE_single
 * @param[in] n
 * @param[in] E_MeV_u (array of size n)
 * @param[in] particle_no
 * @param[in] m
 * @param[in] D0_Gy
 * @param[in] sigma0_um2
 * @param[in] kappa
 * @param[in] a0_um
 * @param[in] katz_model_flavour
 * @param[in] approximate
 * @param[in] stopping_power_source_no
 * @param[in] level
 * @param[out] rbe (array of size n)
 * @return status code
 */
int AT_KatzModel_RBE(
        const long n,
        const double E_MeV_u[],
        const long particle_no,
        const double m,
        const double D0_Gy,
        const double sigma0_um2,
        const double kappa,
        const double a0_um,
        const long katz_model_flavour,
        const bool approximate,
        const long stopping_power_source_no,
        const double level,
        double rbe[]);

#endif //AMTRACK_AT_KATZMODEL_H
