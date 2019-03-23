/**
 * @brief Demo
 */

/*
 *    AT_demo.c
 *    ===================
 *
 *    Created on: 2010-10-06
 *    Creator: kongruencja
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

#include <stdio.h>
#include <stdlib.h>

#include "AT_PhysicsRoutines.h"
#include "AT_ProtonAnalyticalModels.h"
#include "AT_ProtonAnalyticalBeamParameters.h"

int main(int argc, char *argv[]) {

    if (argc != 1) {
        printf("Usage: %s\n", argv[0]);
        return EXIT_FAILURE;
    }

    const double E_MeV_u = 150.0;
    double beta = AT_beta_from_E_single(E_MeV_u);
    printf("Relative speed of particle with energy %4.2f is equal %1.3f\n", E_MeV_u, beta);

    const double z_cm = 10;
    const double fluence_cm2 = 1e9;
    const double sigma_E_MeV_u = 1.0;
    const long material_no = 1;
    const double eps = 0.02;

    double dose_Gy = AT_dose_Bortfeld_Gy_single(z_cm,
                                                E_MeV_u,
                                                fluence_cm2,
                                                sigma_E_MeV_u,
                                                material_no,
                                                eps);

    printf("Dose: %g Gy\n", dose_Gy);

    double LETt_keV_um = AT_LET_t_Wilkens_keV_um_single(z_cm,
                                                        E_MeV_u,
                                                        sigma_E_MeV_u,
                                                        material_no);

    double LETd_keV_um = AT_LET_d_Wilkens_keV_um_single(z_cm,
                                                        E_MeV_u,
                                                        sigma_E_MeV_u,
                                                        material_no);

    printf("LETt: %g keV/um\n", LETt_keV_um);
    printf("LETd: %g keV/um\n", LETd_keV_um);

    double max_location_cm = AT_max_location_Bortfeld_cm(E_MeV_u,
                                                         sigma_E_MeV_u,
                                                         material_no,
                                                         eps);

    printf("max location: %g cm\n", max_location_cm);

    double range_cm = AT_range_Bortfeld_cm(E_MeV_u,
                                           sigma_E_MeV_u,
                                           material_no,
                                           eps,
                                           0.8,
                                           1);

    printf("range: %g cm\n", range_cm);

    double fwhm_cm = AT_fwhm_Bortfeld_cm(E_MeV_u,
                                         sigma_E_MeV_u,
                                         material_no,
                                         eps);

    printf("FWHM: %g cm\n", fwhm_cm);

    double energy_MeV = AT_energy_Bortfeld_MeV_u(range_cm,
                                                 sigma_E_MeV_u,
                                                 material_no,
                                                 eps,
                                                 0.8);

    printf("energy: %g MeV\n", energy_MeV);


    double max_plateau = 4.0;
    double dose_drop = 0.8;

    double fit_E_MeV_u;
    double fit_sigma_E_MeV_u;
    double fit_eps;

    AT_fit_Bortfeld(range_cm,
                    fwhm_cm,
                    max_plateau,
                    material_no,
                    dose_drop,
                    &fit_E_MeV_u,
                    &fit_sigma_E_MeV_u,
                    &fit_eps);

    printf("E = %g, deltaE = %g, eps = %g\n", fit_E_MeV_u, fit_sigma_E_MeV_u, fit_eps);

    AT_fit_Bortfeld(2.9,
                    0.3,
                    4.8,
                    material_no,
                    0.9,
                    &fit_E_MeV_u,
                    &fit_sigma_E_MeV_u,
                    &fit_eps);
    printf("E = %g, deltaE = %g, eps = %g\n", fit_E_MeV_u, fit_sigma_E_MeV_u, fit_eps);

    AT_fit_Bortfeld(3.25,
                    0.54,
                    3.7,
                    material_no,
                    0.9,
                    &fit_E_MeV_u,
                    &fit_sigma_E_MeV_u,
                    &fit_eps);
    printf("E = %g, deltaE = %g, eps = %g\n", fit_E_MeV_u, fit_sigma_E_MeV_u, fit_eps);

    AT_fit_Bortfeld(32.0,
                    2.708,
                    3.76,
                    material_no,
                    0.9,
                    &fit_E_MeV_u,
                    &fit_sigma_E_MeV_u,
                    &fit_eps);
    printf("E = %g, deltaE = %g, eps = %g\n", fit_E_MeV_u, fit_sigma_E_MeV_u, fit_eps);

    AT_fit_Bortfeld(3.2,
                    0.4,
                    4.48,
                    material_no,
                    0.9,
                    &fit_E_MeV_u,
                    &fit_sigma_E_MeV_u,
                    &fit_eps);
    printf("E = %g, deltaE = %g, eps = %g\n", fit_E_MeV_u, fit_sigma_E_MeV_u, fit_eps);

    AT_fit_Bortfeld(15.799712435,
                    1.5986288328329863,
                    5.198384978576822,
                    material_no,
                    0.9,
                    &fit_E_MeV_u,
                    &fit_sigma_E_MeV_u,
                    &fit_eps);
    printf("E = %g, deltaE = %g, eps = %g\n", fit_E_MeV_u, fit_sigma_E_MeV_u, fit_eps);

    return EXIT_SUCCESS;
}
