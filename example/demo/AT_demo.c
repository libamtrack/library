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

    const double E_MeV = 150.0;

    const double E_MeV_u = AT_E_MeV_u_from_E_MeV(E_MeV, PARTICLE_PROTON_NUMBER);
    printf("Proton energy %4.4f [MeV] corresponds to %4.4f [MeV]/u\n", E_MeV, E_MeV_u);

    const double beta = AT_beta_from_E_single(E_MeV_u);
    printf("Relative speed of proton with energy %4.2f [MeV] is equal %1.3f\n", E_MeV, beta);

    const double z_cm = 10;
    const double fluence_cm2 = 1e9;
    const double sigma_E_MeV = 1.5;
    const long material_no = 1;
    const double eps = 0.02;

    double dose_Gy = AT_dose_Bortfeld_Gy_single(z_cm,
                                                E_MeV,
                                                fluence_cm2,
                                                sigma_E_MeV,
                                                material_no,
                                                eps);

    printf("Dose: %4.3f [Gy] (at fluence %g [1/cm2])\n", dose_Gy, fluence_cm2);

    double LETt_keV_um = AT_LET_t_Wilkens_keV_um_single(z_cm,
                                                        E_MeV,
                                                        sigma_E_MeV,
                                                        material_no);

    double LETd_keV_um = AT_LET_d_Wilkens_keV_um_single(z_cm,
                                                        E_MeV,
                                                        sigma_E_MeV,
                                                        material_no);

    printf("LETt: %4.3f [keV/um]\n", LETt_keV_um);
    printf("LETd: %4.3f [keV/um]\n", LETd_keV_um);

    double max_location_cm = AT_max_location_Bortfeld_cm(E_MeV,
                                                         sigma_E_MeV,
                                                         material_no,
                                                         eps);

    printf("maximum dose located at: %g [cm]\n", max_location_cm);

    const double dose_drop_factor = 0.8;
    const double range_cm = AT_range_Bortfeld_cm(E_MeV,
                                                 sigma_E_MeV,
                                                 material_no,
                                                 eps,
                                                 dose_drop_factor,
                                                 1);

    printf("beam range (at %3.2f of max dose): %g cm\n", dose_drop_factor, range_cm);

    const double fwhm_cm = AT_fwhm_Bortfeld_cm(E_MeV,
                                               sigma_E_MeV,
                                               material_no,
                                               eps);

    printf("FWHM: %g [cm]\n", fwhm_cm);

    const double energy_MeV = AT_energy_Bortfeld_MeV(range_cm,
                                                     sigma_E_MeV,
                                                     material_no,
                                                     eps,
                                                     0.8);

    printf("energy calculated from range: %g [MeV], original energy: %g [MeV]\n", energy_MeV, E_MeV);

    double fit_E_MeV;
    double fit_sigma_E_MeV;
    double fit_eps;

    const double measured_range_cm = 15.799712435;
    const double measured_fwhm_cm = 1.5986288328329863;
    const double measured_max_to_plateau = 5.198384978576822;
    const double measured_dose_drop_factor = 0.9;

    AT_fit_Bortfeld(measured_range_cm,
                    measured_fwhm_cm,
                    measured_max_to_plateau,
                    material_no,
                    measured_dose_drop_factor,
                    &fit_E_MeV,
                    &fit_sigma_E_MeV,
                    &fit_eps);
    printf("measured range = %g [cm], FWHM = %g [cm], max/plateau = %g\n", measured_range_cm, measured_fwhm_cm,
           measured_max_to_plateau);
    printf("fitted E = %g [MeV], deltaE = %g [MeV], eps = %g\n", fit_E_MeV, fit_sigma_E_MeV, fit_eps);

    return EXIT_SUCCESS;
}
