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
#include "AT_RDD.h"
#include "AT_KatzModel_Implementation.h"

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
    const long material_no = Water_Liquid;
    const double eps = 0.02;

    const double er_model = ER_Tabata;
    const double el_Rex_m = AT_max_electron_range_m( E_MeV_u, material_no, er_model);
    printf("Max range of delta-ray emitted by ion with energy %4.2f [MeV] is equal %g [mm]\n", E_MeV, el_Rex_m * 1e3);

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

    const int N = 2;
    double r_m_tab[] = {1e-13, 1e-8};
    double rdd_E_MeV_u = 150.0;
    long particle_no = AT_particle_no_from_particle_name_single("1H");
    long rdd_model = 6;
    double KatzPoint_r_min_m = 1e-10;
    double a0_m = 1e-8;
    double d_min_Gy = 1e-80;
    double rdd_parameter[] = {KatzPoint_r_min_m, a0_m, d_min_Gy, 0.0};
    long stopping_power_source_no = 2;
    double D_RDD_Gy[2] = {0.0};

    long res = AT_D_RDD_Gy( N,
                            r_m_tab,
                            rdd_E_MeV_u,
                            particle_no,
                            material_no,
                            rdd_model,
                            rdd_parameter,
                            ER_Waligorski,
                            stopping_power_source_no,
                            D_RDD_Gy);

    for(int i=0; i<N; ++i)
    {
        printf("%e : %e\n", r_m_tab[i], D_RDD_Gy[i]);
    }

    const double RBE_E_MeV_u = 10.0;
    particle_no = AT_particle_no_from_particle_name_single("12C");
    const double D0_Gy = 1.1;
    const double m = 2.36;
    const double sigma0_m2 = 140.8 * (1e-6) * (1e-6);
    const bool use_approximation = true;
    const double kappa = 1204.;
    const double survival = 0.1;

    double rbe = AT_KatzModel_single_field_rbe(RBE_E_MeV_u,
                                               particle_no,
                                               RDD_KatzExtTarget,
                                               rdd_parameter,
                                               ER_Waligorski,
                                               D0_Gy,
                                               m,
                                               sigma0_m2,
                                               use_approximation,
                                               kappa,
                                               stopping_power_source_no,
                                               survival);
    printf("RBE = %g\n", rbe);
    return EXIT_SUCCESS;
}
