/**
 * @brief TODO
 */

/*
 *    AT_KatzModel.c
 *    ===========================
 *
 *    Created on: 01.03.2010
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

#include "AT_KatzModel.h"

double AT_KatzModel_sigma_um2_single(
        const double E_MeV_u,
        const long particle_no,
        const double a0_um,
        const double m,
        const double D0_Gy,
        const long katz_model_flavour,
        const long stop_power_source) {

    double sigma_um2 = -1.; // result

    const double gamma_parameters[5] = {0., D0_Gy, 1., m, 0.}; // S_max, D0_Gy, c, m, None
    long rdd_model = RDD_Test;
    const double a0_m = 1e-6 * a0_um;  // um -> m conversion
    double rdd_parameters[3] = {0., a0_m, 0.};
    long er_model = ER_Test;
    double inactivation_cross_section_m2;

    if (katz_model_flavour == Katz_Linear) {
        rdd_model = RDD_KatzExtTarget;
        er_model = ER_ButtsKatz;
    } else if (katz_model_flavour == Katz_PowerLaw) {
        rdd_model = RDD_KatzExtTarget;
        er_model = ER_Waligorski;
    } else if (katz_model_flavour == Katz_Cucinotta) {
        rdd_model = RDD_CucinottaExtTarget;
        er_model = ER_Tabata;
    } else {
        return sigma_um2;
    }
    rdd_parameters[0] = AT_RDD_Data.parameter_default[AT_RDD_index_from_RDD_number(rdd_model)][0];
    rdd_parameters[2] = AT_RDD_Data.parameter_default[AT_RDD_index_from_RDD_number(rdd_model)][2];

    AT_KatzModel_inactivation_cross_section_m2(
            1,
            &E_MeV_u,
            particle_no,
            rdd_model,
            rdd_parameters,
            er_model,
            gamma_parameters,
            stop_power_source,
            &inactivation_cross_section_m2);

    sigma_um2 = (1e6) * (1e6) * inactivation_cross_section_m2;  // m2 -> um2 conversion
    return sigma_um2;
}


int AT_KatzModel_sigma_um2(
        const long n,
        const double E_MeV_u[],
        const long particle_no,
        const double a0_um,
        const double m,
        const double D0_Gy,
        const long katz_model_flavour,
        const long stop_power_source,
        double sigma_um2[]) {

    long i;
    for (i = 0; i < n; i++) {
        sigma_um2[i] = AT_KatzModel_sigma_um2_single(E_MeV_u[i], particle_no, a0_um, m, D0_Gy, katz_model_flavour,
                                                     stop_power_source);
    }

    return EXIT_SUCCESS;
}


double AT_KatzModel_sigma_approx_um2_single(
        const double E_MeV_u,
        const long particle_no,
        const double kappa,
        const double m,
        const double sigma0_um2,
        const long katz_model_flavour) {
    double sigma_um2 = -1.; // result

    if (katz_model_flavour == Katz_Cucinotta) {
        return -1.0;
    } else {

        double beta = AT_beta_from_E_single(E_MeV_u);
        double zeff = AT_effective_charge_from_beta_single(beta, AT_Z_from_particle_no_single(particle_no));
        double z2kappabeta2 = gsl_pow_2(zeff / beta) / kappa;
        double scaling_envelope = pow(1.0 - exp(-z2kappabeta2), m);

        double factor = 0.0;
        if (scaling_envelope < 0.98) { // grain-count regime
            factor = scaling_envelope;
        } else { // track-width regime
            if (katz_model_flavour == Katz_Linear) {
                factor = AT_KatzModel_KatzExtTarget_ButtsKatz_TrackWidth(z2kappabeta2, m);
            } else if (katz_model_flavour == Katz_PowerLaw) {
                factor = AT_KatzModel_KatzExtTarget_Zhang_TrackWidth(z2kappabeta2, m);
            }
        }
        sigma_um2 = sigma0_um2 * factor;
    }

    return sigma_um2;

}

int AT_KatzModel_sigma_approx_um2(
        const long n,
        const double E_MeV_u[],
        const long particle_no,
        const double kappa,
        const double m,
        const double sigma0_um2,
        const long katz_model_flavour,
        double sigma_um2[]){
    long i;
    for (i = 0; i < n; i++) {
        sigma_um2[i] = AT_KatzModel_sigma_approx_um2_single(E_MeV_u[i], particle_no, kappa, m, sigma0_um2, katz_model_flavour);
    }

    return EXIT_SUCCESS;
}
