/*
 *    AT_EnergyLoss.c
 *    ===============
 *
 *    Copyright 2006, 2014 The libamtrack team
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

#include "AT_EnergyLoss.h"

double AT_xi_keV(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {

    const double Z = AT_average_Z_from_material_no(material_no);
    const double A = AT_average_A_from_material_no(material_no);
    const double density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);

    const double beta = AT_beta_from_E_single(E_MeV_u);
    double z = AT_Z_from_particle_no_single(particle_no);

    return (153.5 * (Z / A) * (z * z) / (beta * beta) * density_g_cm3 * (slab_thickness_um / 1e4));

}

double AT_mean_energy_loss_keV( const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {
    
    double stopping_power_keV_um = 0.0;
    AT_Stopping_Power( "PSTAR",
		1,
		&E_MeV_u,
		&particle_no,
		material_no,
		&stopping_power_keV_um);
    return ( stopping_power_keV_um * slab_thickness_um);
}

double AT_kappa_single(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {
    return ( AT_xi_keV(E_MeV_u, particle_no, material_no, slab_thickness_um) /
            (AT_max_E_transfer_MeV_single(E_MeV_u)*1000));
}

void AT_kappa_multi(const long n,
        const double E_MeV_u[],
        const long particle_no[],
        const long material_no,
        const double slab_thickness_um[],
        double kappa[]) {
    long i;
    for (i = 0; i < n; i++) {
        kappa[i] = AT_kappa_single(E_MeV_u[i], particle_no[i], material_no, slab_thickness_um[i]);
    }
}

void AT_Landau_PDF( const long n, 
        const double lambda_landau[], 
        double density[] ) {
    for (int i = 0; i < n; i++) {
        density[i] = CL_denlan(lambda_landau[i]);
    }
}

void AT_Landau_IDF( const long n, 
        const double rnd[], 
        double lambda_landau[] ) {
    for (int i = 0; i < n; i++) {
        lambda_landau[i] = CL_ranlan(rnd[i]);
    }
}

void AT_lambda_mean_multi(const long n,
        const double E_MeV_u[],
        const long particle_no[],
        const long material_no,
        const double slab_thickness_um[],
        double lambda_mean[]) {
    long i;
    for (i = 0; i < n; i++) {
        lambda_mean[i] = AT_lambda_mean_single(E_MeV_u[i], particle_no[i], material_no, slab_thickness_um[i]);
    }

}

double AT_lambda_mean_single(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {
    double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double beta = AT_beta_from_E_single(E_MeV_u);
    return (-0.42278433509 - beta * beta - log(kappa));
}

void AT_lambda_max_multi(const long n,
        const double E_MeV_u[],
        const long particle_no[],
        const long material_no,
        const double slab_thickness_um[],
        double lambda_max[]) {
    long i;
    for (i = 0; i < n; i++) {
        double lambda_mean = AT_lambda_mean_single(E_MeV_u[i], particle_no[i], material_no, slab_thickness_um[i]);
        lambda_max[i] = AT_lambda_max_single(lambda_mean);
    }
}

double AT_lambda_max_single(double lambda_mean) {
    return (0.60715 + 1.1934 * lambda_mean + (0.67794 + 0.052382 * lambda_mean) * exp(0.94753 + 0.74442 * lambda_mean));
}

double AT_lambda_landau_from_energy_loss_single(const double energy_loss_keV,
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {

    double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double beta = AT_beta_from_E_single(E_MeV_u);
    double xi_keV = AT_xi_keV(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double mean_loss_keV = AT_mean_energy_loss_keV(E_MeV_u,
            particle_no,
            material_no,
            slab_thickness_um);

    return ( (energy_loss_keV - mean_loss_keV) / xi_keV - 0.42278433509 - beta * beta - log(kappa));
}

void AT_lambda_landau_from_energy_loss_multi(const long n,
        const double energy_loss_keV[],
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um,
        double lambda_landau[]) {

    for (long i = 0; i < n; i++) {
        lambda_landau[i] = AT_lambda_landau_from_energy_loss_single(energy_loss_keV[i], E_MeV_u, particle_no, material_no, slab_thickness_um);
    }

    return;
}

double AT_energy_loss_from_lambda_landau_single(const double lambda_landau,
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {

    double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double beta = AT_beta_from_E_single(E_MeV_u);
    double xi_keV = AT_xi_keV(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double mean_loss_keV = AT_mean_energy_loss_keV(E_MeV_u,
            particle_no,
            material_no,
            slab_thickness_um);

    return (xi_keV * (lambda_landau + 0.42278433509 + beta * beta + log(kappa)) + mean_loss_keV);
}

void AT_energy_loss_from_lambda_landau_multi(const long n,
        const double lambda_landau[],
        const double E_MeV_u[],
        const long particle_no[],
        const long material_no,
        const double slab_thickness_um[],
        double energy_loss_keV[]) {


    for (long i = 0; i < n; i++) {
        energy_loss_keV[i] = AT_energy_loss_from_lambda_landau_single(lambda_landau[i],
                E_MeV_u[i],
                particle_no[i],
                material_no,
                slab_thickness_um[i]);
    }

    return;
}

double AT_energy_loss_from_lambda_vavilov_single(const double lambda_vavilov,
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {

    double kappa         = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double beta          = AT_beta_from_E_single(E_MeV_u);
    double xi_keV        = AT_xi_keV(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double mean_loss_keV = AT_mean_energy_loss_keV(E_MeV_u,
            particle_no,
            material_no,
            slab_thickness_um);

    return (xi_keV * (lambda_vavilov / kappa + 0.42278433509 + beta * beta) + mean_loss_keV);
}

void AT_energy_loss_from_lambda_vavilov_multi(const long n,
        const double lambda_vavilov[],
        const double E_MeV_u[],
        const long particle_no[],
        const long material_no,
        const double slab_thickness_um[],
        double energy_loss_keV[]) {


    for (long i = 0; i < n; i++) {
        energy_loss_keV[i] = AT_energy_loss_from_lambda_vavilov_single(lambda_vavilov[i],
                E_MeV_u[i],
                particle_no[i],
                material_no,
                slab_thickness_um[i]);
    }

    return;
}

void AT_Landau_energy_loss_distribution(const long n,
        const double energy_loss_keV[],
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um,
        double fDdD[]) {
    /*
            double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
            double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
            double* lambda        = (double*)calloc(n, sizeof(double));

            AT_lambda_from_energy_loss_multi( n,
                            energy_loss_keV,
                            E_MeV_u,
                            particle_no,
                            material_no,
                            slab_thickness_um,
                            lambda);

            AT_Landau_PDF( n,
                            lambda,
                            fDdD);

            for(int i = 0; i < n; i++){
                    fDdD[i] /= xi;
            }
     */
    return;
}

void AT_Vavilov_PDF(const long n, const double lambda_vavilov[], const double kappa, const double beta,
        double density[]) {
    double beta2 = beta * beta;
    CL_vavset(kappa, beta2);

    for (int i = 0; i < n; i++) {
        density[i] = CL_vavden(lambda_vavilov[i]);
    }
}

void AT_Vavilov_IDF(const long n, const double rnd[], const double kappa[], const double beta[],
        double lambda_vavilov[]) {

    for (int i = 0; i < n; i++) {
        lambda_vavilov[i] = CL_vavran(kappa[i], beta[i]*beta[i], rnd[i]);
    }
}

double AT_lambda_vavilov_from_energy_loss_single(const double energy_loss_keV,
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {

    double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double beta = AT_beta_from_E_single(E_MeV_u);
    double xi_keV = AT_xi_keV(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double mean_loss_keV = AT_mean_energy_loss_keV(E_MeV_u,
            particle_no,
            material_no,
            slab_thickness_um);

    return ( kappa * ((energy_loss_keV - mean_loss_keV) / xi_keV - 0.42278433509 - beta * beta));
}

void AT_lambda_vavilov_from_energy_loss_multi(const long n,
        const double energy_loss_keV[],
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um,
        double lambda_vavilov[]) {

    for (long i = 0; i < n; i++) {
        lambda_vavilov[i] = AT_lambda_vavilov_from_energy_loss_single(energy_loss_keV[i],
                E_MeV_u,
                particle_no,
                material_no,
                slab_thickness_um);
    }

    return;
}

void AT_Vavilov_energy_loss_distribution(const long n,
        const double energy_loss_keV[],
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um,
        double fDdD[]) {
    /*
            double 	kappa	      = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
            double 	beta          = AT_beta_from_E_single(E_MeV_u);
            double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
            double* lambda        = (double*)calloc(n, sizeof(double));

            AT_lambda_from_energy_loss( n,
                            energy_loss_keV,
                            E_MeV_u,
                            particle_no,
                            material_no,
                            slab_thickness_um,
                            lambda);

            AT_Vavilov_PDF( n,
                            lambda,
                        kappa,
                        beta,
                            fDdD);

            for(int i = 0; i < n; i++){
                    fDdD[i] /= xi;
            }
     */
    return;
}

void AT_energy_loss_distribution(const long n,
        const double energy_loss_keV[],
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um,
        double fDdD[]) {

    double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    if (kappa <= 0.01) {
        AT_Landau_energy_loss_distribution(n,
                energy_loss_keV,
                E_MeV_u,
                particle_no,
                material_no,
                slab_thickness_um,
                fDdD);
        return;
    }
    if (kappa < 10) {
        AT_Vavilov_energy_loss_distribution(n,
                energy_loss_keV,
                E_MeV_u,
                particle_no,
                material_no,
                slab_thickness_um,
                fDdD);
        return;
    }
}

/**
 * This is code from ROOT (CERN) 5.3.1
 */
double AT_lambda_Vavilov_Mode(const double kappa, const double beta) {

    CL_vavset(kappa, beta);
    double x = -4.22784335098467134e-01 - log(kappa) - beta*beta;
    if (x>-0.223172) x = -0.223172;
    double eps = 0.01;
    double dx;

    do {
        double p0 = CL_vavden(x - eps);
        double p1 = CL_vavden(x);
        double p2 = CL_vavden(x + eps);
        double y1 = 0.5 * (p2 - p0) / eps;
        double y2 = (p2 - 2 * p1 + p0) / (eps * eps);
        if (y2 != 0) {
            dx = -y1 / y2;
        } else {
            dx = 0.0;
        }
        x += dx;
        if (fabs(dx) < eps) eps = 0.1 * fabs(dx);
    } while (fabs(dx) > 1E-5);
    return x;
}

double AT_lambda_Vavilov_FWHM_left(const double kappa, const double beta) {

    double x = AT_lambda_Vavilov_Mode(kappa, beta);
    CL_vavset(kappa, beta);
    double p = CL_vavden(x) * 0.5;

    x -= 1.3637;
    double eps = 0.01;
    double dx;

    do {
        double p0 = CL_vavden(x);
        double p1 = CL_vavden(x - eps);
        double p2 = CL_vavden(x + eps);
        double y1 = p0 - p;
        double y2 = 0.5 * (p2 - p1) / eps;
        if (y2 != 0) {
            dx = -y1 / y2;
        } else {
            dx = 0.0;
        }
        x += dx;
        if (fabs(dx) < eps) eps = 0.1 * fabs(dx);
    } while (fabs(dx) > 1E-5);
    return x;
}

double AT_lambda_Vavilov_FWHM_right(const double kappa, const double beta) {

    double x = AT_lambda_Vavilov_Mode(kappa, beta);
    CL_vavset(kappa, beta);
    double p = CL_vavden(x) * 0.5;

    x += 2.655;
    double eps = 0.01;
    double dx;

    do {
        double p0 = CL_vavden(x);
        double p1 = CL_vavden(x - eps);
        double p2 = CL_vavden(x + eps);
        double y1 = p0 - p;
        double y2 = 0.5 * (p2 - p1) / eps;
        if (y2 != 0) {
            dx = -y1 / y2;
        } else {
            dx = 0.0;
        }
        x += dx;
        if (fabs(dx) < eps) eps = 0.1 * fabs(dx);
    } while (fabs(dx) > 1E-5);
    return x;
}

double AT_lambda_Vavilov_FWHM(const double kappa, const double beta) {
    return (AT_lambda_Vavilov_FWHM_right(kappa, beta) - AT_lambda_Landau_FWHM_left(kappa, beta));
}

double AT_energy_loss_keV_Vavilov_FWHM(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {
    /*
            double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
            double beta  = AT_beta_from_E_single(E_MeV_u);

            return AT_energy_loss_from_lambda_single(AT_lambda_Vavilov_FWHM_right(kappa, beta),
                            E_MeV_u,
                            particle_no,
                            material_no,
                            slab_thickness_um) -
                            AT_energy_loss_from_lambda_single(AT_lambda_Vavilov_FWHM_left(kappa, beta),
                                                    E_MeV_u,
                                                    particle_no,
                                                    material_no,
                                                    slab_thickness_um);
     */
    return 0.0;
}

double AT_energy_loss_keV_Vavilov_Mode(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {
    /*
            double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
            double beta  = AT_beta_from_E_single(E_MeV_u);

            return AT_energy_loss_from_lambda_single(AT_lambda_Vavilov_Mode(kappa, beta),
                            E_MeV_u,
                            particle_no,
                            material_no,
                            slab_thickness_um);
     */
    return 0.0;
}

double AT_lambda_Vavilov_Mean(const double kappa, const double beta) {
    return -0.42278433509 - 1.0 - log(kappa) - beta*beta;
}

double AT_lambda_Vavilov_Variance(const double kappa, const double beta) {
    return (1.0 - 0.5 * beta * beta) / kappa;
}

double AT_Vavilov_FWHM() {
    return 0.0;
}

double AT_lambda_Vavilov_Skewness(const double kappa, const double beta) {
    return (0.5 - (beta * beta) / 3) / (kappa * kappa) * pow((1 - 0.5 * beta * beta) / kappa, -1.5);
}

double AT_lambda_Landau_Mode() {
    return -0.2258;
}

double AT_energy_loss_keV_Landau_Mode(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {
    /*
            return AT_energy_loss_from_lambda_single(AT_lambda_Landau_Mode(),
                            E_MeV_u,
                            particle_no,
                            material_no,
                            slab_thickness_um);
     */
    return 0.0;
}

double AT_lambda_Landau_FWHM_left() {
    return -0.2258 - 1.3637;
}

double AT_lambda_Landau_FWHM_right() {
    return -0.2258 + 2.655;
}

double AT_lambda_Landau_FWHM() {
    return (AT_lambda_Landau_FWHM_right() - AT_lambda_Landau_FWHM_left());
}

double AT_energy_loss_keV_Landau_FWHM(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {
    /*	double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
            double  xi            = kappa * AT_max_E_transfer_MeV_single(E_MeV_u);
            return(4.02 * xi * 1000);
     */
    return 0.0;
}

double AT_lambda_Landau_Mean(const double kappa, const double beta) {
    return -0.42278433509 - log(kappa) - beta*beta;
}

//////////////////////////////
// GAUSS
//////////////////////////////

void AT_Gauss_PDF(const long n, const double lambda_gauss[], double density[]) {
    for (int i = 0; i < n; i++) {
        density[i] = gsl_ran_ugaussian_pdf(lambda_gauss[i]);
    }
}

void AT_Gauss_IDF(const long n, const double rnd[], double lambda_gauss[]) {
    for (int i = 0; i < n; i++) {
        lambda_gauss[i] = gsl_cdf_ugaussian_Pinv(rnd[i]);
    }
}

double AT_energy_loss_from_lambda_gauss_single(const double lambda_gauss,
        const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {

    double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double beta = AT_beta_from_E_single(E_MeV_u);
    double xi_keV = AT_xi_keV(E_MeV_u, particle_no, material_no, slab_thickness_um);
    double mean_loss_keV = AT_mean_energy_loss_keV(E_MeV_u,
            particle_no,
            material_no,
            slab_thickness_um);

    return (lambda_gauss * sqrt(xi_keV * xi_keV / kappa * (1 - beta * beta / 2)) + mean_loss_keV);
}

void AT_energy_loss_from_lambda_gauss_multi(const long n,
        const double lambda_gauss[],
        const double E_MeV_u[],
        const long particle_no[],
        const long material_no,
        const double slab_thickness_um[],
        double energy_loss_keV[]) {


    for (long i = 0; i < n; i++) {
        energy_loss_keV[i] = AT_energy_loss_from_lambda_gauss_single(lambda_gauss[i],
                E_MeV_u[i],
                particle_no[i],
                material_no,
                slab_thickness_um[i]);
    }

    return;
}

double AT_Gauss_Mode() {
    return 0;
}

double AT_Gauss_Mean() {
    return 0;
}

double AT_Gauss_FWHM() {
    return 0;
}


//////////////////////////////
// GENERAL
//////////////////////////////

double AT_energy_loss_mode(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {

    double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    if (kappa <= 0.01) {
        return AT_energy_loss_keV_Landau_Mode(E_MeV_u,
                particle_no,
                material_no,
                slab_thickness_um);
    }
    if (kappa < 10) {
        return AT_energy_loss_keV_Vavilov_Mode(E_MeV_u,
                particle_no,
                material_no,
                slab_thickness_um);
    }
    // kappa >= 10
    return AT_mean_energy_loss_keV(E_MeV_u,
            particle_no,
            material_no,
            slab_thickness_um);
}

double AT_energy_loss_FWHM(const double E_MeV_u,
        const long particle_no,
        const long material_no,
        const double slab_thickness_um) {

    double kappa = AT_kappa_single(E_MeV_u, particle_no, material_no, slab_thickness_um);
    if (kappa <= 0.01) {
        return AT_energy_loss_keV_Landau_FWHM(E_MeV_u,
                particle_no,
                material_no,
                slab_thickness_um);
    }
    if (kappa < 10) {
        return AT_energy_loss_keV_Vavilov_FWHM(E_MeV_u,
                particle_no,
                material_no,
                slab_thickness_um);
    }
    //	// kappa >= 10
    //	return 	AT_Bethe_mean_energy_loss_MeV( E_MeV_u,
    //			particle_no,
    //			material_no,
    //			slab_thickness_um);
    return 0.0;
}


