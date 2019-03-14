/**
 * @brief Bortfeld Model
 */

/*
 *    AT_NumericalRoutines.c
 *    ==============
 *
 *    Created on: 11.02.2019
 *    Creator: grzanka
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

#include "AT_BortfeldModel.h"
#include "AT_DataMaterial.h"

double AT_dose_Bortfeld_Gy_single(const double z_cm,
                                  const double fluence_cm2,
                                  const double E_MeV_u,
                                  const double sigma_E_MeV_u,
                                  const long material_no,
                                  const double eps) {


    const double beta_cm = 0.012; // slope parameter of fluence reduction relation [1/cm]
    const double gamma = 0.6; // fraction of locally absorbed energy released in non-elastic nuclear interactions

    double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])


    double range_cm = alpha * pow(E_MeV_u, p);  // range in [cm]

    double ni1 = 1.0 / p;
    double ni2 = 1.0 + (1.0 / p);

    // assign default value to sigma being 1% of kinetic energy
    double tmp_sigma_E_MeV_u = sigma_E_MeV_u;
    if (sigma_E_MeV_u < 0)
        tmp_sigma_E_MeV_u = 0.01 * E_MeV_u;

    // assign default value to eps being 0.03 (following Bortfled original paper)
    double tmp_eps = eps;
    if (eps < 0)
        tmp_eps = 0.03;

    double sigma_mono_cm = 0.012 * pow(range_cm, 0.935); // width of Gaussian range straggling
    double sigma_cm = sqrt(sigma_mono_cm * sigma_mono_cm + pow((tmp_sigma_E_MeV_u * alpha * p * pow(E_MeV_u, p - 1)), 2));

    double dose_Gy = fluence_cm2; // returned value

    double factor1 = AT_range_straggling_convolution(z_cm, range_cm, sigma_cm, ni1);

    double factor2 = AT_range_straggling_convolution(z_cm, range_cm, sigma_cm, ni2);
    factor2 *= ((beta_cm / p) + (gamma * beta_cm) + (tmp_eps / range_cm));
    factor2 *= sigma_cm;

    // calculate Gamma(1/p)
    double tmp1;
    AT_gamma_(&ni1, &tmp1);
    factor2 *= tmp1;

    // calculate Gamma(1+1/p)
    double tmp2;
    AT_gamma_(&ni2, &tmp2);
    factor2 /= tmp2;

    dose_Gy /= AT_density_g_cm3_from_material_no(material_no);
    dose_Gy /= (p * pow(alpha, ni1));
    dose_Gy /= (1.0 + beta_cm * range_cm);

    dose_Gy *= (factor1 + factor2);

    dose_Gy *= 1.6021766e-10;  // convert MeV/g into Gy

    return dose_Gy;
}


void AT_dose_Bortfeld_Gy_multi(const long n,
                               const double z_cm[],
                               const double fluence_cm2,
                               const double E_MeV_u,
                               const double sigma_E_MeV_u,
                               const long material_no,
                               const double eps,
                               double dose_Gy[]) {
    long i;
    for (i = 0; i < n; i++) {
        dose_Gy[i] = AT_dose_Bortfeld_Gy_single(z_cm[i], E_MeV_u, fluence_cm2, sigma_E_MeV_u, material_no, eps);
    }
}

double AT_LET_t_Wilkens_keV_um_single(const double z_cm,
                                      const double E_MeV_u,
                                      const double sigma_E_MeV_u,
                                      const long material_no){

    double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])

    double range_cm = alpha * pow(E_MeV_u, p);  // range in [cm]
    double regul_factor_cm = 2 * 1e-4; // regularization factor in [cm] ( 1cm = 1e4 um)

    double ni1 = 1.0 + (1.0 / p);

    // assign default value to sigma being 1% of kinetic energy
    double tmp_sigma_E_MeV_u = sigma_E_MeV_u;
    if (sigma_E_MeV_u < 0)
        tmp_sigma_E_MeV_u = 0.01 * E_MeV_u;

    double sigma_mono_cm = 0.012 * pow(range_cm, 0.935); // width of Gaussian range straggling
    double sigma_cm = sqrt(sigma_mono_cm * sigma_mono_cm + pow((tmp_sigma_E_MeV_u * alpha * p * pow(E_MeV_u, p - 1)), 2));

    double xi = (z_cm - range_cm) / sigma_cm;
    double zeta = (z_cm - range_cm - regul_factor_cm) / sigma_cm;


    double LET_t_keV_um = 1.0; // returned value

    LET_t_keV_um /= (sigma_cm * range_cm * pow(alpha, 1.0 / p));

    double nominator = 0.0;

    nominator += M_SQRTPI * sqrt(2.0) * (AT_range_straggling_convolution(z_cm, range_cm + regul_factor_cm, sigma_cm, ni1) -
            AT_range_straggling_convolution(z_cm, range_cm, sigma_cm, ni1));

    nominator += range_cm * pow(range_cm / 2.0, 1.0 / p) * exp( -(xi +zeta)*(xi+zeta) / 8.0 );

    double denominator = M_SQRTPI * sqrt(2.0) * AT_range_straggling_convolution(z_cm, range_cm, sigma_cm, 1.0);

    LET_t_keV_um *= nominator;
    LET_t_keV_um /= denominator;

    LET_t_keV_um *= 0.1; // from MeV/cm to keV/um

    return LET_t_keV_um;

}


void AT_LET_t_Wilkens_keV_um_multi(const long n,
                                   const double z_cm[],
                                   const double E_MeV_u,
                                   const double sigma_E_MeV_u,
                                   const long material_no,
                                   double LET_keV_um[]){
    long i;
    for (i = 0; i < n; i++) {
        LET_keV_um[i] = AT_LET_t_Wilkens_keV_um_single(z_cm[i], E_MeV_u, sigma_E_MeV_u, material_no);
    }
}