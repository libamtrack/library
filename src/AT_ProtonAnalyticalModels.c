/**
 * @brief Proton analytical models of dose, LET and RBE
 */

/*
 *    AT_ProtonAnalyticalModels.c
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

#include "AT_ProtonAnalyticalModels.h"

/****************************************** Bortfeld dose model *******************************************************/

double AT_dose_Bortfeld_Gy_single(const double z_cm,
                                  const double fluence_cm2,
                                  const double E_MeV,
                                  const double sigma_E_MeV,
                                  const long material_no,
                                  const double eps) {


    const double beta_cm = 0.012; // slope parameter of fluence reduction relation [1/cm]
    const double gamma = 0.6; // fraction of locally absorbed energy released in non-elastic nuclear interactions

    const double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    const double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])

    const double range_cm = alpha * pow(E_MeV, p);  // range in [cm]

    const double ni1 = 1.0 / p;
    const double ni2 = 1.0 + (1.0 / p);

    // assign default value to sigma being 1% of kinetic energy
    double tmp_sigma_E_MeV = sigma_E_MeV;
    if (sigma_E_MeV < 0)
        tmp_sigma_E_MeV = 0.01 * E_MeV;

    // assign default value to eps being 0.03 (following Bortfled original paper)
    double tmp_eps = eps;
    if (eps < 0)
        tmp_eps = 0.03;

    const double sigma_mono_cm = 0.012 * pow(range_cm, 0.935); // width of Gaussian range straggling
    const double sigma_cm = sqrt(
            sigma_mono_cm * sigma_mono_cm + pow((tmp_sigma_E_MeV * alpha * p * pow(E_MeV, p - 1)), 2));

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
                               const double E_MeV,
                               const double sigma_E_MeV,
                               const long material_no,
                               const double eps,
                               double dose_Gy[]) {
    long i;
    for (i = 0; i < n; i++) {
        dose_Gy[i] = AT_dose_Bortfeld_Gy_single(z_cm[i], fluence_cm2, E_MeV, sigma_E_MeV, material_no, eps);
    }
}

/****************************************** Wilkens LET model *******************************************************/

double AT_LET_t_Wilkens_keV_um_single(const double z_cm,
                                      const double E_MeV,
                                      const double sigma_E_MeV,
                                      const long material_no) {

    double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])

    double range_cm = alpha * pow(E_MeV, p);  // range in [cm]

    double regul_factor_cm = 2.0 * 1e-4; // regularization factor in [cm] ( 1cm = 1e4 um)

    double ni1 = 1.0 + (1.0 / p);

    // assign default value to sigma being 1% of kinetic energy
    double tmp_sigma_E_MeV = sigma_E_MeV;
    if (sigma_E_MeV < 0)
        tmp_sigma_E_MeV = 0.01 * E_MeV;

    double sigma_mono_cm = 0.012 * pow(range_cm, 0.935); // width of Gaussian range straggling
    double sigma_cm = sqrt(
            sigma_mono_cm * sigma_mono_cm + pow((tmp_sigma_E_MeV * alpha * p * pow(E_MeV, p - 1)), 2));

    double xi = (z_cm - range_cm) / sigma_cm;
    double zeta = (z_cm - range_cm - regul_factor_cm) / sigma_cm;


    double LET_t_keV_um = 1.0; // returned value

    LET_t_keV_um /= (sigma_cm * regul_factor_cm * pow(alpha, 1.0 / p));

    double nominator = M_SQRTPI * M_SQRT2 * sigma_cm;

    nominator *= (AT_range_straggling_convolution(z_cm, range_cm + regul_factor_cm, sigma_cm, ni1) -
                  AT_range_straggling_convolution(z_cm, range_cm, sigma_cm, ni1));

    nominator -= regul_factor_cm * pow(regul_factor_cm / 2.0, 1.0 / p) * exp(-(xi + zeta) * (xi + zeta) / 8.0);

    if (nominator < 0)
        nominator = 0.0;

    double denominator = M_SQRTPI * M_SQRT2 * AT_range_straggling_convolution(z_cm, range_cm, sigma_cm, 1.0);

    if (denominator < 0)
        denominator = 0.0;

    LET_t_keV_um *= nominator;
    LET_t_keV_um /= denominator;

    LET_t_keV_um *= 0.1; // from MeV/cm to keV/um

    return LET_t_keV_um;

}


void AT_LET_t_Wilkens_keV_um_multi(const long n,
                                   const double z_cm[],
                                   const double E_MeV,
                                   const double sigma_E_MeV,
                                   const long material_no,
                                   double LET_keV_um[]) {
    long i;
    for (i = 0; i < n; i++) {
        LET_keV_um[i] = AT_LET_t_Wilkens_keV_um_single(z_cm[i], E_MeV, sigma_E_MeV, material_no);
    }
}

double AT_LET_d_Wilkens_keV_um_single(const double z_cm,
                                      const double E_MeV,
                                      const double sigma_E_MeV,
                                      const long material_no) {

    double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])

    double range_cm = alpha * pow(E_MeV, p);  // range in [cm]

    double regul_factor_cm = 2.0 * 1e-4; // regularization factor in [cm] ( 1cm = 1e4 um)

    double ni2 = 1.0 + (1.0 / p);
    double ni3 = 2.0 / p;

    // assign default value to sigma being 1% of kinetic energy
    double tmp_sigma_E_MeV = sigma_E_MeV;
    if (sigma_E_MeV < 0)
        tmp_sigma_E_MeV = 0.01 * E_MeV;

    double sigma_mono_cm = 0.012 * pow(range_cm, 0.935); // width of Gaussian range straggling
    double sigma_cm = sqrt(
            sigma_mono_cm * sigma_mono_cm + pow((tmp_sigma_E_MeV * alpha * p * pow(E_MeV, p - 1)), 2));

    double xi = (z_cm - range_cm) / sigma_cm;
    double zeta = (z_cm - range_cm - regul_factor_cm) / sigma_cm;


    double LET_d_keV_um = 1.0; // returned value

    LET_d_keV_um /= (pow(alpha, 1.0 / p) * p * (2.0 - p));

    double nominator = M_SQRTPI * M_SQRT2 * sigma_cm;

    nominator *= (AT_range_straggling_convolution(z_cm, range_cm + regul_factor_cm, sigma_cm, ni3) -
                  AT_range_straggling_convolution(z_cm, range_cm, sigma_cm, ni3));

    nominator -= 2.0 * regul_factor_cm * pow(regul_factor_cm / 2.0, ni3) * exp(-(xi + zeta) * (xi + zeta) / 8.0);

    if (nominator < 0)
        nominator = 0.0;

    double denominator = M_SQRTPI * M_SQRT2 * sigma_cm;

    denominator *= (AT_range_straggling_convolution(z_cm, range_cm + regul_factor_cm, sigma_cm, ni2) -
                    AT_range_straggling_convolution(z_cm, range_cm, sigma_cm, ni2));

    denominator -= regul_factor_cm * pow(regul_factor_cm / 2.0, 1.0 / p) * exp(-(xi + zeta) * (xi + zeta) / 8.0);

    if (denominator < 0)
        denominator = 0.0;

    LET_d_keV_um *= nominator;
    LET_d_keV_um /= denominator;

    LET_d_keV_um *= 0.1; // from MeV/cm to keV/um

    return LET_d_keV_um;

}


void AT_LET_d_Wilkens_keV_um_multi(const long n,
                                   const double z_cm[],
                                   const double E_MeV,
                                   const double sigma_E_MeV,
                                   const long material_no,
                                   double LET_keV_um[]) {
    long i;
    for (i = 0; i < n; i++) {
        LET_keV_um[i] = AT_LET_d_Wilkens_keV_um_single(z_cm[i], E_MeV, sigma_E_MeV, material_no);
    }
}

/****************************************** RBE models *******************************************************/

double AT_proton_RBE_single(const double z_cm,
                            const double entrance_dose_Gy,
                            const double E_MeV,
                            const double sigma_E_MeV,
                            const double eps,
                            const double ref_alpha_beta_ratio,
                            const int rbe_model_no) {

    double rbe = 0.0;
    double let_keV_um = 0.0;
    double dose_Gy = 0.0;
    double fluence_cm2 = 1.0;

    double alpha_proton_to_alpha_ref = 0.0;
    double sqrt_beta_prot_to_beta_ref = 0.0;

    if (rbe_model_no == RBE_One) {
        rbe = 1.0;
    } else if (rbe_model_no == RBE_OnePointOne) {
        rbe = 1.1;
    } else {
        fluence_cm2 = entrance_dose_Gy;
        fluence_cm2 /= AT_dose_Bortfeld_Gy_single(z_cm, 1.0, E_MeV, sigma_E_MeV, Water_Liquid, eps);

        dose_Gy = AT_dose_Bortfeld_Gy_single(z_cm, fluence_cm2, E_MeV, sigma_E_MeV, Water_Liquid, eps);
        let_keV_um = AT_LET_d_Wilkens_keV_um_single(z_cm, E_MeV, sigma_E_MeV, Water_Liquid);

        switch (rbe_model_no) {
            case RBE_Carabe :
                alpha_proton_to_alpha_ref = 0.843 + 0.154 * 2.686 * let_keV_um / ref_alpha_beta_ratio;
                sqrt_beta_prot_to_beta_ref = 1.090 + 0.006 * 2.686 * let_keV_um / ref_alpha_beta_ratio;
                break;
            case RBE_Wedenberg :
                alpha_proton_to_alpha_ref = 1.000 + 0.434 * let_keV_um / ref_alpha_beta_ratio;
                sqrt_beta_prot_to_beta_ref = 1.000;
                break;
            case RBE_McNamara :
                alpha_proton_to_alpha_ref = 0.99064 + 0.35605 * let_keV_um / ref_alpha_beta_ratio;
                sqrt_beta_prot_to_beta_ref = 1.1012 - 0.0038703 * sqrt(ref_alpha_beta_ratio) * let_keV_um;
                break;
        }

        rbe = 1.0 / (2.0 * dose_Gy);
        rbe *= (sqrt(ref_alpha_beta_ratio * ref_alpha_beta_ratio +
                     4.0 * ref_alpha_beta_ratio * alpha_proton_to_alpha_ref * dose_Gy +
                     4.0 * sqrt_beta_prot_to_beta_ref * sqrt_beta_prot_to_beta_ref * dose_Gy * dose_Gy) -
                ref_alpha_beta_ratio);
    }

    return rbe;
}


void AT_proton_RBE_multi(const long n,
                         const double z_cm[],
                         const double entrance_dose_Gy,
                         const double E_MeV,
                         const double sigma_E_MeV,
                         const double eps,
                         const double ref_alpha_beta_ratio,
                         const int rbe_model_no,
                         double rbe[]) {
    long i;
    for (i = 0; i < n; i++) {
        rbe[i] = AT_proton_RBE_single(z_cm[i], entrance_dose_Gy, E_MeV, sigma_E_MeV, eps, ref_alpha_beta_ratio,
                                      rbe_model_no);
    }
}
