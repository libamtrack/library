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

double AT_dose_Bortfeld_Gy( const double z_cm,
                            const double E_MeV_u,
                            const double fluence_cm2,
                            const double sigma_E_MeV_u,
                            const long material_no,
                            const double eps){

    // [1] Bortfeld, 1997, An analytical approximation of the Bragg curve for therapeutic proton beams, Med. Phys. 24(12), 2024ff

    const double beta_cm = 0.012; // slope parameter of fluence reduction relation [1/cm]
    const double gamma = 0.6; // TODO

    double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])
    double range_cm = alpha * pow(E_MeV_u, p);  // range in [cm]

    double ni1 = 1.0 / p;
    double ni2 = 1.0 + (1.0 / p);

    double sigma_mono_cm = 0.012 * pow(range_cm, 0.935);
    double sigma_cm = sqrt( sigma_mono_cm * sigma_mono_cm + pow((sigma_E_MeV_u * alpha * p * pow(E_MeV_u, p - 1)),2) );

    double dose_Gy = fluence_cm2; // returned value

    double factor1 = AT_range_straggling_convolution( z_cm, range_cm, sigma_cm, ni1);

    double factor2 = AT_range_straggling_convolution( z_cm, range_cm, sigma_cm, ni2);
    factor2 *= ((beta_cm / p) + (gamma * beta_cm) + (eps / range_cm));
    factor2 *= sigma_cm;

    double  tmp1;
    AT_gamma_(&ni1, &tmp1);
    factor2  *=  tmp1;

    double  tmp2;
    AT_gamma_(&ni2, &tmp2);
    factor2  /=  tmp2;

    dose_Gy /= AT_density_g_cm3_from_material_no(material_no);
    dose_Gy /= (p * pow(alpha, ni1));
    dose_Gy /= (1.0 + beta_cm * range_cm);

    dose_Gy *= (factor1 + factor2);

    dose_Gy *= 1.6021766e-10;  // MeV/g into Gy or J/kg

    return dose_Gy;
}
