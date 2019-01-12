/**
 * @brief A simple file to check the library and enable debugging
 */

/*
 *    AT_test.c
 *    ===================
 *
 *    Created on: 2009-06-08
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
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

#include "config.h"

#include "AT_StoppingPower.h"
#include "AT_CernlibFuns.h"
#include "AT_DataRange.h"
#include "AT_NumericalRoutines.h"
#include "AT_ProbabilityDistributions.h"

int test_cernlib() {
  double lambda, kappa, beta2;
  double cernlib_density, gsl_density, root_density;

  // compare two methods of calculating Landau PDF over wide range of lambda parameter :
  // fortran-based from CERNLIB and C-based from GSL
  for (lambda = -30.5; lambda < 1000.0; lambda += 0.1) {
    cernlib_density = CL_denlan(lambda);
    gsl_density = gsl_ran_landau_pdf(lambda);

    // for large numbers compare ratio and difference
    if (cernlib_density > 1e-6 && gsl_density > 1e-6) {
      if (fabs(cernlib_density - gsl_density) > 1e-6) {
        printf("%g %g %g (difference CERNLIB vs GSL too big)\n", lambda, cernlib_density, gsl_density);
        return EXIT_FAILURE;
      }
      if (fabs(cernlib_density / gsl_density - 1.0) > 2e-6) {
        printf("%g %g %g (ratio %g CERNLIB vs GSL too far from 1.0)\n",
               lambda,
               cernlib_density,
               gsl_density,
               cernlib_density / gsl_density - 1.0);
        return EXIT_FAILURE;
      }
    } else { // for small numbers compare only difference
      if (fabs(cernlib_density - gsl_density) > 1e-6) {
        printf("%g %g %g (difference CERNLIB vs GSL too big)\n", lambda, cernlib_density, gsl_density);
        return EXIT_FAILURE;
      }
    }
  }

  // compare two methods of calculating Vavilov PDF over some combination of kappa, beta2 and lambda parameters
  // fortran-based from CERNLIB and C-based from ROOT

  for (kappa = 0.01; kappa < 0.5; kappa += 0.01) {
    for (beta2 = 0.1; beta2 < 0.9; beta2 += 0.1) {
      for (lambda = -10.0; lambda < 8.0; lambda += 0.01) {

        CL_vavset(kappa, beta2);
        cernlib_density = CL_vavden(lambda);

        ROOT_GXXXC1 init = ROOT_vavset(kappa, beta2);
        root_density = ROOT_vav_pdf(lambda, &init);

        // for large numbers compare ratio and difference
        if (cernlib_density > 1e-4 && root_density > 1e-4) {
          if (fabs(cernlib_density - root_density) > 1e-4) {
            printf("kappa=%g beta2=%g lambda=%g %g %g (difference CERNLIB vs ROOT too big)\n",
                   kappa,
                   beta2,
                   lambda,
                   cernlib_density,
                   root_density);
            return EXIT_FAILURE;
          }
          if (fabs(cernlib_density / root_density - 1.0) > 1e-3) {
            printf("kappa=%g beta2=%g lambda=%g %g %g (ratio %g CERNLIB vs ROOT too far from 1.0)\n",
                   kappa,
                   beta2,
                   lambda,
                   cernlib_density,
                   root_density,
                   cernlib_density / root_density - 1.0);

            return EXIT_FAILURE;
          }
        } else { // for small numbers compare only difference
          if (fabs(cernlib_density - root_density) > 5e-3) {
            printf("kappa=%g beta2=%g lambda=%g  %g %g (difference CERNLIB vs ROOT too big)\n",
                   kappa,
                   beta2,
                   lambda,
                   cernlib_density,
                   root_density);
            return EXIT_FAILURE;
          }
        }

      }
    }
  }

  return EXIT_SUCCESS;
}

int main() {

  double Test = AT_CSDA_range_g_cm2_single(270.55,
                                           0.1,
                                           6012,
                                           1);

  printf("Range of 270.55 C-12 in water: %3.2f cm\n\n", Test);

  const double E_MeV_u[3] = {1, 10, 100};
  const long particle_no[3] = {6012, 6012, 6012};
  const long material_no = 2;

  double test[3] = {0., 0., 0.};

  AT_Mass_Stopping_Power(
      "/home/greilich/R/x86_64-pc-linux-gnu-library/3.2/libamtrack/extdata/FLUKA_DEDX_WATER_76.8eV.txt",
      3,
      E_MeV_u,
      particle_no,
      2,
      test);

  printf("Ergebnis FromFile: %e, %e, %e\n", test[0], test[1], test[2]);

  AT_Mass_Stopping_Power("Bethe",
                         3,
                         E_MeV_u,
                         particle_no,
                         material_no,
                         test);
  printf("Ergebnis Bethe: %e, %e, %e\n", test[0], test[1], test[2]);

  AT_Mass_Stopping_Power("PSTAR",
                         3,
                         E_MeV_u,
                         particle_no,
                         material_no,
                         test);
  printf("Ergebnis PSTAR: %e, %e, %e\n", test[0], test[1], test[2]);

  AT_Mass_Stopping_Power("ICRU",
                         3,
                         E_MeV_u,
                         particle_no,
                         material_no,
                         test);
  printf("Ergebnis ICRU: %e, %e, %e\n", test[0], test[1], test[2]);

  double l = 0.0, kappa = 1.0, beta = 0.5;
  CL_vavset(kappa, beta * beta);

  printf("Ergebnis Landau / Vavilov: %e, %e\n", CL_denlan(l), CL_vavden(l));

  if (test_cernlib() == EXIT_FAILURE)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
};

