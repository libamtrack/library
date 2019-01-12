#ifndef _AT_PROBABILITYDISTRIBUTIONS_H
#define _AT_PROBABILITYDISTRIBUTIONS_H

/**
 * @brief Probability Distributions (Landau, Vavilov)
 */


/*
 *    AT_ProbabilityDistributions.h
 *    ==================
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

#include <math.h>
#include <stdio.h>

#include <gsl/gsl_randist.h>


/**
 * IDF of Landau distribution
 * Adapted from gsl_ran_landau by explicitly specifying X variable (from uniform distribution [0,1])
 * The code is based on CERNLIB G110 RANLAN method (see http://hep.fi.infn.it/cernlib.pdf)
 * @param[in] X - random number from uniform distribution [0,1]
 * @return Landau random number
 */
double CL_ranlan_idf(const double X);


/**
 * IDF of Landau distribution
 * code copied from http://git.savannah.gnu.org/cgit/gsl.git/tree/randist/landau.c
 * and compatible with https://github.com/root-project/root/blob/master/math/mathcore/src/ProbFuncMathCore.cxx
* The code is based on CERNLIB G110 DISLAN method (see http://hep.fi.infn.it/cernlib.pdf)
 * @param X
 * @return
 */
double CL_ranlan_cdf(const double X);

/**
 * Common block filled by VAVSET in CERBLIB implementation
 * or by VavilovFast::SetKappaBeta2 in ROOT
 *
 * @struct ROOT_GXXXC1
 */
typedef struct {
  double fAC[14];
  double fHC[9];
  double fWCM[201];
  int    fItype;
  int    fNpt;
} ROOT_GXXXC1;

/**
 * Initialisation of Vavilov distribution
 * The code is based on ROOT VavilovFast class
 * @param[in] kappa - The parameter \f$\kappa\f$, which should be in the range \f$0.01 \le \kappa \le 10 \f$
 * @param[in] beta2 - The parameter \f$\beta^2\f$, which must be in the range \f$0 \le \beta^2 \le 1 \f$
 * @return structure representing common block with precalculated data
 */
ROOT_GXXXC1 ROOT_vavset(const double kappa, const double beta2);

/**
 * PDF of Vavilov distribution
 * The code is based on ROOT VavilovFast class
 * @param[in] x - The Landau parameter \f$x = \lambda_L\f$
 * @param[in] init - The precalculated data
 * @return PDF value
 */
double ROOT_vav_pdf(const double x, const ROOT_GXXXC1 * init);

#endif //_AT_PROBABILITYDISTRIBUTIONS_H
