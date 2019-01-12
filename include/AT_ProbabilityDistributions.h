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

/**
 * IDF of Landau distribution
 * Adapted from gsl_ran_landau by explicitly specifying X variable (from uniform distribution [0,1])
 * The code is based on CERNLIB G110 method (see http://hep.fi.infn.it/cernlib.pdf)
 * @param[in] X - random number from uniform distribution [0,1]
 * @return Landau random number
 */
double CL_ranlan_idf(const double X);

#endif //_AT_PROBABILITYDISTRIBUTIONS_H
