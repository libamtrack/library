#ifndef AT_CERNLIBFUNS_H_
#define AT_CERNLIBFUNS_H_

/**
 * @brief Algorithms from CERNLIB
 */

/*
 *    AT_CernlibFuns.h
 *    ================
 *
 *    Created on: 15.06.2015
 *    Creator: greilich
 *
 *    Copyright 2006, 2015 The libamtrack team
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

double CL_denlan( double lambda_landau );
double CL_ranlan( double rnd );

void   CL_vavset( double kappa, double beta2 );
double CL_vavden( double lambda_vavilov );
double CL_vavran( double kappa, double beta2, double rnd );

#endif /* AT_CERNLIBFUNS_H_ */
