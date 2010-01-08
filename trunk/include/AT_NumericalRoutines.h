#ifndef AT_NUMERICALROUTINES_H_
#define AT_NUMERICALROUTINES_H_

/**
 *    AT_GammaResponse.h
 *    ==================
 *
 *    Copyright 2006, 2009 Steffen Greilich / the libamtrack team
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


#include "AT_Constants.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern int indent_counter;
extern char isp[];
extern FILE * debf;

void AT_Funs(  float*  fz,
        float*  fR0,
        float*  fsigma,
        float*  fni,
        float*  funs);

void AT_fDyx(  float*  fy,
        float*  fx,
        float*  fDyx);

double   d_sign(  double *a, double *b);
int   pbdv_(  double *v, double *x, double *dv, double *dp, double *pdf, double *pdd);
int   dvsa_(  double *va, double *x, double *pd);
int   dvla_(  double *va, double *x, double *pd);
int   vvla_(  double *va, double *x, double *pv);
int   gamma_(  double *x, double *ga);

float gammln(float xx);
void gcf(float *gammcf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);
float gammp(float a, float x);
float erff(float x);
void nrerror(char error_text[]);

/*   From Numerical Recipes in C, 2nd ed., 1992:
  Using Ridders' method, return the root of a function func known to lie between x1 and x2.
  The root, returned as zriddr, will be refined to an approximate accuracy xacc.
 */
float zriddr(float (*func)(float,void*), void * params, float x1, float x2, float xacc);


#endif /* AT_NUMERICALROUTINES_H_ */
