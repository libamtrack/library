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
#include <stdbool.h>
#include <float.h>
#include <malloc.h>
#include <string.h>

extern int indent_counter;
extern char isp[];
extern FILE * debf;

long int lminl(const long int x, const long int y);

long int lmaxl(const long int x, const long int y);

void trim_long_to_int(long * n);

void AT_Funs(  const float*  fz,
    const float*  fR0,
    const float*  fsigma,
    const float* fni,
    float* funs);

void AT_fDyx(  const float*  fy,
    const float* fx,
    float* fDyx);

double  d_sign( const double *a, const double *b);
int pbdv_(  double *v, double *x, double *dv, double *dp, double *pdf, double *pdd);
int dvsa_(  double *va, double *x, double *pd);
int dvla_(  double *va, double *x, double *pd);
int vvla_(  double *va, double *x, double *pv);
int gamma_( const double *x, double *ga);

float gammln(const float xx);
void gcf(float *gammcf, const float a, const float x, float *gln);
void gser(float *gamser, const float a, float const x, float *gln);
float gammp(const float a, const float x);
float erff(const float x);
void nrerror(const char error_text[]);

/*   From Numerical Recipes in C, 2nd ed., 1992:
  Using Ridders' method, return the root of a function func known to lie between x1 and x2.
  The root, returned as zriddr, will be refined to an approximate accuracy xacc.
 */
float zriddr(float (*func)(float,void*), void * params, const float x1, const float x2, const float xacc);

// finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
// is reported per element
// a vector "matches" of length n_elements has to be provided
void pmatchi(const long* elements, const long* n_elements, const long* set, const long* n_set, long* matches);

// finds character elements in a set (n elements) and returns indices - only one (the first) match
// is reported per element
// a vector "matches" of length n_elements has to be provided
void pmatchc(const char** elements, const long* n_elements, const char* const* set, const long* n_set, long* matches);

// finds a character element in a set and returns boolean match vector
// a vector "matches" of length n_set has to be provided
void matchc(const char* element, const char* const * set, const long* n_set, bool* matches);

// finds a integer element in a set and returns boolean match vector
// a vector "matches" of length n_set has to be provided
void matchi(const long* element, const long* set, const long* n_set, bool* matches);

////////////////////////////////////////////////////////////////////////////////////////////////////
// interpolation on a table: code (w/ adapted indices) from Numerical Recipes, 2rd ed., chapter 3.1
// added wrapping function interp which allows to chose degree of interpolation polynomial
// (1 = linear, 2 = quadratic, etc.)
void locate(const float* xx, const long* n, const float* x, long* j);
void polint(const float* xa, const float* ya, const long* n, const float* x, float *y, float *dy);
void interp(const float* xa, const float* ya, const long* n, const long* n_pol, const float* x, float *y, float *dy);


#endif /* AT_NUMERICALROUTINES_H_ */
