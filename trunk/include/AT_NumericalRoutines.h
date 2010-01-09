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

long int lminl(long int x, long int y);

long int lmaxl(long int x, long int y);

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

// finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
// is reported per element
// a vector "matches" of length n_elements has to be provided
void pmatchi(long* elements, long* n_elements, long* set, long* n_set, long* matches);

// finds character elements in a set (n elements) and returns indices - only one (the first) match
// is reported per element
// a vector "matches" of length n_elements has to be provided
void pmatchc(char** elements, long* n_elements, char** set, long* n_set, long* matches);

// finds a character element in a set and returns boolean match vector
// a vector "matches" of length n_set has to be provided
void matchc(char* element, char** set, long* n_set, bool* matches);

// finds a integer element in a set and returns boolean match vector
// a vector "matches" of length n_set has to be provided
void matchi(long* element, long* set, long* n_set, bool* matches);

////////////////////////////////////////////////////////////////////////////////////////////////////
// interpolation on a table: code (w/ adapted indices) from Numerical Recipes, 2rd ed., chapter 3.1
// added wrapping function interp which allows to chose degree of interpolation polynomial
// (1 = linear, 2 = quadratic, etc.)
void locate(float* xx, long* n, float* x, long* j);
void polint(float* xa, float* ya, long* n, float* x, float *y, float *dy);
void interp(float* xa, float* ya, long* n, long* n_pol, float* x, float *y, float *dy);


#endif /* AT_NUMERICALROUTINES_H_ */
