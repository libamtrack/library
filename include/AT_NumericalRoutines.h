#ifndef AT_NUMERICALROUTINES_H_
#define AT_NUMERICALROUTINES_H_

/**
 * @file
 * @brief Numerical Routines
 */


/*
 *    AT_NumericalRoutines.h
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


/**
 * TODO
 * @param x
 * @param y
 * @return
 */
long int lminl(const long int x, const long int y);


/**
 * TODO
 * @param x
 * @param y
 * @return
 */
long int lmaxl(const long int x, const long int y);


/**
 * TODO
 * @param fz
 * @param fR0
 * @param fsigma
 * @param fni
 * @param funs
 */
void AT_Funs(  const float*  fz,
    const float*  fR0,
    const float*  fsigma,
    const float* fni,
    float* funs);


/**
 * TODO
 * @param fy
 * @param fx
 * @param fDyx
 */
void AT_fDyx(  const float*  fy,
    const float* fx,
    float* fDyx);


/**
 * TODO
 * @param a
 * @param b
 * @return
 */
double  d_sign( const double *a, const double *b);


/**
 * TODO
 * @param v
 * @param x
 * @param dv
 * @param dp
 * @param pdf
 * @param pdd
 * @return
 */
int pbdv_(  double *v, double *x, double *dv, double *dp, double *pdf, double *pdd);


/**
 * Compute parabolic cylinder function Dv(x) for small argument
 * routines called: GAMMA
 * @param x argument
 * @param va order
 * @param pd output Dv(x)
 */
int dvsa_(  double *va, double *x, double *pd);


/**
 * Compute parabolic cylinder function Dv(x) for large argument
 * Routines called:
 *             (1) VVLA for computing Vv(x) for large |x|
 *             (2) GAMMA for computing �(x)
 * @param x argument
 * @param va order
 * @param pd output Dv(x)
 */
int dvla_(  double *va, double *x, double *pd);


/**
 * Compute parabolic cylinder function Vv(x) for large argument
 * Routines called:
 *             (1) DVLA for computing Dv(x) for large |x|
 *             (2) GAMMA for computing �(x)
 * @param x argument
 * @param va order
 * @param pv output Vv(x)
 */
int vvla_(  double *va, double *x, double *pv);


/**
 * Compute parabolic gamma function
 * @param x argument (x is not equal to 0,-1,-2,...)
 * @param ga output
 */
int gamma_( const double *x, double *ga);


/**
 * TODO
 * @param xx
 * @return
 */
float gammln(const float xx);


/**
 * TODO
 * @param gammcf
 * @param a
 * @param x
 * @param gln
 */
void gcf(float *gammcf, const float a, const float x, float *gln);


/**
 * TODO
 * @param gamser
 * @param a
 * @param x
 * @param gln
 */
void gser(float *gamser, const float a, float const x, float *gln);


/**
 * TODO
 * @param a
 * @param x
 * @return
 */
float gammp(const float a, const float x);


/**
 * TODO
 * @param x
 * @return
 */
float erff(const float x);


/**
 * Numerical Recipes standard error handler
 * @param error_text
 */
void nrerror(const char error_text[]);


/**
 * From Numerical Recipes in C, 2nd ed., 1992:
 * Using Ridders' method, return the root of a function func known to lie between x1 and x2.
 * The root, returned as zriddr, will be refined to an approximate accuracy xacc.
 * @param func
 * @param params
 * @param x1
 * @param x2
 * @param xacc
 * @return
 */
float zriddr(float (*func)(float,void*), void * params, const float x1, const float x2, const float xacc);


/**
 * finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
 * is reported per element a vector "matches" of length n_elements has to be provided
 * @param elements
 * @param n_elements
 * @param set
 * @param n_set
 * @param matches
 */
void are_elements_int(const int* elements, const int n_elements, const int* set, const int n_set, int* matches);


/**
 * finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
 * is reported per element a vector "matches" of length n_elements has to be provided
 * @param elements
 * @param n_elements
 * @param set
 * @param n_set
 * @param matches
 */
void find_elements_int(const long elements[], const long n_elements, const long set[], const long n_set, long matches[]);


/**
 * finds character elements in a set (n elements) and returns indices - only one (the first) match
 * is reported per element a vector "matches" of length n_elements has to be provided
 * @param elements
 * @param n_elements
 * @param set
 * @param n_set
 * @param matches
 */
void find_elements_char(const char** elements, const long* n_elements, const char* const* set, const long* n_set, long* matches);


/**
 * finds a character element in a set and returns boolean match vector
 * a vector "matches" of length n_set has to be provided
 * @param element
 * @param set
 * @param n_set
 * @param matches
 */
void is_element_char(const char* element, const char* const * set, const long* n_set, bool* matches);


/**
 * finds a integer element in a set and returns boolean match vector
 * a vector "matches" of length n_set has to be provided
 * @param element
 * @param set
 * @param n_set
 * @param matches
 */
void is_element_int(const long* element, const long* set, const long* n_set, bool* matches);


/**
 * interpolation on a table: code (w/ adapted indices) from Numerical Recipes, 2rd ed., chapter 3.1
 * added wrapping function interp which allows to chose degree of interpolation polynomial
 * (1 = linear, 2 = quadratic, etc.)
 * @param xx
 * @param n
 * @param x
 * @param j
 */
void locate(const float* xx, const long* n, const float* x, long* j);


/**
 * TODO
 * @param xa
 * @param ya
 * @param n
 * @param x
 * @param y
 * @param dy
 */
void polint(const float* xa, const float* ya, const long* n, const float* x, float *y, float *dy);


/**
 * TODO
 * @param xa
 * @param ya
 * @param n
 * @param n_pol
 * @param x
 * @param y
 * @param dy
 */
void interp(const float* xa, const float* ya, const long* n, const long* n_pol, const float* x, float *y, float *dy);


#endif /* AT_NUMERICALROUTINES_H_ */
