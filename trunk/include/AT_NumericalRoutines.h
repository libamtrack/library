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
#include <malloc.h>
#include <string.h>

#include "gsl/gsl_math.h"


/**
 * Numerical Recipes: Logarithm of gamma function
 * @param xx argument for gamma function
 * @return
 */
double gammln(const double xx);


/**
 * Numerical Recipes: Continued fraction used by gammp and gammq
 * @param gammcf
 * @param a
 * @param x
 * @param gln
 */
void gcf(double *gammcf, const double a, const double x, double *gln);


/**
 * Numerical Recipes: Series used by gammp and gammq
 * @param gamser
 * @param a
 * @param x
 * @param gln
 */
void gser(double *gamser, const double a, double const x, double *gln);


/**
 * Numerical Recipes: Incomplete gamma function
 * @param a
 * @param x
 * @return
 */
double gammp(const double a, const double x);


/**
 * Numerical Recipes: standard error handler
 * @param error_text
 */
void nrerror(const char error_text[]);


// TODO equation solvers implemented in GSL should be tested and should possibly replace numerical recipes algorithms
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
double zriddr(double (*func)(double,void*), void * params, const double x1, const double x2, const double xacc);


/**
 * finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
 * is reported per element a vector "matches" of length n_elements has to be provided
 * @param elements
 * @param n_elements
 * @param set
 * @param n_set
 * @param matches
 */
void are_elements_int(const int elements[], const int n_elements, const int set[], const int n_set, int matches[]);


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
void find_elements_char(const char** elements, const long n_elements, const char* const * set, const long n_set, long matches[]);


/**
 * finds a character element in a set and returns boolean match vector
 * a vector "matches" of length n_set has to be provided
 * @param element
 * @param set
 * @param n_set
 * @param matches
 */
void is_element_char(const char element[], const char* const * set, const long n_set, bool matches[]);


/**
 * finds a integer element in a set and returns boolean match vector
 * a vector "matches" of length n_set has to be provided
 * @param element
 * @param set
 * @param n_set
 * @param matches
 */
void is_element_int(const long element, const long set[], const long n_set, bool matches[]);


/**
 * interpolation on a table: code (w/ adapted indices) from Numerical Recipes, 2rd ed., chapter 3.1
 * added wrapping function interp which allows to chose degree of interpolation polynomial
 * (1 = linear, 2 = quadratic, etc.)
 * @param xx  TODO
 * @param n   TODO
 * @param x   TODO
 * return j   result
 */
long locate(const double xx[], const long n, const double x);


/**
 * Numerical Recipes: Polynomial interpolation
 * @param xa
 * @param ya
 * @param n
 * @param x
 * @param y
 * @param dy
 */
void polint(const double xa[], const double ya[], const long n, const double x, double *y, double *dy);


/**
 * TODO: ???
 * @param xa
 * @param ya
 * @param n
 * @param n_pol
 * @param x
 * @param y
 * @param dy
 */
void interp(const double xa[], const double ya[], const long n, const long n_pol, const double x, double y[], double dy[]);


#endif /* AT_NUMERICALROUTINES_H_ */
