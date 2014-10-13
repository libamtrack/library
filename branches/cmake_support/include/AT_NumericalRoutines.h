#ifndef AT_NUMERICALROUTINES_H_
#define AT_NUMERICALROUTINES_H_

/**
 * @brief Numerical Routines
 */


/*
 *    AT_NumericalRoutines.h
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


#include "AT_Constants.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
	#include <sys/malloc.h>
#else
	#include <malloc.h>
#endif

#include "gsl/gsl_math.h"


/**
 * Computes the convolution of a term (R0 - z)^(ni - 1) with a Gaussian
 * in z with variance sigma^2, i.e.
 * F(z, R0) = 1/(2*pi*sigma) * int_{-inf}^{R0}[ (R0 - z)^(ni - 1) * exp(-(z - z')^2/(2*sigma^2)) * dz']
 * that can be solved using the gamma function and the parabolic cylinder function:
 * F(z, R0) = 1/(2*pi*sigma) * exp((R0 - z)/(4*sigma^2)) * sigma^ni * gamma(ni) * D[-ni](-(R0-z)/sigma)
 * where D[-ni] is the parabolic cylinder function of order -ni
 *
 * The procedure is elucidated in Bortfeld, 1997, An analytical approximation of the Bragg curve for therapeutic
 * proton beams, Med. Phys. 24(12), 2024ff., Appendix A, Eqs. A1, A6
 *
 * This function uses gamma_ and AT_Dyx.
 *
 * Cave: Be careful to give the correct ni (not ni - 1)!
 *
 * @param[in]   z
 * @param[in]   R0
 * @param[in]   sigma
 * @param[in]   ni
 * @return      funs
 */
double AT_range_straggling_convolution(  const double z,
    const double R0,
    const double sigma,
    const double ni);

/*
 * parabolic cylinder functions can be implemented:
 *  - using Hermite polynomials (no implemented in GSL ?) for order ni, which is integer
 *  - using confluent hypergeometric function (gsl_sf_hyperg_U) for order ni, which is real number
 * more information here:
 * http://mathworld.wolfram.com/ParabolicCylinderFunction.html
 * http://www.gnu.org/software/gsl/manual/html_node/Hypergeometric-Functions.html
 */

/**
 *  Computes parabolic cylinder function Dy(x)
 *  using subroutine pbdv
 *  Original FORTRAN code mpbdv.f by Jianming Jin, Department of Electrical and
 *  Computer Engineering, University of Illinois at Urbana-Champaign
 *  http://jin.ece.uiuc.edu/routines/mpbdv.for
 *  Converted to C using f2c (version 20060506) by
 *  S. Greilich, reworked as subroutine for libamtrack.dll, abandoning
 *  f2c.h and libf2c.lib, as well as computation (returning) of derivatives
 *
 *  @param[in]   x       argument of Dy(x)
 *  @param[in]   y       order of Dy(x)
 *  return  Dyx
 */
double AT_Dyx(  double  y,
		double  x);


/**
 *  Computes parabolic cylinder function Dv(x) and its derivatives
 *  see comments for AT_Dyx
 *  The function calls dvsa for small |x| and dvla for large |x|
 *
 *  @param[in]   v       order of Dv(x)
 *  @param[in]   x       argument of Dv(x)
 *  @param[out]  dv      DV(na) = Dn+v0(x) with na = |n|, v0 = v-n, |v0| < 1, n = 0, +/-1, +/-2, ...
 *  @param[out]  dp      DP(na) = Dn+v0'(x) with na = |n|, v0 = v-n, |v0| < 1, n = 0, +/-1, +/-2, ...
 *  @param[out]  pdf     Dv(x)
 *  @param[out]  pdd     Dv'(x)
 */
int pbdv_(  double* v,
		double* x,
		double* dv,
		double* dp,
		double* pdf,
		double* pdd);


/**
 * Compute parabolic cylinder function Dv(x) for small argument
 * routines called: GAMMA
 * @param[in]  va order
 * @param[in]  x argument
 * @param[out] pd output Dv(x)
 */
int dvsa_(  double* va,
		double* x,
		double* pd);


/**
 * Compute parabolic cylinder function Dv(x) for large argument
 * Routines called:
 *             (1) VVLA for computing Vv(x) for large |x|
 *             (2) GAMMA for computing Gamma(x)
 * @param[in]  va order
 * @param[in]  x argument
 * @param[out] pd output Dv(x)
 */
int dvla_(  double* va,
		double* x,
		double* pd);


/**
 * Compute parabolic cylinder function Vv(x) for large argument
 * Routines called:
 *             (1) DVLA for computing Dv(x) for large |x|
 *             (2) GAMMA for computing Gamma(x)
 * @param[in] va order
 * @param[in] x argument
 * @param[out] pv output Vv(x)
 */
int vvla_(  double* va,
		double* x,
		double* pv);


/**
 * TODO
 * @param[in] a
 * @param[in] b
 * @return
 */
 double d_sign( const double a,
		 const double b);


/**
 * Compute parabolic gamma function
 * @param[in]  x argument (x is not equal to 0,-1,-2,...)
 * @param[out] ga output
 */
int AT_gamma_( const double* x,
		double* ga);


/* TODO it is logarithm of which gamma function ? it is not used */
/**
 * Numerical Recipes: Logarithm of gamma function
 * @param[in] xx argument for gamma function
 * @return
 */
double gammln( const double xx );


/**
 * Numerical Recipes: standard error handler
 * @param[in] error_text
 */
void nrerror( const char error_text[] );


// TODO equation solvers implemented in GSL should be tested and should possibly replace numerical recipes algorithms
/**
 * From Numerical Recipes in C, 2nd ed., 1992:
 * Using Ridders' method, return the root of a function func known to lie between x1 and x2.
 * The root, returned as zriddr, will be refined to an approximate accuracy xacc.
 * @param[in] func
 * @param[in] params
 * @param[in] x1
 * @param[in] x2
 * @param[in] xacc
 * @return
 */
double zriddr(double (*func)(double,void*),
		void * params,
		const double x1,
		const double x2,
		const double xacc);


/**
 * finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
 * is reported per element a vector "matches" of length n_elements has to be provided
 * @param[in] elements
 * @param[in] n_elements
 * @param[in] set
 * @param[in] n_set
 * @param[out] matches
 */
void are_elements_int(const int elements[],
		const int n_elements,
		const int set[],
		const int n_set,
		int matches[]);


/**
 * finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
 * is reported per element a vector "matches" of length n_elements has to be provided
 * @param[in] elements
 * @param[in] n_elements
 * @param[in] set
 * @param[in] n_set
 * @param[out] matches
 */
void find_elements_int(const long elements[],
		const long n_elements,
		const long set[],
		const long n_set,
		long matches[]);


/**
 * finds character elements in a set (n elements) and returns indices - only one (the first) match
 * is reported per element a vector "matches" of length n_elements has to be provided
 * @param[in] elements
 * @param[in] n_elements
 * @param[in] set
 * @param[in] n_set
 * @param[out] matches (array of size n_elements)
 */
void find_elements_char(const char** elements,
		const long n_elements,
		const char* const* set,
		const long n_set,
		long matches[]);


/**
 * finds a character element in a set and returns boolean match vector
 * a vector "matches" of length n_set has to be provided
 * @param[in] element
 * @param[in] set
 * @param[in] n_set
 * @param[out] matches
 */
void is_element_char(const char element[],
		const char* const * set,
		const long n_set,
		bool matches[]);


/**
 * finds a integer element in a set and returns boolean match vector
 * a vector "matches" of length n_set has to be provided
 * @param[in] element
 * @param[in] set
 * @param[in] n_set
 * @param[out] matches
 */
void is_element_int(const long element,
		const long set[],
		const long n_set,
		bool matches[]);


/**
 * Sums all elements in a data array
 * @param[in]  n                number of array elements
 * @param[out]  data            data to sum (array of size n)
 * @return  sum of all elements
 */
 double AT_sum(     const long n,
              const double data[]);


/**
 * Normalizes a data array by dividing each element by the sum of all elements
 * @param[in]  n                number of array elements
 * @param[in]  data             data to normalize
 * @param[out] normalized_data  results
 */
 void AT_normalize(     const long n,
                    const double data[],
                    double normalized_data[]);


/**
 * interpolation on a table: code (w/ adapted indices) from Numerical Recipes, 2rd ed., chapter 3.1
 * added wrapping function interp which allows to chose degree of interpolation polynomial
 * (1 = linear, 2 = quadratic, etc.)
 * @param[in] xx  TODO
 * @param[in] n   TODO
 * @param[in] x   TODO
 * return j   result
 */
long locate( const double xx[],
		const long n,
		const double x);


/**
 * Locates such index i that following inequation holds:
 * xx[i-1] <= x < x[i]
 * @param[in] xx
 * @param[in] lowest_index
 * @param[in] highest_index
 * @param[in] x
 * @param[in] index_in_row
 * @return
 */
long locate_index_in_2d_table( const double xx[][2],
		const long lowest_index,
		const long highest_index,
		const double x,
		int index_in_row);


/**
 * TODO
 * @param[in] input_data_x (array of size length_of_input_data)
 * @param[in] input_data_y (array of size length_of_input_data)
 * @param[in] length_of_input_data
 * @param[in] intermediate_x
 * @return
 */
double AT_get_interpolated_y_from_input_table( const double input_data_x[],
		const double input_data_y[],
		const long   length_of_input_data,
		const double intermediate_x);


/**
 * TODO
 * @param[in] input_data_xy
 * @param[in] length_of_input_data
 * @param[in] intermediate_x
 * @return
 */
double AT_get_interpolated_y_from_input_2d_table( const double input_data_xy[][2],
		const long length_of_input_data,
		const double intermediate_x);


/**
 *TODO
 * @param[in] input_data_xy
 * @param[in] lowest_index
 * @param[in] highest_index
 * @param[in] intermediate_y
 * @return
 */
double AT_get_interpolated_x_from_input_2d_table( const double input_data_xy[][2],
		const long lowest_index,
		const long highest_index,
		const double intermediate_y);


/**
 * TODO
 * @param[in] left_x
 * @param[in] left_y
 * @param[in] right_x
 * @param[in] right_y
 * @param[in] intermediate_x
 * @return
 */
double AT_get_interpolated_y_from_interval( const double left_x,
		const double left_y,
		const double right_x,
		const double right_y,
		const double intermediate_x);


// TODO implement linear interpolation on logarithmic scale


#endif /* AT_NUMERICALROUTINES_H_ */
