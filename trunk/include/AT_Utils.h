#ifndef AT_UTILS_H_
#define AT_UTILS_H_

/**
 *    AT_Utils.h
 *    ==========
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
#include "AT_Data.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

static float maxarg1;
static float maxarg2;
#define FMAX(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static float minarg1;
static float minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))

static long lmaxarg1;
static long lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2))

static long lminarg1;
static long lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ? (lminarg1) : (lminarg2))

int indent_counter;
char isp[30];
FILE * debf;

void indnt_inc();
void indnt_dec();
void indnt_init();


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

// get LET-data for given material
void getPSTARvalue(long* n, float* x, long* material_no, float* x_table, float* y_table, float* y);


#endif // AT_UTILS_H_
