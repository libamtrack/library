#ifndef SGP_UTILS_H_
#define SGP_UTILS_H_

#include <stdbool.h>

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

#ifdef _DEBUG
int indent_counter = 0;
char isp[30] = "\0                            ";
FILE * debf;

void indnt_inc();
void indnt_dec();
void indnt_init();

#endif


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

float gammln(float xx);
void gcf(float *gammcf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);
float gammp(float a, float x);
float erff(float x);
void nrerror(char error_text[]);

/* 	From Numerical Recipes in C, 2nd ed., 1992:
	Using Ridders' method, return the root of a function func known to lie between x1 and x2.
	The root, returned as zriddr, will be refined to an approximate accuracy xacc.
 */
float zriddr(float (*func)(float,void*), void * params, float x1, float x2, float xacc);

#endif // SGP_UTILS_H_
