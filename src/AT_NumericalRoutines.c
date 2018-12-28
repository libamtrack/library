/**
 * @brief Numerical routines
 */

/*
 *    AT_NumericalRoutines.c
 *    ==============
 *
 *    Created on: 8.01.2010
 *    Creator: kongruencja
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

#include "AT_NumericalRoutines.h"

double AT_range_straggling_convolution(  const double z,
    const double R0,
    const double sigma,
    const double ni)
{
  double  F     =  0;
  double  u     =  R0 - z;              // Residual range
  assert( sigma != 0.0 );
  double  zeta  =  u / sigma;           // parameter to divide convolution into domains
                                        // for zeta < -5, F(z,R0) becomes very small
                                        // for zeta > 10, F(z,R0) gets very close to (R0-z)^(ni-1)

  if(zeta > -5.0 && zeta < 10){
    F  =   1.0 / (sqrt(2.0 * M_PI) * sigma);
    F  *=  exp(-1.0 * u * u / (4.0 * sigma * sigma)) * pow(sigma, ni);

    double  tmp;
    AT_gamma_(&ni, &tmp);
    F  *=  tmp;

    tmp = AT_Dyx(-1.0 * ni, -1.0 * zeta);

    F  *=  tmp;
  }

  if(zeta >= 10.0){
    F      =  pow(u, ni - 1.0);
  }
  return F;
}

double AT_Dyx(  double  y,  double  x)
{
  double  Dyx;
  double dp[101], dv[101];
  double pdd;

  pbdv_(&y, &x, dv, dp, &Dyx, &pdd);

  return Dyx;
}

int pbdv_(double *v,
    double *x,
    double *dv,
    double *dp,
    double *pdf,
    double *pdd)
{
  /* System generated locals */
  int i__1;

  /* TODO why those variables are static ???? */
  /* Local variables */
  static double f;
  static int k, l, m;
  static double f0, f1, s0, v0, v1, v2;
  static int ja, na;
  static double ep, pd, xa;
  static int nk;
  static double vh;
  static int nv;
  static double pd0, pd1;

  /* TODO why do we need extern subroutines here ? */
  extern /* Subroutine */ int dvla_(double *, double *, double *
  ), dvsa_(double *, double *, double *);

  static double c_b31 = 1.;

  xa = fabs(*x);
  vh = *v;
  *v += d_sign(c_b31, *v);
  nv = (int) (*v);
  v0 = *v - nv;
  na = abs(nv);
  ep = exp(*x * -.25 * *x);
  if (na >= 1) {
    ja = 1;
  }
  if (*v >= 0.f) {
    if (v0 == 0.f) {
      pd0 = ep;
      pd1 = *x * ep;
    } else {
      i__1 = ja;
      for (l = 0; l <= i__1; ++l) {
        v1 = v0 + l;
        if (xa <= 5.8f) {
          dvsa_(&v1, x, &pd1);
        }
        if (xa > 5.8f) {
          dvla_(&v1, x, &pd1);
        }
        if (l == 0) {
          pd0 = pd1;
        }
        /* L10: */
      }
    }
    dv[0] = pd0;
    dv[1] = pd1;
    i__1 = na;
    for (k = 2; k <= i__1; ++k) {
      *pdf = *x * pd1 - (k + v0 - 1.) * pd0;
      dv[k] = *pdf;
      pd0 = pd1;
      /* L15: */
      pd1 = *pdf;
    }
  } else {
    if (*x <= 0.f) {
      if (xa <= 5.8) {
        dvsa_(&v0, x, &pd0);
        v1 = v0 - 1.;
        dvsa_(&v1, x, &pd1);
      } else {
        dvla_(&v0, x, &pd0);
        v1 = v0 - 1.;
        dvla_(&v1, x, &pd1);
      }
      dv[0] = pd0;
      dv[1] = pd1;
      i__1 = na;
      for (k = 2; k <= i__1; ++k) {
        pd = (-(*x) * pd1 + pd0) / (k - 1. - v0);
        dv[k] = pd;
        pd0 = pd1;
        /* L20: */
        pd1 = pd;
      }
    } else if (*x <= 2.f) {
      v2 = nv + v0;
      if (nv == 0) {
        v2 += -1.;
      }
      nk = (int) (-v2);
      dvsa_(&v2, x, &f1);
      v1 = v2 + 1.;
      dvsa_(&v1, x, &f0);
      dv[nk] = f1;
      dv[nk - 1] = f0;
      for (k = nk - 2; k >= 0; --k) {
        f = *x * f0 + (k - v0 + 1.) * f1;
        dv[k] = f;
        f1 = f0;
        /* L25: */
        f0 = f;
      }
    } else {
      if (xa <= 5.8f) {
        dvsa_(&v0, x, &pd0);
      }
      if (xa > 5.8f) {
        dvla_(&v0, x, &pd0);
      }
      dv[0] = pd0;
      m = na + 100;
      f1 = 0.;
      f0 = 1e-30;
      for (k = m; k >= 0; --k) {
        f = *x * f0 + (k - v0 + 1.) * f1;
        if (k <= na) {
          dv[k] = f;
        }
        f1 = f0;
        /* L30: */
        f0 = f;
      }
      s0 = pd0 / f;
      i__1 = na;
      for (k = 0; k <= i__1; ++k) {
        /* L35: */
        dv[k] = s0 * dv[k];
      }
    }
  }
  i__1 = na - 1;
  for (k = 0; k <= i__1; ++k) {
    v1 = fabs(v0) + k;
    if (*v >= 0.) {
      dp[k] = *x * .5 * dv[k] - dv[k + 1];
    } else {
      dp[k] = *x * -.5 * dv[k] - v1 * dv[k + 1];
    }
    /* L40: */
  }
  *pdf = dv[na - 1];
  *pdd = dp[na - 1];
  *v = vh;
  return 0;
} /* pbdv_ */


 double d_sign(const double a, const double b)
{
  double x = (a >= 0 ? a : - a);
  return( b >= 0 ? x : -x);
}


int dvsa_(double *va,
    double *x,
    double *pd)
{
  /* System generated locals */
  double d__1;

  /* TODO why those variables are static ???? */
  /* Local variables */
  static int m;
  static double r__, a0, g0, g1, r1, ep, gm, pi, vm, vt, ga0, va0, sq2,
  eps;
  static double c_b40 = 2.;

  eps = 1e-15;
  pi = M_PI;
  sq2 = sqrt(2.);
  ep = exp(*x * -.25 * *x);
  va0 = (1. - *va) * .5;
  if (*va == 0.f) {
    *pd = ep;
  } else {
    if (*x == 0.f) {
      if (va0 <= 0.f && va0 == (double) ((int) va0)) {
        *pd = 0.;
      } else {
        AT_gamma_(&va0, &ga0);
        d__1 = *va * -.5;
        *pd = sqrt(pi) / (pow(c_b40, d__1) * ga0);
      }
    } else {
      d__1 = -(*va);
      AT_gamma_(&d__1, &g1);
      d__1 = *va * -.5 - 1.;
      a0 = pow(c_b40, d__1) * ep / g1;
      vt = *va * -.5;
      AT_gamma_(&vt, &g0);
      *pd = g0;
      r__ = 1.;
      for (m = 1; m <= 250; ++m) {
        vm = (m - *va) * .5;
        AT_gamma_(&vm, &gm);
        r__ = -r__ * sq2 * *x / m;
        r1 = gm * r__;
        *pd += r1;
        if (fabs(r1) < fabs(*pd) * eps) {
          goto L15;
        }
        /* L10: */
      }
      L15:
      *pd = a0 * *pd;
    }
  }
  return 0;
} /* dvsa_ */


int dvla_(double *va, double *x, double *pd)
{
  /* System generated locals */
  double d__1;

  /* TODO why those variables are static ???? */
  /* Local variables */
  static int k;
  static double r__, a0, x1, gl, ep, pi, vl, eps;

  pi = M_PI;
  eps = 1e-12;
  ep = exp(*x * -.25f * *x);
  d__1 = fabs(*x);
  a0 = pow(d__1, *va) * ep;
  r__ = 1.;
  *pd = 1.;
  for (k = 1; k <= 16; ++k) {
    r__ = r__ * -.5 * (k * 2.f - *va - 1.f) * (k * 2.f - *va - 2.f) / (k *
        *x * *x);
    *pd += r__;
    if ((d__1 = r__ / *pd, fabs(d__1)) < eps) {
      goto L15;
    }
    /* L10: */
  }
  L15:
  *pd = a0 * *pd;
  if (*x < 0.) {
    x1 = -(*x);
    vvla_(va, &x1, &vl);
    d__1 = -(*va);
    AT_gamma_(&d__1, &gl);
    *pd = pi * vl / gl + cos(pi * *va) * *pd;
  }
  return 0;
} /* dvla_ */


int vvla_(double *va, double *x, double *pv)
{
  /* System generated locals */
  double d__1, d__2;

  /* TODO why those variables are static ???? */
  /* Local variables */
  static int k;
  static double r__, a0, x1, gl, qe, pi, pdl, dsl, eps;
  pi = M_PI;
  eps = 1e-12;
  qe = exp(*x * .25f * *x);
  d__1 = fabs(*x);
  d__2 = -(*va) - 1.;
  a0 = pow(d__1, d__2) * sqrt(2. / pi) * qe;
  r__ = 1.;
  *pv = 1.;
  for (k = 1; k <= 18; ++k) {
    r__ = r__ * .5 * (k * 2.f + *va - 1.f) * (k * 2.f + *va) / (k * *x * *
        x);
    *pv += r__;
    if ((d__1 = r__ / *pv, fabs(d__1)) < eps) {
      goto L15;
    }
    /* L10: */
  }
  L15:
  *pv = a0 * *pv;
  if (*x < 0.) {
    x1 = -(*x);
    dvla_(va, &x1, &pdl);
    d__1 = -(*va);
    AT_gamma_(&d__1, &gl);
    dsl = sin(pi * *va) * sin(pi * *va);
    *pv = dsl * gl / pi * pdl - cos(pi * *va) * *pv;
  }
  return 0;
} /* vvla_ */


int AT_gamma_(const double *x, double *ga)
{
  /* TODO why those variables are static ???? */
  /* Initialized data */
  static double g[26] = { 1.,.5772156649015329,-.6558780715202538,
      -.0420026350340952,.1665386113822915,-.0421977345555443,
      -.009621971527877,.007218943246663,-.0011651675918591,
      -2.152416741149e-4,1.280502823882e-4,-2.01348547807e-5,
      -1.2504934821e-6,1.133027232e-6,-2.056338417e-7,6.116095e-9,
      5.0020075e-9,-1.1812746e-9,1.043427e-10,7.7823e-12,-3.6968e-12,
      5.1e-13,-2.06e-14,-5.4e-15,1.4e-15,1e-16 };

  /* System generated locals */
  int i__1;

  /* Local variables */
  static int k, m;
  static double r__, z__;
  static int m1;
  static double pi, gr;

  pi = M_PI;
  if (*x == (double) ((int) (*x))) {
    if (*x > 0.) {
      *ga = 1.;
      m1 = (int) (*x - 1);
      i__1 = m1;
      for (k = 2; k <= i__1; ++k) {
        /* L10: */
        *ga *= k;
      }
    } else {
      *ga = 1e300;
    }
  } else {
    if (fabs(*x) > 1.) {
      z__ = fabs(*x);
      m = (int) z__;
      r__ = 1.;
      i__1 = m;
      for (k = 1; k <= i__1; ++k) {
        /* L15: */
        r__ *= z__ - k;
      }
      z__ -= m;
    } else {
      z__ = *x;
    }
    gr = g[25];
    for (k = 25; k >= 1; --k) {
      /* L20: */
      gr = gr * z__ + g[k - 1];
    }
    *ga = 1. / (gr * z__);
    if (fabs(*x) > 1.) {
      *ga *= r__;
      if (*x < 0.) {
        *ga = -pi / (*x * *ga * sin(pi * *x));
      }
    }
  }
  return 0;
} /* gamma_ */




double gammln(const double xx)
{
  /* TODO why those variables are static ???? */
  double x,y,tmp,ser;
  static double cof[6]=  {  76.18009172947146,-86.50532032941677,
      24.01409824083091,-1.231739572450155,
      0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return (-tmp+log(2.5066282746310005*ser/x));
}


void nrerror(const char error_text[])
{
#ifndef NDEBUG
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
#endif
 }


/* TODO could SIGN be replaced by some system call ?*/
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAXIT 60
#define UNUSED (-1.11e30)


double zriddr(double (*func)(double,void*), void * params, const double x1, const double x2, const double xacc)
{
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
  fl=(*func)(x1,params);
  fh=(*func)(x2,params);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;
    ans=UNUSED;                         //  Any highly unlikely value, to simplify logic below.
    for (j=1;j<=MAXIT;j++) {
      xm=0.5*(xl+xh);
      fm=(*func)(xm,params);                     // First of two function evaluations per iteration
      s=sqrt(fm*fm-fl*fh);
      if (s == 0.0) return ans;
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);   // Updating formula.
      if (fabs(xnew-ans) <= xacc) return ans;
      ans=xnew;
      fnew=(*func)(ans,params);                   // Second of two function evaluations per
      if (fnew == 0.0) return ans;             // iteration.
      if (SIGN(fm,fnew) != fm) {               // Bookkeeping to keep the root bracketed
        xl=xm;                       // on next iteration.
        fl=fm;
        xh=ans;
        fh=fnew;
      } else if (SIGN(fl,fnew) != fl) {
        xh=ans;
        fh=fnew;
      } else if (SIGN(fh,fnew) != fh) {
        xl=ans;
        fl=fnew;
      } else nrerror("never get here.");
      if (fabs(xh-xl) <= xacc) return ans;
    }
    nrerror("zriddr exceed maximum iterations");
  }
  else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    nrerror("root must be bracketed in zriddr");
  }
  return 0.0;                         // Never get here.
}


void are_elements_int(const int elements[], const int n_elements, const int set[], const int n_set, int matches[]){
  long  i;
  for (i = 0; i < n_elements; i++){
    matches[i] = 0;

    while ((set[matches[i]] != elements[i]) && (matches[i] < n_set)){
      matches[i]++;
    }

    if (matches[i] == n_set) {
      matches[i] = -1;
    }
  }
}


void find_elements_int(const long elements[], const long n_elements, const long set[], const long n_set, long matches[]){
  long  i;
  for (i = 0; i < n_elements; i++){
    matches[i] = 0;

    while ((set[matches[i]] != elements[i]) && (matches[i] < n_set)){
      matches[i]++;
    }

    if (matches[i] == n_set) {
      matches[i] = -1;
    }
  }
}


void find_elements_char(const char** elements, const long n_elements, const char* const * set, const long n_set, long matches[]){

  long  i;
  for (i = 0; i < n_elements; i++){
    matches[i] = 0;

    while ( strcmp( set[matches[i]], elements[i]) != 0 ){
      matches[i]++;
      if( matches[i] == n_set ){
    	  matches[i] = -1;
    	  break;
      }
    }
  }
}


void is_element_char(const char element[], const char* const * set, const long n_set, bool matches[]){

  long  i;
  for (i = 0; i < n_set; i++){
    if(strcmp(element, set[i])==0){
      matches[i]  = true;
    }else{
      matches[i]  = false;
    }
  }
}


int is_element_int(const long element, const long set[], const long n_set, bool matches[]){

  long  i, n = 0;
  for (i = 0; i < n_set; i++){
    if(element == set[i]){
      matches[i]  = true;
      n++;
    } else{
      matches[i]  = false;
    }
  }
  return n;
}


 double AT_sum(     const long n,
              const double data[])
{
  long i;
  double sum = 0.0;
  for (i = 0; i < n; i++){
    sum += data[i];
  }
  return sum;
}


 void AT_normalize(     const long n,
                    const double data[],
                    double normalized_data[])
{
  long i;
  double sum = AT_sum(n, data);

  for (i = 0; i < n; i++){
    normalized_data[i] = data[i] / sum;
  }
}


long locate(const double xx[], const long n, const double x)
{
  long  ju, jm, jl, j;
  int    ascnd;

  jl    =  0;
  ju    =  n + 1;
  ascnd  =  (xx[n-1] >= xx[0]);
  while (ju - jl > 1){
    jm    =  (ju + jl) >> 1;
    if ((x >= xx[jm-1]) == ascnd)
      jl  =  jm;
    else
      ju  =  jm;
  }
  if ( x == xx[0]) j = 1;
  else if (x == xx[n - 1]) j = n - 1;
  else j  =  jl;
  return j;
}


long locate_index_in_2d_table(const double xx[][2], const long lowest_index, const long highest_index, const double x, int index_in_row)
{
  assert( index_in_row >= 0 );
  assert( index_in_row <= 2 );
  assert( lowest_index >= 0 );
  assert( highest_index >= lowest_index );

  long  j_upper, jm, j_lower, j;
  int    ascnd;

  j_lower    =  lowest_index;
  j_upper    =  highest_index;
  ascnd  =  (xx[highest_index][index_in_row] >= xx[lowest_index][index_in_row]);

  while (j_upper - j_lower > 1){
    jm    =  (j_upper + j_lower) >> 1;
    if ((x >= xx[jm-1][index_in_row]) == ascnd)
      j_lower  =  jm;
    else
      j_upper  =  jm;
  }
  if ( x == xx[lowest_index][index_in_row]) j = lowest_index + 1;
  else if (x == xx[highest_index][index_in_row]) j = highest_index;
  else j  =  j_lower;
  return j;
}


double AT_get_interpolated_y_from_input_table(const double input_data_x[], const double input_data_y[], const long length_of_input_data, const double intermediate_x){
	int i = locate( input_data_x, length_of_input_data, intermediate_x );

	assert( i >= 0 );
	assert( i < length_of_input_data);

	return AT_get_interpolated_y_from_interval( input_data_x[i-1], input_data_y[i-1], input_data_x[i], input_data_y[i], intermediate_x);
}


double AT_get_interpolated_y_from_input_2d_table(const double input_data_xy[][2], const long length_of_input_data, const double intermediate_x){
	int i = locate_index_in_2d_table( input_data_xy, 0, length_of_input_data-1, intermediate_x, 0 );

	assert( i >= 0 );
	assert( i < length_of_input_data);

	return AT_get_interpolated_y_from_interval( input_data_xy[i-1][0], input_data_xy[i-1][1], input_data_xy[i][0], input_data_xy[i][1], intermediate_x);
}


double AT_get_interpolated_x_from_input_2d_table(const double input_data_xy[][2], const long lowest_index, const long highest_index, const double intermediate_y){
	long i = locate_index_in_2d_table( input_data_xy, lowest_index, highest_index, intermediate_y, 1 );

	assert( i >= lowest_index );
	assert( i <= highest_index );

	if( input_data_xy[i][1] > input_data_xy[i-1][1]  ){
		return AT_get_interpolated_y_from_interval( input_data_xy[i-1][1], input_data_xy[i-1][0], input_data_xy[i][1], input_data_xy[i][0], intermediate_y);
	} else {
		return AT_get_interpolated_y_from_interval( input_data_xy[i][1], input_data_xy[i][0], input_data_xy[i-1][1], input_data_xy[i-1][0], intermediate_y);
	}
	return -1;
}


double AT_get_interpolated_y_from_interval(const double left_x, const double left_y, const double right_x, const double right_y, const double intermediate_x){
	// (x - left_x) / (right_x - left_x ) = (y - left_y) / (right_y - left_y)

//	assert( right_x > left_x); // Does that make sense? Is it necessary?
	if(right_x > left_x){
		assert( intermediate_x >= left_x);
		assert( intermediate_x <= right_x);
	}else{
		assert( intermediate_x <= left_x);
		assert( intermediate_x >= right_x);
	}
	if(right_x == left_x){
		return left_y;
	}else{
		return  left_y + (right_y - left_y)*((intermediate_x - left_x) / (right_x - left_x));
	}
}

double CL_ranlan_idf(const double X)
  {
    static double F[983] =
        {
            0.0000000,   /* Add empty element [0] to account for difference
                      between C and Fortran convention for lower bound. */
            00.000000, 00.000000, 00.000000, 00.000000, 00.000000,
            -2.244733, -2.204365, -2.168163, -2.135219, -2.104898,
            -2.076740, -2.050397, -2.025605, -2.002150, -1.979866,
            -1.958612, -1.938275, -1.918760, -1.899984, -1.881879,
            -1.864385, -1.847451, -1.831030, -1.815083, -1.799574,
            -1.784473, -1.769751, -1.755383, -1.741346, -1.727620,
            -1.714187, -1.701029, -1.688130, -1.675477, -1.663057,
            -1.650858, -1.638868, -1.627078, -1.615477, -1.604058,
            -1.592811, -1.581729, -1.570806, -1.560034, -1.549407,
            -1.538919, -1.528565, -1.518339, -1.508237, -1.498254,
            -1.488386, -1.478628, -1.468976, -1.459428, -1.449979,
            -1.440626, -1.431365, -1.422195, -1.413111, -1.404112,
            -1.395194, -1.386356, -1.377594, -1.368906, -1.360291,
            -1.351746, -1.343269, -1.334859, -1.326512, -1.318229,
            -1.310006, -1.301843, -1.293737, -1.285688, -1.277693,
            -1.269752, -1.261863, -1.254024, -1.246235, -1.238494,
            -1.230800, -1.223153, -1.215550, -1.207990, -1.200474,
            -1.192999, -1.185566, -1.178172, -1.170817, -1.163500,
            -1.156220, -1.148977, -1.141770, -1.134598, -1.127459,
            -1.120354, -1.113282, -1.106242, -1.099233, -1.092255,
            -1.085306, -1.078388, -1.071498, -1.064636, -1.057802,
            -1.050996, -1.044215, -1.037461, -1.030733, -1.024029,
            -1.017350, -1.010695, -1.004064, -0.997456, -0.990871,
            -0.984308, -0.977767, -0.971247, -0.964749, -0.958271,
            -0.951813, -0.945375, -0.938957, -0.932558, -0.926178,
            -0.919816, -0.913472, -0.907146, -0.900838, -0.894547,
            -0.888272, -0.882014, -0.875773, -0.869547, -0.863337,
            -0.857142, -0.850963, -0.844798, -0.838648, -0.832512,
            -0.826390, -0.820282, -0.814187, -0.808106, -0.802038,
            -0.795982, -0.789940, -0.783909, -0.777891, -0.771884,
            -0.765889, -0.759906, -0.753934, -0.747973, -0.742023,
            -0.736084, -0.730155, -0.724237, -0.718328, -0.712429,
            -0.706541, -0.700661, -0.694791, -0.688931, -0.683079,
            -0.677236, -0.671402, -0.665576, -0.659759, -0.653950,
            -0.648149, -0.642356, -0.636570, -0.630793, -0.625022,
            -0.619259, -0.613503, -0.607754, -0.602012, -0.596276,
            -0.590548, -0.584825, -0.579109, -0.573399, -0.567695,
            -0.561997, -0.556305, -0.550618, -0.544937, -0.539262,
            -0.533592, -0.527926, -0.522266, -0.516611, -0.510961,
            -0.505315, -0.499674, -0.494037, -0.488405, -0.482777,
            -0.477153, -0.471533, -0.465917, -0.460305, -0.454697,
            -0.449092, -0.443491, -0.437893, -0.432299, -0.426707,
            -0.421119, -0.415534, -0.409951, -0.404372, -0.398795,
            -0.393221, -0.387649, -0.382080, -0.376513, -0.370949,
            -0.365387, -0.359826, -0.354268, -0.348712, -0.343157,
            -0.337604, -0.332053, -0.326503, -0.320955, -0.315408,
            -0.309863, -0.304318, -0.298775, -0.293233, -0.287692,
            -0.282152, -0.276613, -0.271074, -0.265536, -0.259999,
            -0.254462, -0.248926, -0.243389, -0.237854, -0.232318,
            -0.226783, -0.221247, -0.215712, -0.210176, -0.204641,
            -0.199105, -0.193568, -0.188032, -0.182495, -0.176957,
            -0.171419, -0.165880, -0.160341, -0.154800, -0.149259,
            -0.143717, -0.138173, -0.132629, -0.127083, -0.121537,
            -0.115989, -0.110439, -0.104889, -0.099336, -0.093782,
            -0.088227, -0.082670, -0.077111, -0.071550, -0.065987,
            -0.060423, -0.054856, -0.049288, -0.043717, -0.038144,
            -0.032569, -0.026991, -0.021411, -0.015828, -0.010243,
            -0.004656, 00.000934, 00.006527, 00.012123, 00.017722,
            00.023323, 00.028928, 00.034535, 00.040146, 00.045759,
            00.051376, 00.056997, 00.062620, 00.068247, 00.073877,
            00.079511, 00.085149, 00.090790, 00.096435, 00.102083,
            00.107736, 00.113392, 00.119052, 00.124716, 00.130385,
            00.136057, 00.141734, 00.147414, 00.153100, 00.158789,
            00.164483, 00.170181, 00.175884, 00.181592, 00.187304,
            00.193021, 00.198743, 00.204469, 00.210201, 00.215937,
            00.221678, 00.227425, 00.233177, 00.238933, 00.244696,
            00.250463, 00.256236, 00.262014, 00.267798, 00.273587,
            00.279382, 00.285183, 00.290989, 00.296801, 00.302619,
            00.308443, 00.314273, 00.320109, 00.325951, 00.331799,
            00.337654, 00.343515, 00.349382, 00.355255, 00.361135,
            00.367022, 00.372915, 00.378815, 00.384721, 00.390634,
            00.396554, 00.402481, 00.408415, 00.414356, 00.420304,
            00.426260, 00.432222, 00.438192, 00.444169, 00.450153,
            00.456145, 00.462144, 00.468151, 00.474166, 00.480188,
            00.486218, 00.492256, 00.498302, 00.504356, 00.510418,
            00.516488, 00.522566, 00.528653, 00.534747, 00.540850,
            00.546962, 00.553082, 00.559210, 00.565347, 00.571493,
            00.577648, 00.583811, 00.589983, 00.596164, 00.602355,
            00.608554, 00.614762, 00.620980, 00.627207, 00.633444,
            00.639689, 00.645945, 00.652210, 00.658484, 00.664768,
            00.671062, 00.677366, 00.683680, 00.690004, 00.696338,
            00.702682, 00.709036, 00.715400, 00.721775, 00.728160,
            00.734556, 00.740963, 00.747379, 00.753807, 00.760246,
            00.766695, 00.773155, 00.779627, 00.786109, 00.792603,
            00.799107, 00.805624, 00.812151, 00.818690, 00.825241,
            00.831803, 00.838377, 00.844962, 00.851560, 00.858170,
            00.864791, 00.871425, 00.878071, 00.884729, 00.891399,
            00.898082, 00.904778, 00.911486, 00.918206, 00.924940,
            00.931686, 00.938446, 00.945218, 00.952003, 00.958802,
            00.965614, 00.972439, 00.979278, 00.986130, 00.992996,
            00.999875, 01.006769, 01.013676, 01.020597, 01.027533,
            01.034482, 01.041446, 01.048424, 01.055417, 01.062424,
            01.069446, 01.076482, 01.083534, 01.090600, 01.097681,
            01.104778, 01.111889, 01.119016, 01.126159, 01.133316,
            01.140490, 01.147679, 01.154884, 01.162105, 01.169342,
            01.176595, 01.183864, 01.191149, 01.198451, 01.205770,
            01.213105, 01.220457, 01.227826, 01.235211, 01.242614,
            01.250034, 01.257471, 01.264926, 01.272398, 01.279888,
            01.287395, 01.294921, 01.302464, 01.310026, 01.317605,
            01.325203, 01.332819, 01.340454, 01.348108, 01.355780,
            01.363472, 01.371182, 01.378912, 01.386660, 01.394429,
            01.402216, 01.410024, 01.417851, 01.425698, 01.433565,
            01.441453, 01.449360, 01.457288, 01.465237, 01.473206,
            01.481196, 01.489208, 01.497240, 01.505293, 01.513368,
            01.521465, 01.529583, 01.537723, 01.545885, 01.554068,
            01.562275, 01.570503, 01.578754, 01.587028, 01.595325,
            01.603644, 01.611987, 01.620353, 01.628743, 01.637156,
            01.645593, 01.654053, 01.662538, 01.671047, 01.679581,
            01.688139, 01.696721, 01.705329, 01.713961, 01.722619,
            01.731303, 01.740011, 01.748746, 01.757506, 01.766293,
            01.775106, 01.783945, 01.792810, 01.801703, 01.810623,
            01.819569, 01.828543, 01.837545, 01.846574, 01.855631,
            01.864717, 01.873830, 01.882972, 01.892143, 01.901343,
            01.910572, 01.919830, 01.929117, 01.938434, 01.947781,
            01.957158, 01.966566, 01.976004, 01.985473, 01.994972,
            02.004503, 02.014065, 02.023659, 02.033285, 02.042943,
            02.052633, 02.062355, 02.072110, 02.081899, 02.091720,
            02.101575, 02.111464, 02.121386, 02.131343, 02.141334,
            02.151360, 02.161421, 02.171517, 02.181648, 02.191815,
            02.202018, 02.212257, 02.222533, 02.232845, 02.243195,
            02.253582, 02.264006, 02.274468, 02.284968, 02.295507,
            02.306084, 02.316701, 02.327356, 02.338051, 02.348786,
            02.359562, 02.370377, 02.381234, 02.392131, 02.403070,
            02.414051, 02.425073, 02.436138, 02.447246, 02.458397,
            02.469591, 02.480828, 02.492110, 02.503436, 02.514807,
            02.526222, 02.537684, 02.549190, 02.560743, 02.572343,
            02.583989, 02.595682, 02.607423, 02.619212, 02.631050,
            02.642936, 02.654871, 02.666855, 02.678890, 02.690975,
            02.703110, 02.715297, 02.727535, 02.739825, 02.752168,
            02.764563, 02.777012, 02.789514, 02.802070, 02.814681,
            02.827347, 02.840069, 02.852846, 02.865680, 02.878570,
            02.891518, 02.904524, 02.917588, 02.930712, 02.943894,
            02.957136, 02.970439, 02.983802, 02.997227, 03.010714,
            03.024263, 03.037875, 03.051551, 03.065290, 03.079095,
            03.092965, 03.106900, 03.120902, 03.134971, 03.149107,
            03.163312, 03.177585, 03.191928, 03.206340, 03.220824,
            03.235378, 03.250005, 03.264704, 03.279477, 03.294323,
            03.309244, 03.324240, 03.339312, 03.354461, 03.369687,
            03.384992, 03.400375, 03.415838, 03.431381, 03.447005,
            03.462711, 03.478500, 03.494372, 03.510328, 03.526370,
            03.542497, 03.558711, 03.575012, 03.591402, 03.607881,
            03.624450, 03.641111, 03.657863, 03.674708, 03.691646,
            03.708680, 03.725809, 03.743034, 03.760357, 03.777779,
            03.795300, 03.812921, 03.830645, 03.848470, 03.866400,
            03.884434, 03.902574, 03.920821, 03.939176, 03.957640,
            03.976215, 03.994901, 04.013699, 04.032612, 04.051639,
            04.070783, 04.090045, 04.109425, 04.128925, 04.148547,
            04.168292, 04.188160, 04.208154, 04.228275, 04.248524,
            04.268903, 04.289413, 04.310056, 04.330832, 04.351745,
            04.372794, 04.393982, 04.415310, 04.436781, 04.458395,
            04.480154, 04.502060, 04.524114, 04.546319, 04.568676,
            04.591187, 04.613854, 04.636678, 04.659662, 04.682807,
            04.706116, 04.729590, 04.753231, 04.777041, 04.801024,
            04.825179, 04.849511, 04.874020, 04.898710, 04.923582,
            04.948639, 04.973883, 04.999316, 05.024942, 05.050761,
            05.076778, 05.102993, 05.129411, 05.156034, 05.182864,
            05.209903, 05.237156, 05.264625, 05.292312, 05.320220,
            05.348354, 05.376714, 05.405306, 05.434131, 05.463193,
            05.492496, 05.522042, 05.551836, 05.581880, 05.612178,
            05.642734, 05.673552, 05.704634, 05.735986, 05.767610,
            05.799512, 05.831694, 05.864161, 05.896918, 05.929968,
            05.963316, 05.996967, 06.030925, 06.065194, 06.099780,
            06.134687, 06.169921, 06.205486, 06.241387, 06.277630,
            06.314220, 06.351163, 06.388465, 06.426130, 06.464166,
            06.502578, 06.541371, 06.580553, 06.620130, 06.660109,
            06.700495, 06.741297, 06.782520, 06.824173, 06.866262,
            06.908795, 06.951780, 06.995225, 07.039137, 07.083525,
            07.128398, 07.173764, 07.219632, 07.266011, 07.312910,
            07.360339, 07.408308, 07.456827, 07.505905, 07.555554,
            07.605785, 07.656608, 07.708035, 07.760077, 07.812747,
            07.866057, 07.920019, 07.974647, 08.029953, 08.085952,
            08.142657, 08.200083, 08.258245, 08.317158, 08.376837,
            08.437300, 08.498562, 08.560641, 08.623554, 08.687319,
            08.751955, 08.817481, 08.883916, 08.951282, 09.019600,
            09.088889, 09.159174, 09.230477, 09.302822, 09.376233,
            09.450735, 09.526355, 09.603118, 09.681054, 09.760191,
            09.840558, 09.922186, 10.005107, 10.089353, 10.174959,
            10.261958, 10.350389, 10.440287, 10.531693, 10.624646,
            10.719188, 10.815362, 10.913214, 11.012789, 11.114137,
            11.217307, 11.322352, 11.429325, 11.538283, 11.649285,
            11.762390, 11.877664, 11.995170, 12.114979, 12.237161,
            12.361791, 12.488946, 12.618708, 12.751161, 12.886394,
            13.024498, 13.165570, 13.309711, 13.457026, 13.607625,
            13.761625, 13.919145, 14.080314, 14.245263, 14.414134,
            14.587072, 14.764233, 14.945778, 15.131877, 15.322712,
            15.518470, 15.719353, 15.925570, 16.137345, 16.354912,
            16.578520, 16.808433, 17.044929, 17.288305, 17.538873,
            17.796967, 18.062943, 18.337176, 18.620068, 18.912049,
            19.213574, 19.525133, 19.847249, 20.180480, 20.525429,
            20.882738, 21.253102, 21.637266, 22.036036, 22.450278,
            22.880933, 23.329017, 23.795634, 24.281981, 24.789364,
            25.319207, 25.873062, 26.452634, 27.059789, 27.696581,
            28.365274, 29.068370, 29.808638, 30.589157, 31.413354,
            32.285060, 33.208568, 34.188705, 35.230920, 36.341388,
            37.527131, 38.796172, 40.157721, 41.622399, 43.202525,
            44.912465, 46.769077, 48.792279, 51.005773, 53.437996,
            56.123356, 59.103894
        };
    double U, V, RANLAN;
    int I;

    U = 1000.0 * X;
    I = U;
    U = U - I;

    if (I >= 70 && I <= 800)
    {
      RANLAN = F[I] + U * (F[I + 1] - F[I]);
    }
    else if (I >= 7 && I <= 980)
    {
      RANLAN = F[I]
          + U * (F[I + 1] - F[I]
              - 0.25 * (1 - U) * (F[I + 2] - F[I + 1] - F[I] + F[I - 1]));
    }
    else if (I < 7)
    {
      V = log(X);
      U = 1 / V;
      RANLAN = ((0.99858950 + (3.45213058E1 + 1.70854528E1 * U) * U) /
          (1 + (3.41760202E1 + 4.01244582 * U) * U)) *
          ( -log( -0.91893853 - V) - 1);
    }
    else
    {
      U = 1 - X;
      V = U * U;
      if (X <= 0.999)
      {
        RANLAN = (1.00060006 + 2.63991156E2 * U + 4.37320068E3 * V) /
            ((1 + 2.57368075E2 * U + 3.41448018E3 * V) * U);
      }
      else
      {
        RANLAN = (1.00001538 + 6.07514119E3 * U + 7.34266409E5 * V) /
            ((1 + 6.06511919E3 * U + 6.94021044E5 * V) * U);
      }
    }

    return RANLAN;
  }