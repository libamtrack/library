/**
 * @file
 * @brief Numerical routines
 */

/*
*    AT_NumericalRoutines.c
*    ==============
*
*    Created on: 8.01.2010
*    Author: kongruencja
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

#include "AT_NumericalRoutines.h"

long int lminl(const long int x, const long int y)
{
  return (x < y) ? x : y;
}

long int lmaxl(const long int x, const long int y)
{
  return (x > y) ? x : y;
}


/*       ========================================================= */
/*       Purpose: This program computes the parabolic cylinder */
/*                functions Dv(x) and their derivatives using */
/*                subroutine PBDV */
/*       Input:   x --- Argument of Dv(x) */
/*                v --- Order of Dv(x) */
/*       Output:  DV(na) --- Dn+v0(x) */
/*                DP(na) --- Dn+v0'(x) */
/*                ( na = |n|, n = int(v), v0 = v-n, |v0| < 1 */
/*                  n = 0,+_1,+-2,..., |n| ? 100 ) */
/*                PDF --- Dv(x) */
/*                PDD --- Dv'(x) */
/*       Example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,...,5 */

/*                  n+v0      Dv(x)           Dv'(x) */
/*                --------------------------------------- */
/*                  0.5   .43971930D-10  -.21767183D-09 */
/*                  1.5   .43753148D-09  -.21216995D-08 */
/*                  2.5   .43093569D-08  -.20452956D-07 */
/*                  3.5   .41999741D-07  -.19491595D-06 */
/*                  4.5   .40491466D-06  -.18355745D-05 */
/*                  5.5   .38601477D-05  -.17073708D-04 */

/*                Dv(x)= .38601477D-05,  Dv'(x)=-.17073708D-04 */
/*       ========================================================= */
/*       FORTRAN code by Jianming Jin, Department of Electrical and */
/*       Computer Engineering, University of Illinois at Urbana-Champaign */
/*       http://jin.ece.uiuc.edu/routines/mpbdv.for */
/*       Downloaded and converted to C using f2c by */
/*       S. Greilich, Risoe National Laboratory, Denmark */
/*       Reworked as subroutine for AmTrack.dll, abandoning */
/*       f2c.h and libf2c.lib */
/*       ========================================================= */

/* mpbdv.f -- translated by f2c (version 20060506).*/

/* Table of constant values */

//static int c__9 = 9;
//static int c__1 = 1;
//static int c__5 = 5;
static double c_b31 = 1.;
static double c_b40 = 2.;

void AT_Dyx(  double*  y,  double*  x,  double*  Dyx)
{
  /* Local variables */
  //    static int na;
  static double dp[101], dv[101];
  static double pdd, pdf;
  extern /* Subroutine */ int pbdv_(double *, double *, double *
      , double *, double *, double *);

  pbdv_(y, x, dv, dp, &pdf, &pdd);

  *Dyx  =  pdf;
}

// TODO this is not used !
void AT_fDyx(  const float*  fy,
    const float* fx,
    float* fDyx)
{
  double  y, x, Dyx;
  y  = (double)*fy;
  x  = (double)*fx;
  AT_Dyx(  &y, &x, &Dyx);
  *fDyx  = Dyx;
}


double d_sign(const double *a, const double *b)
{
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
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
  extern /* Subroutine */ int dvla_(double *, double *, double *
  ), dvsa_(double *, double *, double *);


  /*       ==================================================== */
  /*       Purpose: Compute parabolic cylinder functions Dv(x) */
  /*                and their derivatives */
  /*       Input:   x --- Argument of Dv(x) */
  /*                v --- Order of Dv(x) */
  /*       Output:  DV(na) --- Dn+v0(x) */
  /*                DP(na) --- Dn+v0'(x) */
  /*                ( na = |n|, v0 = v-n, |v0| < 1, */
  /*                  n = 0,�1,�2,��� ) */
  /*                PDF --- Dv(x) */
  /*                PDD --- Dv'(x) */
  /*       Routines called: */
  /*             (1) DVSA for computing Dv(x) for small |x| */
  /*             (2) DVLA for computing Dv(x) for large |x| */
  /*       ==================================================== */

  xa = fabs(*x);
  vh = *v;
  *v += d_sign(&c_b31, v);
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



int dvsa_(double *va,
    double *x,
    double *pd)
{
  /* System generated locals */
  double d__1;

  /* Local variables */
  static int m;
  static double r__, a0, g0, g1, r1, ep, gm, pi, vm, vt, ga0, va0, sq2,
  eps;

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
        gamma_(&va0, &ga0);
        d__1 = *va * -.5;
        *pd = sqrt(pi) / (pow(c_b40, d__1) * ga0);
      }
    } else {
      d__1 = -(*va);
      gamma_(&d__1, &g1);
      d__1 = *va * -.5 - 1.;
      a0 = pow(c_b40, d__1) * ep / g1;
      vt = *va * -.5;
      gamma_(&vt, &g0);
      *pd = g0;
      r__ = 1.;
      for (m = 1; m <= 250; ++m) {
        vm = (m - *va) * .5;
        gamma_(&vm, &gm);
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
    gamma_(&d__1, &gl);
    *pd = pi * vl / gl + cos(pi * *va) * *pd;
  }
  return 0;
} /* dvla_ */


int vvla_(double *va, double *x, double *pv)
{
  /* System generated locals */
  double d__1, d__2;

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
    gamma_(&d__1, &gl);
    dsl = sin(pi * *va) * sin(pi * *va);
    *pv = dsl * gl / pi * pdl - cos(pi * *va) * *pv;
  }
  return 0;
} /* vvla_ */


int gamma_(const double *x, double *ga)
{
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


//TODO this is not used anywhere
void AT_Funs(  const float*  fz,
    const float*  fR0,
    const float*  fsigma,
    const float* fni,
    float* funs)
{
  double  z    =  (double)*fz;
  double  R0    =  (double)*fR0;
  double  sigma  =  (double)*fsigma;
  double  ni    =  (double)*fni + 1.0;

  double  u    =  R0 - z;
  double  zeta  =  u / sigma;

  *funs = 0.0f;

  if(zeta > -5.0 && zeta < 10){
    double  tmp1  =  1.0 / (sqrt(2.0 * M_PI) * sigma);
    double  tmp2  =  exp(-1.0 * u * u / (4.0 * sigma * sigma)) * pow(sigma, ni);
    double  tmp3;
    gamma_(&ni, &tmp3);
    double  tmp4;
    double  y    =  -1.0 * ni;
    double  x    =  -1.0 * zeta;
    AT_Dyx(&y, &x, &tmp4);

    double  result  =  tmp1 * tmp2 * tmp3 * tmp4;
    *funs      =  (float)result;
  }

  if(zeta >= 10.0){
    *funs      =  pow(u, ni - 1.0);
  }
}


float gammln(const float xx)
{
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
  return (float)(-tmp+log(2.5066282746310005*ser/x));
}


#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(float *gammcf, const float a, const float x, float *gln)
{
  int i;
  float an,b,c,d,del,h;
  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  //  if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}


void gser(float *gamser, const float a, const float x, float *gln)
{
  int n;
  float sum,del,ap;
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) return;
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
        *gamser=sum*exp(-x+a*log(x)-(*gln));
        return;
      }
    }
    //  nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}


float gammp(const float a, const float x)
{
  float gamser, gammcf, gln;

  if (x < 0.0f || a <= 0.0f) return 0;
  if (x < (a + 1.0f)) {
    gser(&gamser, a, x, &gln);
    return gamser;
  } else {
    gcf(&gammcf, a, x, &gln);
    return 1.0f - gammcf;
  }
}


float erff(const float x)
{
  return x < 0.0f ? -gammp(0.5f, x*x) : gammp(0.5f, x*x);
}


void nrerror(const char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}


#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAXIT 60
#define UNUSED (-1.11e30)


float zriddr(float (*func)(float,void*), void * params, const float x1, const float x2, const float xacc)
{
  int j;
  float ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
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
    nrerror("root must be bracketed in zriddr.");
  }
  return 0.0;                         // Never get here.
}


void are_elements_int(const int* elements, const int n_elements, const int* set, const int n_set, int* matches){
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


void find_elements_char(const char** elements, const long* n_elements, const char* const * set, const long* n_set, long* matches){

  long  i;
  for (i = 0; i < *n_elements; i++){
    matches[i] = 0;

    while ((strcmp( set[matches[i]], elements[i]) != 0) && (matches[i] < *n_set)){
      matches[i]++;
    }

    if (matches[i] == *n_set) {
      matches[i] = -1;
    }
  }
}


void is_element_char(const char* element, const char* const * set, const long* n_set, bool* matches){

  long  i;
  for (i = 0; i < *n_set; i++){
    if(strcmp(element, set[i])==0){
      matches[i]  = true;
    }else{
      matches[i]  = false;
    }
  }
}


void is_element_int(const long* element, const long* set, const long* n_set, bool* matches){

  long  i;
  for (i = 0; i < *n_set; i++){
    if(*element == set[i]){
      matches[i]  = true;
    } else{
      matches[i]  = false;
    }
  }
}


void locate(const float* xx, const long* n, const float* x, long* j)
{
  long  ju, jm, jl;
  int    ascnd;

  jl    =  0;
  ju    =  *n + 1;
  ascnd  =  (xx[*n-1] >= xx[1-1]);
  while (ju - jl > 1){
    jm    =  (ju + jl) >> 1;
    if (*x >= xx[jm-1] == ascnd)
      jl  =  jm;
    else
      ju  =  jm;
  }
  if ( *x == xx[1 - 1]) *j = 1;
  else if (*x == xx[*n - 1]) *j = *n - 1;
  else *j  =  jl;
  return;
}


void polint(const float* xa, const float* ya, const long* n, const float* x, float *y, float *dy)
{
  long  i, m, ns=1;
  float  den, dif, dift, ho, hp, w;
  float  *c,*d;

  dif    =  (float)fabs(*x-xa[1-1]);
  c    =  (float*)calloc(*n, sizeof(float));
  d    =  (float*)calloc(*n, sizeof(float));
  for (i = 1; i <= *n; i++) {
    if ( (dift = (float)fabs(*x - xa[i-1])) < dif) {
      ns    =  i;
      dif    =  dift;
    }
    c[i-1]  =  ya[i-1];
    d[i-1]  =  ya[i-1];
  }

  *y  =  ya[(ns--)-1];
  for (m = 1; m < *n; m++) {
    for (i = 1; i <= *n - m; i++) {
      ho  =  xa[i-1] - *x;
      hp  =  xa[i+m-1] - *x;
      w  =  c[i+1-1] - d[i-1];
      den  =  ho - hp;
      if ( den == 0.0) return;
      den  =  w / den;
      d[i-1]=  hp * den;
      c[i-1]=  ho * den;

    }
    *y += (*dy=(2*ns < (*n-m) ? c[ns+1-1] : d[(ns--)-1]));
  }
  free(d);
  free(c);
}


//TODO change pointer to single variable where necessary: n_pol for example
void interp(const float* xa, const float* ya, const long* n, const long* n_pol, const float* x, float *y, float *dy)
{
  long  j;
  locate(  xa,          // find index nearest to x
      n,
      x,
      &j);
  long  k  =  lminl(lmaxl(j - (*n_pol-1) / 2, 1), *n + 1 - *n_pol);
  polint(  &xa[k-1 -1],
      &ya[k-1 -1],
      n_pol,
      x,
      y,
      dy);
  return;
}
