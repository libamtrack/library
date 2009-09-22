/**
 *    AT_ParabolicCylinderFunction.c
 *    ==============================
 *
 *    Created on: 28.07.2009
 *    Author: greilich
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

#include "AT_ParabolicCylinderFunction.h"

#include "AT_Constants.h"

#include <math.h>
#include <stdlib.h>

/*       ========================================================= */
/*       Purpose: This program computes the parabolic cylinder */
/*                functions Dv(x) and their derivatives using */
/*                subroutine PBDV */
/*       Input:   x --- Argument of Dv(x) */
/*                v --- Order of Dv(x) */
/*       Output:  DV(na) --- Dn+v0(x) */
/*                DP(na) --- Dn+v0'(x) */
/*                ( na = |n|, n = int(v), v0 = v-n, |v0| < 1 */
/*                  n = 0,�1,�2,���, |n| � 100 ) */
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

void AT_fDyx(  float*  fy,  float*  fx,  float*  fDyx)
{
  double  y, x, Dyx;
  y  = (double)*fy;
  x  = (double)*fx;

  AT_Dyx(  &y, &x, &Dyx);

  *fDyx  = Dyx;
}


double d_sign(double *a, double *b)
{
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}


int pbdv_(double *v, double *x, double *dv,
  double *dp, double *pdf, double *pdd)
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

int dvsa_(double *va, double *x, double *pd)
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    static int m;
    static double r__, a0, g0, g1, r1, ep, gm, pi, vm, vt, ga0, va0, sq2,
      eps;
    extern /* Subroutine */ int gamma_(double *, double *);


/*       =================================================== */
/*       Purpose: Compute parabolic cylinder function Dv(x) */
/*                for small argument */
/*       Input:   x  --- Argument */
/*                va --- Order */
/*       Output:  PD --- Dv(x) */
/*       Routine called: GAMMA for computing �(x) */
/*       =================================================== */

    eps = 1e-15;
    pi = 3.141592653589793;
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
    extern /* Subroutine */ int vvla_(double *, double *, double *
      ), gamma_(double *, double *);


/*       ==================================================== */
/*       Purpose: Compute parabolic cylinder functions Dv(x) */
/*                for large argument */
/*       Input:   x  --- Argument */
/*                va --- Order */
/*       Output:  PD --- Dv(x) */
/*       Routines called: */
/*             (1) VVLA for computing Vv(x) for large |x| */
/*             (2) GAMMA for computing �(x) */
/*       ==================================================== */

    pi = 3.141592653589793;
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
    extern /* Subroutine */ int dvla_(double *, double *, double *
      ), gamma_(double *, double *);


/*       =================================================== */
/*       Purpose: Compute parabolic cylinder function Vv(x) */
/*                for large argument */
/*       Input:   x  --- Argument */
/*                va --- Order */
/*       Output:  PV --- Vv(x) */
/*       Routines called: */
/*             (1) DVLA for computing Dv(x) for large |x| */
/*             (2) GAMMA for computing �(x) */
/*       =================================================== */

    pi = 3.141592653589793;
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

int gamma_(double *x, double *ga)
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


/*       ================================================== */
/*       Purpose: Compute gamma function �(x) */
/*       Input :  x  --- Argument of �(x) */
/*                       ( x is not equal to 0,-1,-2,���) */
/*       Output:  GA --- �(x) */
/*       ================================================== */

    pi = 3.141592653589793;
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



void AT_Funs(  float*  fz,  float*  fR0,  float*  fsigma, float* fni, float* funs)
{
  double  z    =  (double)*fz;
  double  R0    =  (double)*fR0;
  double  sigma  =  (double)*fsigma;
  double  ni    =  (double)*fni + 1.0;

  double  u    =  R0 - z;
  double  zeta  =  u / sigma;

  *funs = 0.0f;

  if(zeta > -5.0 && zeta < 10){
    double  tmp1  =  1.0 / (sqrt(2.0 * pi) * sigma);
    double  tmp2  =  exp(-1.0 * u * u / (4.0 * sigma * sigma)) * pow(sigma, ni);
    double  tmp3;
    gamma_(&ni, &tmp3);
    double  tmp4;
    double  y    =  -1.0 * ni;
    double  x    =  -1.0 * zeta;
    AT_Dyx(&y, &x, &tmp4);

    double  result  =  tmp1 * tmp2 * tmp3 * tmp4;
    *funs      =  (float)result;}

  if(zeta >= 10.0){
    *funs      =  pow(u, ni - 1.0);}
}

