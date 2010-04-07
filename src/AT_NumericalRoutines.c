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


double gammln(const double xx)
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
  return (-tmp+log(2.5066282746310005*ser/x));
}


#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30


void gcf(double *gammcf, const double a, const double x, double *gln)
{
  int i;
  double an,b,c,d,del,h;
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


void gser(double *gamser, const double a, const double x, double *gln)
{
  int n;
  double sum,del,ap;
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


double gammp(const double a, const double x)
{
  double gamser, gammcf, gln;

  if (x < 0.0 || a <= 0.0) return 0;
  if (x < (a + 1.0)) {
    gser(&gamser, a, x, &gln);
    return gamser;
  } else {
    gcf(&gammcf, a, x, &gln);
    return 1.0 - gammcf;
  }
}

// TODO in the standard math.h there is already implemented function erff which calculates error function
double erff_custom(const double x)
{
  return x < 0.0 ? -gammp(0.5, x*x) : gammp(0.5, x*x);
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
    nrerror("root must be bracketed in zriddr.");
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

    while ((strcmp( set[matches[i]], elements[i]) != 0) && (matches[i] < n_set)){
      matches[i]++;
    }

    if (matches[i] == n_set) {
      matches[i] = -1;
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


void is_element_int(const long element, const long set[], const long n_set, bool matches[]){

  long  i;
  for (i = 0; i < n_set; i++){
    if(element == set[i]){
      matches[i]  = true;
    } else{
      matches[i]  = false;
    }
  }
}


long locate(const double xx[], const long n, const double x)
{
  long  ju, jm, jl, j;
  int    ascnd;

  jl    =  0;
  ju    =  n + 1;
  ascnd  =  (xx[n-1] >= xx[1-1]);
  while (ju - jl > 1){
    jm    =  (ju + jl) >> 1;
    if (x >= xx[jm-1] == ascnd)
      jl  =  jm;
    else
      ju  =  jm;
  }
  if ( x == xx[1 - 1]) j = 1;
  else if (x == xx[n - 1]) j = n - 1;
  else j  =  jl;
  return j;
}


void polint(const double xa[], const double ya[], const long n, const double x, double *y, double *dy)
{
  long  i, m, ns=1;
  double  den, dif, dift, ho, hp, w;
  double  *c,*d;

  dif  =  fabs(x-xa[0]);
  c    =  (double*)calloc(n, sizeof(double));
  d    =  (double*)calloc(n, sizeof(double));
  for (i = 1; i <= n; i++) {
    if ( (dift = fabs(x - xa[i-1])) < dif) {
      ns    =  i;
      dif    =  dift;
    }
    c[i-1]  =  ya[i-1];
    d[i-1]  =  ya[i-1];
  }

  *y  =  ya[(ns--)-1];
  for (m = 1; m < n; m++) {
    for (i = 1; i <= n - m; i++) {
      ho  =  xa[i-1] - x;
      hp  =  xa[i+m-1] - x;
      w  =  c[i+1-1] - d[i-1];
      den  =  ho - hp;
      if ( den == 0.0) return;
      den  =  w / den;
      d[i-1]=  hp * den;
      c[i-1]=  ho * den;

    }
    // TODO do we really have "=" inside ?
    // underline code looks really magically
    *y += (*dy=(2*ns < (n-m) ? c[ns] : d[(ns--)-1]));
  }
  free(d);
  free(c);
}


void interp(const double xa[], const double ya[], const long n, const long n_pol, const double x, double y[], double dy[])
{
  long  j = locate(  xa,          // find index nearest to x
      n,
      x);
  long  k  =  GSL_MIN(GSL_MAX(j - (n_pol-1) / 2, 1), n + 1 - n_pol);
  polint(  &xa[k-1 -1],
      &ya[k-1 -1],
      n_pol,
      x,
      y,
      dy);
}
