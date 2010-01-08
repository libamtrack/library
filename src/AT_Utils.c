/**
 *    AT_Utils.c
 *    ==========
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

#include "AT_Utils.h"


void indnt_init(){
  debf = stderr;
  indent_counter = 0;
  char tmp[30] = "\0                            ";
  strcpy( isp , tmp);
};

void indnt_inc(){
   indent_counter++;
   isp[indent_counter] = '\0';
   isp[indent_counter-1] = ' ';
}

void indnt_dec(){
   indent_counter--;
   isp[indent_counter+1] = ' ';
   isp[indent_counter] = '\0';
}


// finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
// is reported per element
// a vector "matches" of length n_elements has to be provided
void pmatchi(long* elements, long* n_elements, long* set, long* n_set, long* matches){
  long  i;
  for (i = 0; i < *n_elements; i++){
    matches[i] = 0;

    while ((set[matches[i]] != elements[i]) && (matches[i] < *n_set))
    {
      matches[i]++;
    }

    if (matches[i] == *n_set) {matches[i] = -1;}
  }
}

// finds character elements in a set (n elements) and returns indices - only one (the first) match
// is reported per element
// a vector "matches" of length n_elements has to be provided
void pmatchc(char** elements, long* n_elements, char** set, long* n_set, long* matches){

  long  i;
  for (i = 0; i < *n_elements; i++){
    matches[i] = 0;

    while ((strcmp( set[matches[i]], elements[i]) != 0) && (matches[i] < *n_set))
    {
      matches[i]++;
    }

    if (matches[i] == *n_set) {matches[i] = -1;}
  }
}

// finds a character element in a set and returns boolean match vector
// a vector "matches" of length n_set has to be provided
void matchc(char* element, char** set, long* n_set, bool* matches){

  long  i;
  for (i = 0; i < *n_set; i++){
    if(strcmp(element, set[i])==0){
      matches[i]  = true;}
    else{
      matches[i]  = false;}
  }
}

// finds a integer element in a set and returns boolean match vector
// a vector "matches" of length n_set has to be provided
void matchi(long* element, long* set, long* n_set, bool* matches){

  long  i;
  for (i = 0; i < *n_set; i++){
    if(*element == set[i]){
      matches[i]  = true;}
    else{
      matches[i]  = false;}
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// interpolation on a table: code (w/ adapted indices) from Numerical Recipes, 2rd ed., chapter 3.1
// added wrapping function interp which allows to chose degree of interpolation polynomial
// (1 = linear, 2 = quadratic, etc.)
void locate(float* xx, long* n, float* x, long* j)
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

void polint(float* xa, float* ya, long* n, float* x, float *y, float *dy)
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

void interp(float* xa, float* ya, long* n, long* n_pol, float* x, float *y, float *dy)
{
  long  j;
  locate(  xa,          // find index nearest to x
      n,
      x,
      &j);
  long  k  =  LMIN(LMAX(j - (*n_pol-1) / 2, 1), *n + 1 - *n_pol);
  polint(  &xa[k-1 -1],
      &ya[k-1 -1],
      n_pol,
      x,
      y,
      dy);
  return;
}

// get LET-data for given material
void getPSTARvalue(long* n, float* x, long* material_no, float* x_table, float* y_table, float* y)
{
  // first: find those PSTAR entries that match the material name
  bool*    matches    =  (bool*)calloc(AT_PSTAR_Data.n, sizeof(bool));
  matchi(    material_no,
        AT_PSTAR_Data.material_no,
        &AT_PSTAR_Data.n,
        matches);

  long    n_matches  = 0;
  long    i;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){n_matches++;}
//    printf(debf,"idx: %i, match: %d\n",i, matches[i]);
  }

  // allocate vectors for extracted LET entries
  float*  x_c  =  (float*)calloc(n_matches, sizeof(float));
  float*  y_c  =  (float*)calloc(n_matches, sizeof(float));

  // and get the values
  long     j  = 0;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      x_c[j]  = x_table[i];
      y_c[j]  = y_table[i];
      j++;
    }
  }
  long  n_pol      = 4 + 1;
  for (i = 0; i < *n; i++){
    // Get proton-LET for scaled energy from table E, L using 4th degree polynomial (n_pol - 1 = 2) interpolation
    float  err_y_tmp  = 0.0f;    // dummy
    interp(    x_c,
          y_c,
          &n_matches,
          &n_pol,
          &x[i],
          &y[i],
          &err_y_tmp);
  }

  free(x_c);
  free(y_c);
  free(matches);
}


