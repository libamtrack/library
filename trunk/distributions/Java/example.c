/* File : example.c */
#include "example.h"
 
int fact(int n) {
  double x = 5.0;
  double y = gsl_sf_bessel_J0(x);
  printf("J0(%g) = %.18e\n", x, y);
  printf("sqrt(%g) = %g\n", x, sqrt(x));
  if (n <= 1) return 1;
  else return n*fact(n-1);
}
 
int my_mod(int x, int y) {
  return (x%y);
}

char *get_time()
{
  time_t ltime;
  time(&ltime);
  return ctime(&ltime);
}
 
