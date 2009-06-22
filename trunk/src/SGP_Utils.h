#ifndef SGP_UTILS_H_
#define SGP_UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <malloc.h>
#include "SGP_Constants.h"
#include "SGP_Data.h"
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

//
#ifdef _DEBUG
int indent_counter = 0;
char isp[30] = "\0                            ";
FILE * debf;

void indnt_inc();
void indnt_dec();
void indnt_init();

void indnt_init(){
	debf = stderr;
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
#endif


// finds integer (32bit) elements in a set (n elements) and returns indices - only one (the first) match
// is reported per element
// a vector "matches" of length n_elements has to be provided
void pmatchi(long* elements, long* n_elements, long* set, long* n_set, long* matches){
//	printf("n_elements = %ld\n", *n_elements);
//	printf("n_set = %ld\n", *n_set);
//	long ii;
//	for( ii = 0 ; ii < *n_elements ; ii++){
//		printf("elements[%ld]=%ld\n", ii , elements[ii]);
//		printf("matches[%ld]=%ld\n", ii , matches[ii]);
//	}
//	for( ii = 0 ; ii < *n_set ; ii++){
//		printf("set[%ld]=%ld\n", ii , set[ii]);
//	}

	long	i;
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

	long	i;
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

	long	i;
	for (i = 0; i < *n_set; i++){
		if(strcmp(element, set[i])==0){
			matches[i]	= true;}
		else{
			matches[i]	= false;}
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// interpolation on a table: code (w/ adapted indices) from Numerical Recipes, 2rd ed., chapter 3.1
// added wrapping function interp which allows to chose degree of interpolation polynomial
// (1 = linear, 2 = quadratic, etc.)
void locate(float* xx, long* n, float* x, long* j)
{
	long	ju, jm, jl;
	int		ascnd;

	jl		=	0;
	ju		=	*n + 1;
	ascnd	=	(xx[*n-1] >= xx[1-1]);
	while (ju - jl > 1){
		jm		=	(ju + jl) >> 1;
		if (*x >= xx[jm-1] == ascnd)
			jl	=	jm;
		else
			ju	=	jm;
	}
	if ( *x == xx[1 - 1]) *j = 1;
	else if (*x == xx[*n - 1]) *j = *n - 1;
	else *j	=	jl;
	return;
}

void polint(float* xa, float* ya, long* n, float* x, float *y, float *dy)
{
	long	i, m, ns=1;
	float	den, dif, dift, ho, hp, w;
	float	*c,*d;

	dif		=	(float)fabs(*x-xa[1-1]);
	c		=	(float*)calloc(*n, sizeof(float));
	d		=	(float*)calloc(*n, sizeof(float));
	for (i = 1; i <= *n; i++) {
		if ( (dift = (float)fabs(*x - xa[i-1])) < dif) {
			ns		=	i;
			dif		=	dift;
		}
		c[i-1]	=	ya[i-1];
		d[i-1]	=	ya[i-1];
	}

	*y	=	ya[(ns--)-1];
	for (m = 1; m < *n; m++) {
		for (i = 1; i <= *n - m; i++) {
			ho	=	xa[i-1] - *x;
			hp	=	xa[i+m-1] - *x;
			w	=	c[i+1-1] - d[i-1];
			den	=	ho - hp;
			if ( den == 0.0) return;
			den	=	w / den;
			d[i-1]=	hp * den;
			c[i-1]=	ho * den;

		}
		*y += (*dy=(2*ns < (*n-m) ? c[ns+1-1] : d[(ns--)-1]));
	}
	free(d);
	free(c);
}

void interp(float* xa, float* ya, long* n, long* n_pol, float* x, float *y, float *dy)
{
	long	j;
	locate(	xa,					// find index nearest to x
			n,
			x,
			&j);
	long	k	=	LMIN(LMAX(j - (*n_pol-1) / 2, 1), *n + 1 - *n_pol);
	polint(	&xa[k-1 -1],
			&ya[k-1 -1],
			n_pol,
			x,
			y,
			dy);
	return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]=	{	76.18009172947146,-86.50532032941677,
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

void gcf(float *gammcf, float a, float x, float *gln)
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
//	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void gser(float *gamser, float a, float x, float *gln)
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
//	nrerror("a too large, ITMAX too small in routine gser");
	return;
	}
}

float gammp(float a, float x)
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


float erff(float x)
{
	return x < 0.0f ? -gammp(0.5f, x*x) : gammp(0.5f, x*x);
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAXIT 60
#define UNUSED (-1.11e30)
float zriddr(float (*func)(float,void*), void * params, float x1, float x2, float xacc)
/* 	From Numerical Recipes in C, 2nd ed., 1992:
	Using Ridders' method, return the root of a function func known to lie between x1 and x2.
	The root, returned as zriddr, will be refined to an approximate accuracy xacc.
 */
{
	int j;
	float ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
	fl=(*func)(x1,params);
	fh=(*func)(x2,params);
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl=x1;
		xh=x2;
		ans=UNUSED; 												//	Any highly unlikely value, to simplify logic below.
		for (j=1;j<=MAXIT;j++) {
			xm=0.5*(xl+xh);
			fm=(*func)(xm,params); 										// First of two function evaluations per its=
			s=sqrt(fm*fm-fl*fh); 									// eration.
			if (s == 0.0) return ans;
			xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s); 	// Updating formula.
			if (fabs(xnew-ans) <= xacc) return ans;
			ans=xnew;
			fnew=(*func)(ans,params); 									// Second of two function evaluations per
			if (fnew == 0.0) return ans; 						// iteration.
			if (SIGN(fm,fnew) != fm) { 							// Bookkeeping to keep the root bracketed
				xl=xm; 											// on next iteration.
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
	return 0.0; 												// Never get here.
}


#endif // SGP_UTILS_H_
