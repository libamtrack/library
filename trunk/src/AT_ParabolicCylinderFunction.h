#ifndef SGP_PARABOLICCYLINDER_H_
#define SGP_PARABOLICCYLINDER_H_

#include <stdio.h>

extern int indent_counter;
extern char isp[];
extern FILE * debf;

void SGP_Funs(	float*	fz,
				float*	fR0,
				float*	fsigma,
				float*	fni,
				float*	funs);

void SGP_fDyx(	float*	fy,
				float*	fx,
				float*	fDyx);

double 	d_sign(	double *a, double *b);
int 	pbdv_(	double *v, double *x, double *dv, double *dp, double *pdf, double *pdd);
int 	dvsa_(	double *va, double *x, double *pd);
int 	dvla_(	double *va, double *x, double *pd);
int 	vvla_(	double *va, double *x, double *pv);
int 	gamma_(	double *x, double *ga);

#endif // SGP_PARABOLICCYLINDER_H_
