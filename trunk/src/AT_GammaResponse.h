#ifndef SGP_GAMMARESPONSE_H_
#define SGP_GAMMARESPONSE_H_

#include <stdio.h>

extern int indent_counter;
extern char isp[];
extern FILE * debf;


void SGP_gamma_response(	long*	n,
							float*	d_Gy,
							long*	gamma_model,
							float*	gamma_parameter,
							// return
							float*	S);

void SGP_get_gamma_response(	long*	n,
								float*	d_Gy,
								float*	dd_Gy,
								float*	f,
								float*	f0,
								long*	gamma_model,
								float*	gamma_parameter,
								// return
								float*	S,
								float*	S_HCP,
								float*	S_gamma,
								float*	efficiency);

#endif // SGP_GAMMARESPONSE_H_
