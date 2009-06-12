#ifndef SGP_GAMMARESPONSE_H_
#define SGP_GAMMARESPONSE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void SGP_gamma_response(	long*	n,
							float*	d_Gy,
							long*	n_gamma_parameter,
							long*	gamma_model,
							float*	gamma_parameter,
							// return
							float*	S){
	long	i,j;
	/*
	 * (0) Testmodel, response m*x + c
	 * parameters:	m - gamma_parameter[0]
	 * 				c - gamma_parameter[1]
	 */
	if(*gamma_model == 0){
		float	m	= gamma_parameter[0];
		float	c	= gamma_parameter[1];
		for (i = 0; i < *n; i++){
			S[i]		=	m * d_Gy[i] + c;
		}
		return;
	}
	/*
	 *  (1) General multi-hit, multi-target model
	 *
	 *  N.B.: In this case the array length IS NOT (*n_parameter) but 4*(*n_parameter) !!
	 *
	 *	4 * i parameters:	k		- relative contribution from component i (= Smax_i), preferably adding to 1 for all components (not mandatory though!)
	 *	 					D1		- characteristic dose (one hit per target in average) of component i
	 * 						c		- hittedness of component i
	 * 						m		- number of targets for component i
	 */
	if(*gamma_model == 1){
		for (i = 0; i < *n; i++){
			S[i]	= 0.0f;
		}
		for (j = 0; j < *n_gamma_parameter; j++){						// loop over all components
			float k			=	gamma_parameter[j * 4];
			float D1		=	gamma_parameter[j * 4 + 1];
			float c			=	gamma_parameter[j * 4 + 2];
			float m			=	gamma_parameter[j * 4 + 3];

			float tmp		=	0.0f;

			for (i = 0; i < *n; i++){
				if(c == 1){										// in case of c = 1 use simplified, faster formula
					tmp			=	1.0f - exp(-1.0f * d_Gy[i] / D1);
				}else{											// for c > 1 use incomplete gamma function
					tmp			=	gammp(c, d_Gy[i] / D1);
				}

				if(m == 1){										// in case of m = 1 do not use pow for speed reasons
					tmp			*=	k;
				}else{
					tmp			= k * pow(tmp, m);
				}

				S[i]		+=	tmp;
			}
		}
		return;
	}
	/*
	 *  (2) RL ACCUMULATED COUNTS MODEL
	 *
	 *	parameters:	Smax 	- max. RATE signal (in contrast to model == 1, where Smax is the maximum response signal)
	 * 				chi		- dynamics, ration between Smax and start count rate c0
	 * 				D1		- characteristic dose at which rate signal reaches saturation
	 *
	 * 				are transformed into
	 *
	 * 				c0		- RATE signal at D = 0
	 * 				B		- linear slope of RATE signal below D = D1
	 * 				D1		- see above
	 */

	if(*gamma_model == 2){
		float Smax		=	gamma_parameter[0];
		float D1		=	gamma_parameter[1];
		float chi		=	gamma_parameter[2];
		// transform parameters
		float c0		=	Smax / chi;
		float B			=	(Smax - c0) / D1;

		for (i = 0; i < *n; i++){
			if(d_Gy[i]	<= D1){
				S[i]		=	c0 * d_Gy[i] + 0.5f * B * d_Gy[i] * d_Gy[i];}
			else{
				S[i]		=	c0 * D1 + 0.5f * B * D1 * D1 + Smax * (d_Gy[i] - D1);}
		}
		return;
	}

	/*
	 *  (3) EXPONENTIAL-SATURATION MODEL, obsolete
	 *
	 *	parameters:	Smax 	- max. response signal
	 * 				D1		- characteristic dose at which rate signal reaches saturation
	 */
	if(*gamma_model == 3){
		// exp-saturation model
		float	Smax	=	gamma_parameter[0];
		float	D0		=	gamma_parameter[1];

		for (i = 0; i < *n; i++){
			S[i]		=	Smax * (1.0f - (float)exp(-1.0f * d_Gy[i] / D0));
		}
		return;
	}
	return;
}

void SGP_get_gamma_response(	long*	n,
								float*	d_Gy,
								float*	dd_Gy,
								float*	f,
								float*	f0,
								long*	n_gamma_parameter,
								long*	gamma_model,
								float*	gamma_parameter,
								// return
								float*	S,
								float*	S_HCP,
								float*	S_gamma,
								float*	efficiency)
{
	long i;

	SGP_gamma_response(	n,
						d_Gy,
						n_gamma_parameter,
						gamma_model,
						gamma_parameter,
						// return
						S);

	*S_HCP			=	0.0f;
	float D_gamma	=	0.0f;

	for(i = 0; i < *n; i++){
			D_gamma		+=	d_Gy[i] * dd_Gy[i] * f[i];
			*S_HCP		+=	S[i] * dd_Gy[i] * f[i];
		}

	i	= 1;
	SGP_gamma_response(	&i,
						&D_gamma,
						n_gamma_parameter,
						gamma_model,
						gamma_parameter,
						// return
						S_gamma);

	*efficiency		=	*S_HCP / *S_gamma;
	return;
}

#endif // SGP_GAMMARESPONSE_H_
