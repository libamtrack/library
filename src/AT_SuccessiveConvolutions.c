/**
 * @file
 * @brief ...
 */

/*
 *    AT_SuccessiveConvolutions.c
 *    ===========================
 *
 *    Created on: 28.07.2009
 *    Creator: greilich
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

#include "AT_SuccessiveConvolutions.h"

#define MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT  0.002
#define DEBUG_INTERVALS                      8
#define DEBUG_MEAN                           10.0
#define DEBUG_SIGMA                          1.0


long  AT_n_bins_for_singe_impact_local_dose_distrib(
		const long    n,
		const double  E_MeV_u[],
		const long    particle_no[],
		const long    material_no,
		const long    rdd_model,
		const double  rdd_parameter[],
		const long    er_model,
		const long    N2)
{
	/* get lowest and highest dose */
	double d_max_Gy    =  0.0;
	double d_min_Gy    =  0.0;

	// TODO think if d_min calculations can be done in smarter way. LET is only needed for Geiss RDD

	long  i;
	for (i = 0; i < n; i++){

		double max_electron_range_m = AT_max_electron_range_m( E_MeV_u[i], (int)material_no, (int)er_model);
		double LET_MeV_cm2_g        = AT_LET_MeV_cm2_g_single(E_MeV_u[i], particle_no[i], material_no);
		double norm_constant_Gy     = AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model);
		double current_d_min_Gy     = AT_RDD_d_min_Gy( E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model, norm_constant_Gy);
		double current_d_max_Gy     = AT_RDD_d_max_Gy( E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model);

		if(i == 0){
			d_min_Gy      =  current_d_min_Gy;
			d_max_Gy      =  current_d_max_Gy;
		}
		else{
			d_min_Gy      =  GSL_MIN(d_min_Gy, current_d_min_Gy);
			d_max_Gy      =  GSL_MAX(d_max_Gy, current_d_max_Gy);
		}
	}

	long n_bins_for_singe_impact_local_dose_ditrib = 0;
	// get number of bins needed to span that dose range
	if( (d_min_Gy > 0) && (d_max_Gy >0) ){
		double tmp        =  log10(d_max_Gy/d_min_Gy) / log10(2.0) * ((double)N2);
		n_bins_for_singe_impact_local_dose_ditrib        =  (long)(floor(tmp) + 1.0);
	} else {
		printf("AT_n_bins_for_singe_impact_local_dose_ditrib: problem in evaluating n_bins_for_singe_impact_local_dose_ditrib: d_min = %g [Gy], d_max = %g [Gy] \n", d_min_Gy, d_max_Gy);
		exit(EXIT_FAILURE);
	}
	return n_bins_for_singe_impact_local_dose_ditrib;
}


void  AT_single_impact_local_dose_distrib(
		const long    n,
		const double  E_MeV_u[],
		const long    particle_no[],
		const double  fluence_cm2_or_dose_Gy[],
		const long    material_no,
		const long    rdd_model,
		const double  rdd_parameter[],
		const long    er_model,
		const long    N2,
		const long    n_bins_f1,
		const double  f1_parameters[],
		double        f1_d_Gy[],
		double        f1_dd_Gy[],
		double        f1[])
{
	double*  fluence_cm2    =  (double*)calloc(n, sizeof(double));

	long i;
	if(fluence_cm2_or_dose_Gy[0] < 0){
		double*  dose_Gy        =  (double*)calloc(n, sizeof(double));
		for (i = 0; i < n; i++){
			dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
		}
		AT_fluence_cm2(  n,
				E_MeV_u,
				particle_no,
				dose_Gy,
				material_no,
				fluence_cm2);
		free( dose_Gy );
	}else{
		for (i = 0; i < n; i++){
			fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
		}
	}
	double*  norm_fluence                                 =  (double*)calloc(n, sizeof(double));

	// Normalize fluence vector
	AT_normalize(    n,
			fluence_cm2,
			norm_fluence);

	free( fluence_cm2 );

	if(n_bins_f1 > 0){
		double  d_min      =  f1_parameters[0*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3];
		double  d_max      =  f1_parameters[0*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4];

		for (i = 1; i < n; i++){
			d_min          =  GSL_MIN(f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3], d_min);
			d_max          =  GSL_MAX(f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4], d_max);
		}

		double  U        =  (log(2.0) / (double)N2);

		double*  d_df_low      =  (double*)calloc(n_bins_f1, sizeof(double));
		double*  d_df_mid      =  (double*)calloc(n_bins_f1, sizeof(double));
		double*  d_df_high     =  (double*)calloc(n_bins_f1, sizeof(double));
		double*  dd_df         =  (double*)calloc(n_bins_f1, sizeof(double));

		for (i = 0; i < n_bins_f1; i++){
			// TODO: check if n.bins sufficient

			d_df_low[i]          =   d_min * exp((double)i * U);
			d_df_high[i]         =   d_min * exp(((double)i + 1.0) * U);
			d_df_mid[i]          =   d_min * exp(((double)i + 0.5) * U);
			dd_df[i]             =   d_df_high[i] - d_df_low[i];              // OBS: using Kellerer's bin-width = mid.point * U is not entirely correct

			f1[i]            =  0.0;
		}

		long n_bins_used = 1;

		// loop over all particles and energies, compute contribution to f1
		long   k;
		for (k = 0; k < n; k++){

			double  d_min_k        =  f1_parameters[k*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3];
			double  d_max_k        =  f1_parameters[k*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4];

			// find first and last bin to fit this particle's contribution into the all-over f1-frame
			long  i_low, i_high;
			i_low  = locate(d_df_low, n_bins_f1, d_min_k);
			i_high = locate(d_df_low, n_bins_f1, d_max_k);
			i_low            -=  1;
			i_high           -=  1;

			long  n_bins_df      =  i_high - i_low + 1;  // changed from + 2

			if (n_bins_df > 1){
				double*  d_low       =  (double*)calloc(n_bins_df, sizeof(double));
				double*  d_mid       =  (double*)calloc(n_bins_df, sizeof(double));
				double*  d_high      =  (double*)calloc(n_bins_df, sizeof(double));
				double*  dd          =  (double*)calloc(n_bins_df, sizeof(double));
				double*  r           =  (double*)calloc(n_bins_df, sizeof(double));
				double*  F1_1        =  (double*)calloc(n_bins_df, sizeof(double));
				double*  f1_k        =  (double*)calloc(n_bins_df - 1, sizeof(double));

				// extract the corresponding part from the all-over frame
				for (i = 0; i < n_bins_df; i++){
					d_low[i]           =   d_df_low[i_low + i];
					d_high[i]          =   d_df_high[i_low + i];
					d_mid[i]           =   d_df_mid[i_low + i];
					dd[i]              =   d_high[i] - d_low[i];
				};

				// and adjust the edges
				d_low[0]                =  d_min_k;
				d_low[n_bins_df-1]      =  d_max_k;

				d_mid[0]                =  sqrt(d_low[0] * d_high[0]);
				d_mid[n_bins_df-1-1]    =  sqrt(d_low[n_bins_df - 1] * d_high[n_bins_df - 1]);
				d_mid[n_bins_df-1]      =  0.0;

				d_high[n_bins_df-1-1]   =  d_max_k;
				d_high[n_bins_df-1]     =  0.0;  //TODO ??

				dd[n_bins_df-1]         =  0.0;

				// now compute r, F1, and f1, this could be any RDD if implemented
				int inverse_RDD_status_code = AT_r_RDD_m  (  n_bins_df,
						d_low,
						E_MeV_u[k],
						particle_no[k],
						material_no,
						rdd_model,
						rdd_parameter,
						er_model,
						r);

				if( inverse_RDD_status_code != 0 ){
					printf("Problem in evaluating inverse RDD in AT_SC_get_f1, probably wrong combination of ER and RDD used\n");
					char rdd_model_name[100];
					AT_RDD_name_from_number(rdd_model, rdd_model_name);
					char er_model_name[100];
					getERName( er_model, er_model_name);
					printf("rdd_model: %ld (%s), er_model: %ld (%s)\n", rdd_model, rdd_model_name, er_model, er_model_name);
					exit(EXIT_FAILURE);
				}

				for (i = 0; i < n_bins_df; i++){
					F1_1[i]            = gsl_pow_2(r[i] / f1_parameters[k*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 2]);        // F1 - 1 instead of F1 to avoid numeric cut-off problems
				}

				F1_1[n_bins_df-1]    =  0.0;

				// now compute f1 as the derivative of F1
				for (i = 0; i < (n_bins_df - 1); i++){
					f1_k[i]          =  -1.0 * (F1_1[i + 1] - F1_1[i]) / dd[i];
				}

				// adjust the density in first and last bin, because upper limit is not d.max.Gy and lower not d.min.Gy
				f1_k[0]              *=  dd[0] / dd_df[i_low];
				f1_k[n_bins_df-2]    *=  dd[n_bins_df - 2] / dd_df[i_high - 1];

				// and paste f1 for this energy /particle into the over all data frame according to rel. fluence
				for (i = 0; i < (n_bins_df - 1); i++){
					f1[i_low + i]      +=  norm_fluence[k] * f1_k[i];
				}

				free(d_low);
				free(d_mid);
				free(d_high);
				free(dd);
				free(r);
				free(F1_1);
				free(f1_k);
			}
			else{ // n_bins_df == 1
				f1[i_low ]        +=  norm_fluence[k] * 1.0 / dd_df[i_low];
			}

			// remember highest bin used
			n_bins_used          =  GSL_MAX(n_bins_used, i_high);
		}

		// copy back for the dose axis
		for (i = 0; i < n_bins_f1; i++){
			f1_d_Gy[i]    =  d_df_mid[i];
			f1_dd_Gy[i]   =  dd_df[i];
		}

		free(d_df_low);
		free(d_df_mid);
		free(d_df_high);
		free(dd_df);

		// normalize f1 (should be ok anyway but there could be small round-off errors)
		double  f1_norm    =  0.0;
		for (i = 0; i < n_bins_f1; i++){
			f1_norm    +=    f1[i] * f1_dd_Gy[i];
		}
		for (i = 0; i < n_bins_f1; i++){
			f1[i]    /=    f1_norm;
		}

	} // if(f1_d_Gy != NULL)

	free( norm_fluence );

}


void  AT_n_bins_for_low_fluence_local_dose_distribution( const double  u,
		const double   fluence_factor,
		const long     N2,
		const long     n_bins_f1,
		const double   f1_d_Gy[],
		const double   f1_dd_Gy[],
		const double   f1[],
		long*          n_bins_f,
		double*        u_start,
		long*          n_convolutions)
{
	// Get expectation value of dose from f1
	double  d_f1_Gy    =  0.0;

	long   i;
	for (i = 0; i < n_bins_f1; i++){
		d_f1_Gy    +=  f1_d_Gy[i] * f1_dd_Gy[i] * f1[i];
	}

	// The dose set by the input data is therefore
	double  d_f_Gy    = u * fluence_factor * d_f1_Gy;

	// How many convolution are necessary starting from a small mean
	// impact number that allows linear approximation (e.g. 0.002)
	*n_convolutions    = 0;

	*u_start      =    d_f_Gy  / d_f1_Gy;
	while(*u_start  > MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT){
		*u_start      =    0.5 * (*u_start);
		(*n_convolutions)++;
	}

	// Every convolution will add a factor of two, so the array for f has to
	// be expanded from f1 by N2 * n_convolutions
	*n_bins_f      =  (*n_convolutions + 1) * N2;
	*n_bins_f      +=  n_bins_f1;
}


void  AT_low_fluence_local_dose_distribution(  const long     n_bins_f1,
		const long     N2,
		const double   f1_d_Gy[],
		const double   f1_dd_Gy[],
		const double   f1[],
		const long     n_bins_f,
		double         f_d_Gy[],
		double         f_dd_Gy[],
		double         f_start[])
{
	// temporary arrays
	double*  d_low         =  (double*)calloc(n_bins_f, sizeof(double));
	double*  d_high        =  (double*)calloc(n_bins_f, sizeof(double));

	double  U              =  log(2.0) / (double)N2;

	long   i;
	for (i = 0; i < n_bins_f; i++){
		d_low[i]            =   f1_d_Gy[0] * exp(((double)i - 0.5)* U);
		d_high[i]           =   f1_d_Gy[0] * exp(((double)i + 0.5)* U);
		if (i < n_bins_f1){
			f_d_Gy[i]         =  f1_d_Gy[i];
			f_dd_Gy[i]        =  f1_dd_Gy[i];
			f_start[i]        =  f1[i];
		}else{
			f_d_Gy[i]         =  sqrt(d_low[i] * d_high[i]);
			f_dd_Gy[i]        =  d_high[i] - d_low[i];
			f_start[i]        =  0.0;
		}
	}

	free(d_low);
	free(d_high);
}


void AT_Kellerer_normalize(const long array_size,
		const long MIH,
		const long MIE,
		const long LEH,
		const double H0,
		const double E[],
		const double DE[],
		double* CM1,
		double* CM2,
		double* X,
		double H[]){

	double  Y        =  (*CM1) * 2;
	double  Z        =  (*CM2) * 2;
	double  CM0      =  H0;
	*CM1     =  0;

	long    N        =  MIH - MIE;
	long     L;
	for (L = 1; L <= LEH; L++){
		long    LE     =  L + N;
		double  S      =  H[L-1] * DE[LE-1];

		CM0            =  CM0 + S;
		*CM1   =  *CM1 + S * E[LE-1];
	}

	double  TT       =  (1.0 - H0) / (CM0 - H0);
	*CM1     *=  TT;
	double  R        = (*CM1) * (*CM1);
	*CM2     =  R * H0;

	for (L = 1; L <= LEH; L++){
		long   LE        =    L + N;
		double  EC       =    E[LE-1] - *CM1;
		double  E2       =    EC * EC;
		H[L-1]  =    H[L-1] * TT;
		double  S        =    H[L-1] * DE[LE-1] * E2;
		*CM2     +=   S;
	}

	*X      *= CM0;
	Y          =  *CM1 / Y;
	Z          =  *CM2 / Z;
}

void AT_Kellerer_interpolation(const long N2,
		const long LEF,
		const long array_size,
		double	F[],
		double	A[],
		double	BI[]){

	A[0]             =  F[1] - F[0];
	BI[0]            =  0.0;
	assert(LEF < array_size);
	F[LEF]  		   =  0.0;
	long   i;
	assert(N2 < array_size);
	for (i = 1; i <= N2; i++){
		long L           =  LEF + i;
		A[L-1]  =  0.0;
		BI[L-1] =  0.0;
	}

	long   L;
	for (L = 2; L <= LEF; L++){
		A[L -1]    =  0.5 * (F[L] - F[L - 1 -1]);
		BI[L -1]   =  A[L-1] + F[L - 1 -1] - F[L -1];
	}
}


void AT_Kellerer_reset(long* N2,
		double* U,
		long* array_size,
		long* LEF,
		long* MIE,
		long* MIF,
		const double E0,
		double E[],
		double DE[],
		double F[],
		double A[],
		double BI[],
		double DI[]){

	if (*N2 <= 256){
		if(*LEF <= 64){
			double S              =  log(2.0);
			double TT             =  (double)(*N2);
			//      N2            =  N2 * 2;
			*N2          +=  (long)(0.1 + exp((double)((long)(log(TT) / S - 0.99)) * S));
			TT                   =  TT / (double)(*N2);
			*U           =  S / (double)(*N2);

			//      theKList            =  AT_SC_INTERP(theKList);
			AT_Kellerer_interpolation( *N2,
					*LEF,
					*array_size,
					F,
					A,
					BI);

			F[*LEF]  =  0;
			long N              =  *MIF;
			*MIF          =  (long)((double)(*MIF) / TT) + 1;    /////////////////////
			*LEF          =  (long)((double)(*LEF) / TT) - 1;    // added (SG) : -1 //
			/////////////////////
			long   K;
			for (K = 1; K <= *LEF; K++){
				long   L              =  *LEF - K + 1;
				double  FLF           =  (double)(L + (*MIF)) * TT - (double)N;
				long   LFF            =  (long)(FLF + 0.5);
				double S              =  FLF - (double)LFF;

				////////////////////////////////////////////////////////////////////////////////
				// Replaced Kellerer's original quadratic interpolation by a slower           //
				// but more correct approach. The original produced wrong interpolation &     //
				// negative values of F on the new N2's grid at very steep, irregular parts   //
				// of F  and eventually systematic deviation of moments.                       //
				////////////////////////////////////////////////////////////////////////////////

				// original:  F[L -1]        =  F[LFF -1] + S * (A[LFF - 1] + S * BI[LFF - 1]);

				F[L -1]        =  F[LFF -1];
				if((S < 0 ) && (LFF >= 2)){
					F[L -1]      =  pow(F[LFF - 1 -1], -1.0 * S) * pow(F[LFF -1], 1.0 + S);
				}
				if((S > 0 ) && (LFF <= *LEF - 1)){
					F[L -1]      =  pow(F[LFF -1], 1.0 - S) * pow(F[LFF], S);
				}
			}

			long   L;
			for (L = *N2; L <= *array_size; L++){;
			double S           =  (double)(L - (*N2)) * (*U);
			double tmp         =  -1.0 * log(1.0 - 0.5 * exp(-S)) / (*U);
			DI[L -1]  =  tmp - (double)(*N2);    // type casts necessary to prevent round of errors (that will eventually yield negative H-values in AT_SC_FOLD
			}

			*MIE          =  *MIF;

			long   J;
			for (J = 1; J <= *array_size; J++){
				double S            =  (double)(J + (*MIE));
				E[J -1]    =  exp(S * (*U)) * E0;
				///////////////////////////////////////////////////////////////////////////
				// addition SG: not to use Kellerer's formula for new DE's, as it is     //
				// not exact (but deviation are small)                                   //
				///////////////////////////////////////////////////////////////////////////
				double* high_E      =  (double*)calloc(*array_size, sizeof(double));
				S                   =  (double)(J + *MIE + 1);
				high_E[J - 1]       =  exp(S * (*U)) * E0;
				DE[J -1]   =  high_E[J -1] - E[J -1];
				free(high_E);
			}
		}else{
			return;
		}
	}else{
		*MIE          =  *MIF;

		long   J;
		for (J = 1; J <= *array_size; J++){
			double S             =  (double)(J + *MIE);
			E[J -1]     =  exp(S * (*U)) * E0;
			///////////////////////////////////////////////////////////////////////////
			// addition SG: not to use Kellerer's formula for new DE's, as it is     //
			// not exact (but deviation are small)                                   //
			///////////////////////////////////////////////////////////////////////////
			double* high_E       =  (double*)calloc(*array_size, sizeof(double));
			S                    =  (double)(J + *MIE + 1);
			high_E[J - 1]        =  exp(S * (*U)) * E0;
			DE[J -1]    =  high_E[J -1] - E[J -1];
			free(high_E);
		}
	}
}

void AT_Kellerer_zero(const long MIF,
		const long array_size,
		const long MIE,
		const long LEF,
		const double F0,
		const double F[],
		const double DE[],
		long* MIH,
		long* LEH,
		double* X,
		double H[]){

	*X          =  0;
	long N              =  *MIH - MIE;

	long   L;
	for (L = 1; L <= *LEH; L++){
		long K            =  L + N;
		*X        +=  H[L -1] * DE[K -1];
	}

	double S            =  (1.0 - F0) * (1.0 - F0) / (*X);
	*X          =  2.0 / S;

	for (L = 1; L <= *LEH; L++){;
	H[L -1]    =  H[L -1] * S;
	}

	N                   =  *MIH - MIF;
	*MIH        =  MIF;
	*LEH        += N;

	long   LL;
	for (LL = 1; LL <= *LEH; LL++){
		long L            =  *LEH + 1 - LL;
		long K            =  L + N;
		H[K -1]  =  H[L -1];
	}

	for (L = 1; L <= N; L++){
		H[L -1]  =  0.0;
	}

	S                   =  F0 * 2.0;

	for (L = 1; L <= LEF; L++){
		H[L -1]  =  H[L -1] + F[L -1] * S;
	}
}


void AT_Kellerer_shrink(const long array_size,
		const long MIE,
		const double shrink_tails_under,
		const double DE[],
		long* MIH,
		long* LEH,
		double H[]){

	double  EX          =  shrink_tails_under;
	double  S           =  0.0;
	long  N             =  *MIH - MIE;

	long   L;
	for (L = 1; L <= *LEH; L++){
		long K            =  L + N;
		S                 =  S + H[L -1] * DE[K -1];
		if(S > 1000.0 * EX){
			*MIH    =  *MIH + L - 1;
			break;}
	}

	long    M           =  L - 1;
	S                   =  0;

	long   K;
	for (K = 1; K <= *LEH; K++){
		L                 =  *LEH + 1 - K;
		long KK           =  L + N;
		S                 =  S + H[L - 1] * DE[KK - 1];
		if(S > EX){
			break;
		}
	}

	*LEH        =  L - M;
	for (L = 1; L <= *LEH; L++){
		K                 =  L + M;
		H[L -1]  =  H[K -1];
	}

	K                   =  *LEH + 1;
	long  KK            =  *LEH + M;
	for (L = K; L <= KK; L++){
		H[L -1]  =  0;
	}
}


void AT_Kellerer_folding(		const long N2,
		const long array_size,
		const long LEF,
		const long MIE,
		const long MIF,
		const double DE[],
		const double DI[],
		long* MIH,
		long* LEH,
		const double F0,
		double* H0,
		double F[],
		double H[],
		double A[],
		double BI[]){
	double*  FDE        =  (double*)calloc(array_size, sizeof(double));

	*H0         =  F0 * F0;
	*MIH        =  MIF + N2;
	*LEH        =  LEF;
	long  K             =  LEF + 1;
	long KK             =  K + N2;

	long   L;
	for (L = K; L <= KK; L++){
		F[L -1]  =  0;
	}

	//theKList            =  AT_SC_INTERP(theKList);
	AT_Kellerer_interpolation( N2,
			LEF,
			array_size,
			F,
			A,
			BI);
	long N              =  MIF - MIE;

	for (L = 1; L <= *LEH; L++){
		K                 =  L + N;
		FDE[L -1]         =  F[L -1] * DE[K -1];
	}

	long   LH;
	for (LH = 1; LH <= *LEH; LH++){
		double   HLH      =  0;
		long   LL         =  LH + N2;
		long   LF;
		for (LF = 1; LF <= LH; LF++){
			K               =  LL - LF;
			double FLF      =  (double)LH - DI[K -1];
			long LFF        =  (long)(FLF + 0.5);
			double S        =  FLF - (double)LFF;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Modification SG: if Kellerer's quadratic interpolation fails, use simple estimate
			double tmp      =  F[LFF -1] + S * (A[LFF -1] + S * BI[LFF -1]);
			if (tmp <0){
				tmp = 0.0;        // Very crude - better to replace by interpolation as done in RESET
			}                   // which is time-consuming, however.
			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			HLH             =  HLH + FDE[LF -1] * tmp;
		}
		H[LH -1] =  HLH - FDE[LH -1] * F[LH -1] * 0.5;
	}

	free(FDE);
}

void   AT_SuccessiveConvolutions( const double  u,
		const long    n_bins_f,
		long*         N2,
		long*         n_bins_f_used,
		double        f_d_Gy[],
		double        f_dd_Gy[],
		double        f[],
		double*       f0,
		double        fdd[],
		double        dfdd[],
		double*       d,
		const bool    write_output,
		const bool    shrink_tails,
		const double  shrink_tails_under,
		const bool    adjust_N2)
{
	long     array_size;                  /**< size of function arrays F..BI */
	double   U;                           /**< Logarithmic stepwith corresponding to N2, U = ln(2)/N2 */

	double   X;                           /**< seem to not be used, TODO remove */
	double   FINAL;                       /**< Final mean impact number * @<d@> from single impact distribution */

	double   CN;                          /**< Number of contributing tracks for linearized f_start, corresp. to u_start */
	long     N1;                          /**< Current convolution number */

	double   CM1;                         /**< Central moment 1 of current distribution H, computed in AT_SC_NORMAL */
	double   CM2;                         /**< Central moment 2 of current distribution H, computed in AT_SC_NORMAL */

	double   D1;                          /**< Central moment 1 of single impact (f1) distribution */

	double   F0;                          /**< Zero bin value for source distribution F (= distribution to be convoluted) */
	double*  F;                           /**< Array holding values of source distribution F */
	long     MIF;                         /**< Index of first bin currently used for F */
	long     LEF;                         /**< Number of bins currently used for F */

	double   H0;                          /**< Zero bin value for resulting distribution H */
	double*  H;                           /**< Array holding values of resulting distribution H */
	long     MIH;                         /**< Index of first bin currently used for H */
	long     LEH;                         /**< Number of bins currently used for H */

	double   E0;                          /**< Value for the LEFT limit of the first energy bin */
	double*  E;                           /**< Array holding energy bin limits (midpoint) for distributions F, H */
	double*  DE;                          /**< Array holding energy bin widths for distributions F, H */
	long     MIE;                         /**< Index of first bin currently used for E */

	double*  DI;                          /**< Auxiliary array that enables easy index operations */
	double*  A;                           /**< Auxiliary array used for Kellerer's quadratic interpolation during AT_SC_FOLD */
	double*  BI;                          /**< Auxiliary array used for Kellerer's quadratic interpolation during AT_SC_FOLD */

	array_size  = n_bins_f;
	U           = log(2.0) / (double)(*N2);

	F        = (double*)calloc(array_size, sizeof(double));
	H        = (double*)calloc(array_size, sizeof(double));
	E        = (double*)calloc(array_size, sizeof(double));
	DE       = (double*)calloc(array_size, sizeof(double));
	DI       = (double*)calloc(array_size, sizeof(double));
	A        = (double*)calloc(array_size, sizeof(double));
	BI       = (double*)calloc(array_size, sizeof(double));

	// Some other initializations
	MIH      = 0;
	MIE      = 0;
	N1       = 0;
	H0       = 0;
	X        = 1;
	CN       = 1;
	CM1      = 1;
	CM2      = 1;

	// Added by Leszek
	MIF      = 0;
	LEF      = 0;

	// Copy input data
	E0      = f_d_Gy[0] * exp(-1.0 * U);

	long   L;
	for (L = 1; L <= array_size; L++){
		E[L -1]      = f_d_Gy[L -1];
		DE[L -1]     = f_dd_Gy[L -1];
		H[L -1]      = f[L -1];
	}

	LEH            = *n_bins_f_used;

	///////////////////////////////////////
	// Fill array for auxilary function that enables easy index operations
	for  (L = *N2; L <= array_size; L++){
		double S          =  (double)(L - *N2) * U;
		double tmp        =  -1.0 * log(1.0 - 0.5 * exp(-S)) / U;
		DI[L -1]    =  tmp - (double)(*N2);
	}    // type casts necessary to prevent round of errors (that will eventually yield negative H-values in AT_SC_FOLD

	///////////////////////////////////////
	// Normalize distribution
	///////////////////////////////////////
	AT_Kellerer_normalize(array_size, MIH, MIE, LEH, H0, E, DE, &CM1, &CM2, &X, H);

	///////////////////////////////////////
	// Get moments of single impact f1
	///////////////////////////////////////
	D1    =    CM1;

	///////////////////////////////////////
	// AT_SC_SHRINK distribution
	///////////////////////////////////////
	if(shrink_tails){
		AT_Kellerer_shrink(array_size, MIE, shrink_tails_under, DE, &MIH, &LEH, H);
	}

	///////////////////////////////////////
	// Get approximation for small hit numbers
	///////////////////////////////////////

	FINAL        = u * D1;  // Final mean impact number

	long n_convolutions    = 0;
	CN               =    FINAL  / D1; // i.e. CN = u , TODO remove FINAL
	while(CN > MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT){
		CN             =    0.5 * CN;
		n_convolutions++;
	}

	H0         =    1.0 - CN;

	for (L = 1; L <= LEH; L++){
		H[L -1]  =  H[L -1] * CN;
	}

	///////////////////////////////////////
	// Convolution loop
	///////////////////////////////////////
	long   j;
	for(j = 0; j < n_convolutions; j++){
		N1        =  N1 + 1;
		CN        =  CN * 2.0;

		for (L = 1; L <= LEH; L++){
			F[L -1]      =  H[L -1];
		}

		F0         =  H0;
		LEF        =  LEH;
		MIF        =  MIH;

		if((CN >= 10.0) && (adjust_N2 == true)){
			AT_Kellerer_reset(N2, &U, &array_size, &LEF, &MIE, &MIF, E0, E, DE, F, A, BI, DI);
		}

		AT_Kellerer_folding(*N2,array_size, LEF, MIE, MIF, DE, DI, &MIH, &LEH, F0, &H0, F, H, A, BI);

		if (F0 < 1e-10){
			X        =  2.0;
		}else{
			AT_Kellerer_zero(MIF, array_size, MIE, LEF, F0, F, DE, &MIH, &LEH, &X, H);
		}

		if(shrink_tails){
			AT_Kellerer_shrink(array_size, MIE, shrink_tails_under, DE, &MIH, &LEH, H);
		}
		AT_Kellerer_normalize(array_size, MIH, MIE, LEH, H0, E, DE, &CM1, &CM2, &X, H);
	}


	//////////////////////////////////////////
	// Copy results back to input structure
	// and adjust according to MIH, MIE
	//////////////////////////////////////////

	*d    = 0.0;

	for (L = 1; L <= array_size; L++){
		f_d_Gy[L -1]      =  0.0;
		f_dd_Gy[L -1]     =  0.0;
		f[L -1]           =  0.0;
		fdd[L -1]         =  0.0;
		dfdd[L -1]        =  0.0;
	}

	long  N        = MIH - MIE;
	for (L = 1; L <= LEH; L++){
		long LE          =  L + N;
		f_d_Gy[L -1]     =  E[LE -1];
		f_dd_Gy[L -1]    =  DE[LE -1];
		f[L -1]          =  H[L-1];
		fdd[L -1]        =  f[L -1] * f_dd_Gy[L -1];
		dfdd[L -1]       =  fdd[L -1] * f_d_Gy[L -1];
		*d              +=  dfdd[L -1];
	}

	*n_bins_f_used = LEH;

	*f0            = H0;

	free(F);
	free(H);
	free(E);
	free(DE);
	free(DI);
	free(A);
	free(BI);
}
