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

#define LOW_FLUENCE_APPROX_FOR_MEAN_NUMBER_OF_TRACKS_CONTRIB  0.002
#define DEBUG_INTERVALS                      8
#define DEBUG_MEAN                           10.0
#define DEBUG_SIGMA                          1.0


long  AT_n_bins_for_single_impact_local_dose_distrib(
		const long    n,
		const double  E_MeV_u[],
		const long    particle_no[],
		const long    material_no,
		const long    rdd_model,
		const double  rdd_parameter[],
		const long    er_model,
		const long    N2,
		const long    stopping_power_source_no)
{
	/* get lowest and highest dose */
	double d_max_Gy    =  0.0;
	double d_min_Gy    =  0.0;

	// TODO think if d_min calculations can be done in smarter way. LET is only needed for Geiss RDD

	long  i;
	for (i = 0; i < n; i++){

		double max_electron_range_m = AT_max_electron_range_m( E_MeV_u[i], (int)material_no, (int)er_model);
		double LET_MeV_cm2_g        = AT_Stopping_Power_MeV_cm2_g_single( stopping_power_source_no, E_MeV_u[i], particle_no[i], material_no);
		double norm_constant_Gy     = AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model);
		double current_d_min_Gy     = AT_RDD_d_min_Gy( E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model, norm_constant_Gy);
		double current_d_max_Gy     = AT_RDD_d_max_Gy( E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model, stopping_power_source_no);

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
		// OLD:
		//		double tmp        =  log10(d_max_Gy/d_min_Gy) / log10(2.0) * ((double)N2);
		//		n_bins_for_singe_impact_local_dose_ditrib        =  (long)(floor(tmp) + 1.0);
		AT_histo_n_bins(     d_min_Gy,
				d_max_Gy,
				AT_N2_to_step(N2),
				AT_histo_log,
				&n_bins_for_singe_impact_local_dose_ditrib);
	} else {
#ifndef NDEBUG
		printf("AT_n_bins_for_singe_impact_local_dose_ditrib: problem in evaluating n_bins_for_singe_impact_local_dose_ditrib: d_min = %g [Gy], d_max = %g [Gy] \n", d_min_Gy, d_max_Gy);
		exit(EXIT_FAILURE);
#endif
	}
	return n_bins_for_singe_impact_local_dose_ditrib + 1;
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
		const long    stopping_power_source_no,
		double        f1_d_Gy[],
		double        f1_dd_Gy[],
		double        frequency_1_Gy_f1[])
{
	long i, j;

	/*
	 * Get relative fluence for beam components
	 * Convert dose to fluence if necessary
	 */
	double*  fluence_cm2    =  (double*)calloc(n, sizeof(double));
	if(fluence_cm2_or_dose_Gy[0] < 0){
		double*  dose_Gy        =  (double*)calloc(n, sizeof(double));
		for (i = 0; i < n; i++){
			dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
		}
		AT_fluence_cm2_from_dose_Gy(  n,
				E_MeV_u,
				particle_no,
				dose_Gy,
				material_no,
				stopping_power_source_no,
				fluence_cm2);
		free( dose_Gy );
	}else{
		for (i = 0; i < n; i++){
			fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
		}
	}
	double*  norm_fluence                                 =  (double*)calloc(n, sizeof(double));
	AT_normalize(    n,
			fluence_cm2,
			norm_fluence);
	free( fluence_cm2 );

	/*
	 * Prepare single impact local dose distribution histogram
	 */

	if(n_bins_f1 > 0){
		const double step		= AT_N2_to_step(N2);
		const long   histo_type	= AT_histo_log;

		// Find lowest and highest dose (looking at ALL particles)
		// TODO: redundant, already used in finding number of bins, replace
		double  d_min_f1      =  f1_parameters[0*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3];
		double  d_max_f1      =  f1_parameters[0*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4];
		for (i = 1; i < n; i++){
			d_min_f1          =  GSL_MIN(f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3], d_min_f1);
			d_max_f1          =  GSL_MAX(f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4], d_max_f1);
		}

		double	 lowest_left_limit_f1 = d_min_f1;

		AT_histo_midpoints(n_bins_f1,
				lowest_left_limit_f1,
				step,
				histo_type,
				f1_d_Gy);

		AT_histo_bin_widths(n_bins_f1,
				lowest_left_limit_f1,
				step,
				histo_type,
				f1_dd_Gy);

		for (i = 0; i < n_bins_f1; i++){
			frequency_1_Gy_f1[i] = 0.0;
		}

		/*
		 * Fill histogram with single impact distribution(s) from individual components
		 */

		// loop over all components (i.e. particles and energies), compute contribution to f1
		long n_bins_used = 1;
		for (i = 0; i < n; i++){

			// Find lowest and highest dose for component
			double  d_min_f1_comp   =  f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3];
			double  d_max_f1_comp   =  f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4];

			// Find position and number of bins for component f1 in overall f1
			long lowest_bin_no_comp  		 	= AT_histo_bin_no(n_bins_f1,
					lowest_left_limit_f1,
					step,
					histo_type,
					d_min_f1_comp);
			long highest_bin_no_comp 			= AT_histo_bin_no(n_bins_f1,
					lowest_left_limit_f1,
					step,
					histo_type,
					d_max_f1_comp);
			long n_bins_f1_comp      			=  highest_bin_no_comp - lowest_bin_no_comp + 1;


			if (n_bins_f1_comp > 1){
				/*
				 * Compute component F1 (accumulated single impact density)
				 * Computation is done with bin limits as sampling points and later differential
				 * f1 will be computed (therefore we need one bin more)
				 * The lowest and highest value for F1 have however to be adjusted as they might
				 * not coincide with the actual min/max values for dose, r, and F1 resp.
				 * They bin widths have to be the same though to assure the integral to be 1
				 *
				 * At the limits F1 will be set to 0 and 1, resp. This enable to account for all
				 * dose, e.g. also in the core, where many radii have the same dose. This procedure,
				 * however, will only work with monotonously falling RDDs which we can assume for all
				 * realistic cases.
				 */
				double*  dose_left_limits_Gy_F1_comp =  (double*)calloc(n_bins_f1_comp + 1, sizeof(double));
				double*  r_m_comp           		 =  (double*)calloc(n_bins_f1_comp + 1, sizeof(double));
				double*  F1_comp        			 =  (double*)calloc(n_bins_f1_comp + 1, sizeof(double));

				// left limit of lowest bin for component
				double   lowest_left_limit_f1_comp = 0.0;
				AT_histo_left_limit(n_bins_f1,
						lowest_left_limit_f1,
						step,
						histo_type,
						lowest_bin_no_comp,
						&lowest_left_limit_f1_comp);

				// get all left limits
				AT_histo_left_limits(n_bins_f1_comp + 1,
						lowest_left_limit_f1_comp,
						step,
						histo_type,
						dose_left_limits_Gy_F1_comp);

				// compute radius as function of dose (inverse RDD),
				// but not for lowest and highest value (i.e. 'n_bins_f1_comp - 1'
				// instead of 'n_bins_f1_comp + 1' and &dose_left_limits_Gy_F1_comp[1]
				// as entry point instead of dose_left_limits_Gy_F1_comp
				// exit in case of problems
				int inverse_RDD_status_code = AT_r_RDD_m  (  n_bins_f1_comp - 1,
						&dose_left_limits_Gy_F1_comp[1],
						E_MeV_u[i],
						particle_no[i],
						material_no,
						rdd_model,
						rdd_parameter,
						er_model,
						stopping_power_source_no,
						&r_m_comp[1]);

				if( inverse_RDD_status_code != 0 ){
#ifndef NDEBUG
					printf("Problem in evaluating inverse RDD in AT_SC_get_f1, probably wrong combination of ER and RDD used\n");
#endif
					char rdd_model_name[100];
					AT_RDD_name_from_number(rdd_model, rdd_model_name);
					char er_model_name[100];
					getERName( er_model, er_model_name);
#ifndef NDEBUG
					printf("rdd_model: %ld (%s), er_model: %ld (%s)\n", rdd_model, rdd_model_name, er_model, er_model_name);
					exit(EXIT_FAILURE);
#endif
				}

				// compute F1 as function of radius
				// use F1 - 1 instead of F1 to avoid numeric cut-off problems
				double r_max_m_comp = f1_parameters[i * AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 2];
				for (j = 1; j < n_bins_f1_comp; j++){
					F1_comp[j]            = gsl_pow_2(r_m_comp[j] / r_max_m_comp);
				}

				// Set extreme values of F1
				F1_comp[0]					= 1.0;
				F1_comp[n_bins_f1_comp]		= 0.0;

				FILE* output = fopen("F_output.csv", "w");
				fprintf(output, "bin.no;r.m;d.Gy;F1\n");
				for (j = 0; j < n_bins_f1_comp + 1; j++){
					fprintf(output,
							"%ld;%7.6e;%7.6e;%7.6e\n",
							j, r_m_comp[j], dose_left_limits_Gy_F1_comp[j], F1_comp[j]);
				}
				fclose(output);

				// now compute f1 as the derivative of F1 and add to overall f1
				double f1_comp;
				for (j = 0; j < n_bins_f1_comp; j++){
					f1_comp				  						=  (F1_comp[j] - F1_comp[j + 1]) / (dose_left_limits_Gy_F1_comp[j + 1] - dose_left_limits_Gy_F1_comp[j]);
					frequency_1_Gy_f1[lowest_bin_no_comp + j]   += norm_fluence[i] * f1_comp;
				}

				// adjust the density in first and last bin, because upper limit is not d.max.Gy and lower not d.min.Gy
				free(dose_left_limits_Gy_F1_comp);
				free(r_m_comp);
				free(F1_comp);
			}
			else{ // in case of n_bins_df == 1 (all doses fall into single bin, just add a value of 1.0
				frequency_1_Gy_f1[lowest_bin_no_comp ]        +=  norm_fluence[i] * 1.0 / f1_dd_Gy[lowest_bin_no_comp];
			}

			// remember highest bin used
			n_bins_used          =  GSL_MAX(n_bins_used, highest_bin_no_comp);
		}

		// normalize f1 (should be ok anyway but there could be small round-off errors)
		double  f1_norm    =  0.0;
		for (i = 0; i < n_bins_f1; i++){
			f1_norm    +=    frequency_1_Gy_f1[i] * f1_dd_Gy[i];
		}
		for (i = 0; i < n_bins_f1; i++){
			frequency_1_Gy_f1[i]    /=    f1_norm;
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
	while(*u_start  > LOW_FLUENCE_APPROX_FOR_MEAN_NUMBER_OF_TRACKS_CONTRIB){
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

void AT_Kellerer_normalize(const long values_first_bin,
		const double value_midpoints[],
		const double value_bin_widths[],
		const long frequency_n_bins,
		const double frequency_zero_bin,
		const long frequency_first_bin,
		double frequency[]){

	long    i, j;
	/* compute sum */
	long    bin_shift	=  frequency_first_bin - values_first_bin;
	double  sum			=  frequency_zero_bin;

	//	TODO: why not: for (i = first_bin_frequency; i < first_bin_frequency + n_bins_frequency; i++){
	for (i = 0; i < frequency_n_bins; i++){
		j		    =  i + bin_shift;
		sum			=  sum + frequency[i] * value_bin_widths[j];
	}

	/* normalizing factor */
	double  norm_factor     =  (1.0 - frequency_zero_bin) / (sum - frequency_zero_bin);

	/* normalize values */
	for (i = 0; i < frequency_n_bins; i++){
		frequency[i]  	*=    norm_factor;
	}
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
		const long array_size,
		long* LEF,
		long* MIE,
		long* MIF,
		const double E0,
		double E[],
		double DE[],
		double F[],
		double DI[]){

	double U           = log(2.0) / (double)(*N2);
	if (*N2 <= 256){
		if(*LEF <= 64){
			double S              =  log(2.0);
			double TT             =  (double)(*N2);
			//      N2            =  N2 * 2;
			*N2          +=  (long)(0.1 + exp((double)((long)(log(TT) / S - 0.99)) * S));
			TT                   =  TT / (double)(*N2);
			U           =  S / (double)(*N2);

			/* Compute coefficients for quadratic interpolation */
			double*	A	= (double*)calloc(array_size, sizeof(double));
			double*	BI	= (double*)calloc(array_size, sizeof(double));
			AT_Kellerer_interpolation( *N2,
					*LEF,
					array_size,
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
			for (L = *N2; L <= array_size; L++){;
			double S           =  (double)(L - (*N2)) * U;
			double tmp         =  -1.0 * log(1.0 - 0.5 * exp(-S)) / U;
			DI[L -1]  =  tmp - (double)(*N2);    // type casts necessary to prevent round of errors (that will eventually yield negative H-values in AT_SC_FOLD
			}

			*MIE          =  *MIF;

			long   J;
			for (J = 1; J <= array_size; J++){
				double S            =  (double)(J + (*MIE));
				E[J -1]    =  exp(S * U) * E0;
				///////////////////////////////////////////////////////////////////////////
				// addition SG: not to use Kellerer's formula for new DE's, as it is     //
				// not exact (but deviation are small)                                   //
				///////////////////////////////////////////////////////////////////////////
				double* high_E      =  (double*)calloc(array_size, sizeof(double));
				S                   =  (double)(J + *MIE + 1);
				high_E[J - 1]       =  exp(S * U) * E0;
				DE[J -1]   =  high_E[J -1] - E[J -1];
				free(high_E);
			}
			free(A);
			free(BI);
		}else{
			return;
		}
	}else{
		*MIE          =  *MIF;

		long   J;
		for (J = 1; J <= array_size; J++){
			double S             =  (double)(J + *MIE);
			E[J -1]     =  exp(S * U) * E0;
			///////////////////////////////////////////////////////////////////////////
			// addition SG: not to use Kellerer's formula for new DE's, as it is     //
			// not exact (but deviation are small)                                   //
			///////////////////////////////////////////////////////////////////////////
			double* high_E       =  (double*)calloc(array_size, sizeof(double));
			S                    =  (double)(J + *MIE + 1);
			high_E[J - 1]        =  exp(S * U) * E0;
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
		double H[]){

	double X          =  0;
	long N              =  *MIH - MIE;

	long   L;
	for (L = 1; L <= *LEH; L++){
		long K            =  L + N;
		X        +=  H[L -1] * DE[K -1];
	}

	double S            =  (1.0 - F0) * (1.0 - F0) / X;
	X          =  2.0 / S;

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


void AT_Kellerer_folding(		const long n_bins,
		const long bins_per_factor_2,
		const double delta_i[],
		const long values_first_bin,
		const double values_bin_widths[],
		const long frequency_n_bins_last,
		const long frequency_first_bin_last,
		const double frequency_zero_bin_last,
		double frequency_last[],
		long* frequency_n_bins,
		long* frequency_first_bin,
		double* frequency_zero_bin,
		double frequency[]){

	long   i, j, int_k;
	double k, frac_k;

	/* Convolution shifts support by factor of 2 while the length stays the same */
	*frequency_first_bin	=  frequency_first_bin_last + bins_per_factor_2;
	*frequency_n_bins		=  frequency_n_bins_last;

	/* Clean bins above support */
	for (i = frequency_n_bins_last; i < frequency_n_bins_last + bins_per_factor_2; i++){
		frequency_last[i]  =  0;
	}

	/* Precompute coefficients for quadratic interpolation */
	double*	lin_coeff	= (double*)calloc(n_bins, sizeof(double));
	double*	qua_coeff	= (double*)calloc(n_bins, sizeof(double));
	AT_Kellerer_interpolation( bins_per_factor_2,
			frequency_n_bins_last,
			n_bins,
			frequency_last,
			lin_coeff,
			qua_coeff);

	/* Precompute bin contents */
	long bin_shift	=  frequency_first_bin_last - values_first_bin;
	double*  bin_content	=  (double*)calloc(n_bins, sizeof(double));
	for (i = 0; i < *frequency_n_bins; i++){
		j                 =  i + bin_shift;
		bin_content[i]    =  frequency_last[i] * values_bin_widths[j];
	}

	/* Convolution */
	double 	tmp, sum;
	long	i_shift;
	for (i = 0; i < *frequency_n_bins; i++){
		sum      =  0;
		i_shift  =  i + bins_per_factor_2;
		for (j = 0; j <= i; j++){
			k		=  (double)i - delta_i[i_shift - j - 1];
			int_k	=  (long)(k + 0.5);
			frac_k	=  k - (double)int_k;
			tmp		=  frequency_last[int_k] + frac_k * (lin_coeff[int_k] + frac_k * qua_coeff[int_k]);
			/* If quadratic interpolation fails use simple correction */
			if (tmp <0){ tmp = 0.0;}        // TODO: Very crude - better to replace by interpolation as done in RESET which is time-consuming, however.
			sum     += bin_content[j] * tmp;
		}
		frequency[i] =  sum - bin_content[i] * frequency_last[i] * 0.5;
	}

	free(bin_content);
	free(lin_coeff);
	free(qua_coeff);

	/* Compute new zero bin */
	*frequency_zero_bin		=  frequency_zero_bin_last * frequency_zero_bin_last;
}

void   AT_SuccessiveConvolutions( const double  final_mean_number_of_tracks_contrib,
		const long    n_bins, // TODO: allocate within SC routine
		long*         bins_per_factor_2,
		long*         n_bins_single_impact_local_dose_distrib,
		double        single_impact_local_dose_Gy[],
		double        single_impact_local_dose_bin_width_Gy[],
		double        single_impact_frequency[],
		double*       single_impact_frequency_zero_bin,
		double        fdd[],
		double        dfdd[],
		double*       d,
		const bool    write_output,
		const bool    shrink_tails,
		const double  shrink_tails_under,
		const bool    adjust_dose_spacing)
{
	long	i;

	/* Get number of convolutions */
	long n_convolutions    = 0;
	double current_mean_number_of_tracks_contrib	= final_mean_number_of_tracks_contrib;
	while(current_mean_number_of_tracks_contrib > LOW_FLUENCE_APPROX_FOR_MEAN_NUMBER_OF_TRACKS_CONTRIB){
		current_mean_number_of_tracks_contrib             /=    2;
		n_convolutions++;
	}

	/* Copy input data */
	long 	frequency_first_bin      		= 0;
	long	frequency_n_bins      			= *n_bins_single_impact_local_dose_distrib;
	double*	frequency        				= (double*)calloc(n_bins, sizeof(double));
	double	frequency_zero_bin       		= 0.0f;
	memcpy(frequency, single_impact_frequency, frequency_n_bins * sizeof(double));

	long	local_dose_first_bin			= 0;
	long	local_dose_n_bins		   		= n_bins;
	double	local_dose_lowest_left_limit    = single_impact_local_dose_Gy[0] * exp(-1.0 * log(2.0) / (double)(*bins_per_factor_2));		// lowest midpoint - half bin width
	double*	local_dose_midpoints        	= (double*)calloc(n_bins, sizeof(double));
	double* local_dose_bin_widths       	= (double*)calloc(n_bins, sizeof(double));
	memcpy(local_dose_midpoints,  single_impact_local_dose_Gy,  local_dose_n_bins * sizeof(double));
	memcpy(local_dose_bin_widths, single_impact_local_dose_bin_width_Gy, local_dose_n_bins * sizeof(double));

	/* Normalize distribution */
	AT_Kellerer_normalize(	local_dose_first_bin,
			local_dose_midpoints,
			local_dose_bin_widths,
			frequency_n_bins,
			frequency_zero_bin,
			frequency_first_bin,
			frequency);

	/* Cut tails of distribution */
	if(shrink_tails){
		AT_Kellerer_shrink(	n_bins,
				local_dose_first_bin,
				shrink_tails_under,
				local_dose_bin_widths,
				&frequency_first_bin,
				&frequency_n_bins, frequency);
	}

	/* Small number approximation */
	frequency_zero_bin         =    1.0 - current_mean_number_of_tracks_contrib;
	for (i = 0; i < frequency_n_bins; i++){
		frequency[i]  *=  current_mean_number_of_tracks_contrib;
	}

	/* Precompute non-integer differences in dose indices */
	double*	delta_i		= (double*)calloc(n_bins, sizeof(double));
	double  tmp;
	for  (i = *bins_per_factor_2; i <= n_bins; i++){
		tmp		          =  (double)(i - *bins_per_factor_2) * log(2.0) / (double)(*bins_per_factor_2);
		tmp			      =  -1.0 * log(1.0 - 0.5 * exp(-tmp)) / (log(2.0) / (double)(*bins_per_factor_2));
		delta_i[i-1]      =  tmp - (double)(*bins_per_factor_2);
	}

	/* Actual convolution loop */
	double* frequency_last        			= (double*)calloc(n_bins, sizeof(double));
	double  frequency_zero_bin_last			= 0.0;
	long	frequency_first_bin_last      	= 0;
	long	frequency_n_bins_last      		= 0;

	for(i = 0; i < n_convolutions; i++){
		current_mean_number_of_tracks_contrib	*= 2.0;

		/* Copy last distribution */
		frequency_zero_bin_last         =  frequency_zero_bin;
		frequency_n_bins_last        	=  frequency_n_bins;
		frequency_first_bin_last        =  frequency_first_bin;
		memcpy(frequency_last, frequency, frequency_n_bins * sizeof(double));

		if((current_mean_number_of_tracks_contrib >= 10.0) && (adjust_dose_spacing == true)){
			AT_Kellerer_reset(	bins_per_factor_2,
					n_bins,
					&frequency_n_bins_last,
					&local_dose_first_bin,
					&frequency_first_bin_last,
					local_dose_lowest_left_limit,
					local_dose_midpoints,
					local_dose_bin_widths,
					frequency_last,
					delta_i);
		}

		AT_Kellerer_folding(	n_bins,
				*bins_per_factor_2,
				delta_i,
				local_dose_first_bin,
				local_dose_bin_widths,
				frequency_n_bins_last,
				frequency_first_bin_last,
				frequency_zero_bin_last,
				frequency_last,
				&frequency_n_bins,
				&frequency_first_bin,
				&frequency_zero_bin,
				frequency);

		if (frequency_zero_bin_last >= 1e-10){
			AT_Kellerer_zero(	frequency_first_bin_last,
					n_bins,
					local_dose_first_bin,
					frequency_n_bins_last,
					frequency_zero_bin_last,
					frequency_last,
					local_dose_bin_widths,
					&frequency_first_bin,
					&frequency_n_bins,
					frequency);
		}

		if(shrink_tails){
			AT_Kellerer_shrink(	n_bins,
					local_dose_first_bin,
					shrink_tails_under,
					local_dose_bin_widths,
					&frequency_first_bin,
					&frequency_n_bins,
					frequency);
		}

		AT_Kellerer_normalize(	local_dose_first_bin,
				local_dose_midpoints,
				local_dose_bin_widths,
				frequency_n_bins,
				frequency_zero_bin,
				frequency_first_bin,
				frequency);
	}

	//////////////////////////////////////////
	// Copy results back to input structure
	//////////////////////////////////////////

	*d    = 0.0;

	for (i = 0; i < n_bins; i++){
		single_impact_local_dose_Gy[i]      =  0.0;
		single_impact_local_dose_bin_width_Gy[i]     =  0.0;
		single_impact_frequency[i]           =  0.0;
		fdd[i]         =  0.0;
		dfdd[i]        =  0.0;
	}

	long  N        = frequency_first_bin - local_dose_first_bin;
	for (i = 0; i < frequency_n_bins; i++){
		long j          =  i + N;
		single_impact_local_dose_Gy[i]     =  local_dose_midpoints[j];
		single_impact_local_dose_bin_width_Gy[i]    =  local_dose_bin_widths[j];

		single_impact_frequency[i]          =  frequency[i];
		fdd[i]        =  single_impact_frequency[i] * single_impact_local_dose_bin_width_Gy[i];
		dfdd[i]       =  fdd[i] * single_impact_local_dose_Gy[i];
		*d              +=  dfdd[i];
	}

	*n_bins_single_impact_local_dose_distrib = frequency_n_bins;

	*single_impact_frequency_zero_bin            = frequency_zero_bin;


	free(frequency_last);
	free(frequency);
	free(local_dose_midpoints);
	free(local_dose_bin_widths);
	free(delta_i);
}

