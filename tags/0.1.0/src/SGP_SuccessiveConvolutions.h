#ifndef SGP_SUCCESSIVECONVOLUTIONS_H_
#define SGP_SUCCESSIVECONVOLUTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "SGP_Utils.h"
#include "SGP_GammaResponse.h"

#define MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT		0.002
#define	DEBUG_INTERVALS							8
#define DEBUG_MEAN								10.0f
#define DEBUG_SIGMA								1.0f


void SGP_SC_get_f1_array_size(	long*	n,
								float*	E_MeV_u,
								long*	particle_no,
								char*	material_name,
								float*	parameter,
								long*	N2,
								// from here: return values
								long*	n_bins_f1,
								bool*	debug);

void SGP_SC_get_f1(	long*	n,
					float*	E_MeV_u,
					long*	particle_no,
					float*	fluence_cm2,
					char*	material_name,
					float*	parameter,
					long*	N2,
					long*	n_bins_f1,
					// return values
					float*	norm_fluence,
					float*	LET_MeV_cm2_g,
					float*	r_min_m,
					float*	r_max_m,
					float*	d_min_Gy,
					float*	d_max_Gy,
					float*	k_Gy,
					float*	single_impact_fluence_cm2,
					float*	single_impact_dose_Gy,
					float*	dose_contribution_Gy,
					float*	total_fluence_cm2,
					float*	total_dose_Gy,
					float*	ave_E_MeV,
					float*	dw_E_MeV,
					float*	ave_LET_MeV_cm2_g,
					float*	dw_LET_MeV_cm2_g,
					float*	u,
					float*	f1_d_Gy,
					float*	f1_dd_Gy,
					float*	f1,
					bool*	debug);

void SGP_SC_get_f_array_size(	float*	u,
								float*	fluence_factor,
								long*	N2,
								long*	n_bins_f1,
								float*	f1_d_Gy,
								float*	f1_dd_Gy,
								float*	f1,
								// from here: return values
								long*	n_bins_f,
								float*	u_start,
								long*	n_convolutions);


void	SGP_SC_get_f_start(	float*	u_start,
							long*	n_bins_f1,
							long*	N2,
							float*	f1_d_Gy,
							float*	f1_dd_Gy,
							float*	f1,
							long*	n_bins_f,
							// from here: return values
							float*	f_d_Gy,
							float*	f_dd_Gy,
							float*	f_start);


void SGP_SuccessiveConvolutions(	float*	u,
									long*	n_bins_f,
									long*	N2,
									// input + return values
									long*	n_bins_f_used,
									float*	f_d_Gy,
									float*	f_dd_Gy,
									float*	f,
									// return values
									float*	f0,
									float*	fdd,
									float*	dfdd,
									float*	d,
									bool*	write_output,
									bool*	shrink_tails,
									float*	shrink_tails_under,
									bool*	adjust_N2);


void SGP_SC_Loop(	long*	n,
					float*	E_MeV_u,
					long*	particle_no,
					float*	fluence_cm2,
					long*	slab_no,
					char*	material_name,
					float*	parameter,
					long*	N2,
					long*	n_slabs,
					long*	n_gamma_parameter,
					long*	gamma_model,
					float*	gamma_parameter,
					long*	verbose_level,
					char*	output_fileName,
					// return values
					float*	u,
					float*	total_d_Gy,
					float*	d,
					float*	S_HCP,
					float*	S_gamma,
					float*	efficiency,
					float*	S_HCP_total,
					float*	S_gamma_total,
					float*	efficiency_total);


void SGP_SC_LoopS(	long*	n,
					float*	E_MeV_u,
					long*	particle_no,
					float*	fluence_cm2,
					long*	slab_no,
					char**	material_name,
					float*	parameter,
					long*	N2,
					long*	n_slabs,
					long*	n_gamma_parameter,
					long*	gamma_model,
					float*	gamma_parameter,
					long*	verbose_level,
					char**	output_fileName,
					// return values
					float*	u,
					float*	total_d_Gy,
					float*	d,
					float*	S_HCP,
					float*	S_gamma,
					float*	efficiency,
					float*	S_HCP_total,
					float*	S_gamma_total,
					float*	efficiency_total){
SGP_SC_Loop(	n,
				E_MeV_u,
				particle_no,
				fluence_cm2,
				slab_no,
				*material_name,
				parameter,
				N2,
				n_slabs,
				n_gamma_parameter,
				gamma_model,
				gamma_parameter,
				verbose_level,
				*output_fileName,
				// return values
				u,
				total_d_Gy,
				d,
				S_HCP,
				S_gamma,
				efficiency,
				S_HCP_total,
				S_gamma_total,
				efficiency_total);
};

void	SGP_SC_get_f1_array_size(	long*	n,
									float*	E_MeV_u,
									long*	particle_no,
									char*	material_name,
									float*	parameter,
									long*	N2,
									// from here: return values
									long*	n_bins_f1,
									bool*	debug)
{
	// Allocate memory, also for the return variable not needed here (dummy1 - dummy6)
	float*	dummy1			=	(float*)calloc(*n, sizeof(float));
	float*	dummy2			=	(float*)calloc(*n, sizeof(float));
	float*	dummy3			=	(float*)calloc(*n, sizeof(float));
	float*	dummy4			=	(float*)calloc(*n, sizeof(float));
	float*	dummy5			=	(float*)calloc(*n, sizeof(float));
	float*	dummy6			=	(float*)calloc(*n, sizeof(float));

	float*	d_min_Gy		=	(float*)calloc(*n, sizeof(float));
	float*	d_max_Gy		=	(float*)calloc(*n, sizeof(float));

	// get RDD parameters for all particles and energies
	SGP_RDD_f1_parameters(	n,
							E_MeV_u,
							particle_no,
							material_name,
							parameter,
							// return:
							dummy1,
							dummy2,
							dummy3,
							d_min_Gy,
							d_max_Gy,
							dummy4,
							dummy5,
							dummy6);

	// get lowest and highest dose
	float	d_max			=	d_max_Gy[0];
	float	d_min			=	d_min_Gy[0];

	long 	i;
	for (i = 1; i < *n; i++){
		d_max					=	FMAX(d_max_Gy[i], d_max);
		d_min					=	FMIN(d_min_Gy[i], d_min);}

	// get number of bins needed to span that dose range
	*n_bins_f1				=	(long)floor((float)(log10(d_max/d_min) / log10(2.0f) * ((float)*N2))) + 1;


	// DEBUG //
	if(*debug){
		*n_bins_f1	=	*N2 * DEBUG_INTERVALS;
	}

	return;
}


void	SGP_SC_get_f1(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						float*	fluence_cm2,
						char*	material_name,
						float*	parameter,
						long*	N2,
						long*	n_bins_f1,
						// from here: return values
						float*	norm_fluence,
						float*	LET_MeV_cm2_g,
						float*	r_min_m,
						float*	r_max_m,
						float*	d_min_Gy,
						float*	d_max_Gy,
						float*	k_Gy,
						float*	single_impact_fluence_cm2,
						float*	single_impact_dose_Gy,
						float*	dose_contribution_Gy,
						float*	total_fluence_cm2,
						float*	total_dose_Gy,
						float*	ave_E_MeV,
						float*	dw_E_MeV,
						float*	ave_LET_MeV_cm2_g,
						float*	dw_LET_MeV_cm2_g,
						float*	u,
						float*	f1_d_Gy,
						float*	f1_dd_Gy,
						float*	f1,
						bool*	debug)
{
	// get RDD parameters for all particles and energies
	SGP_RDD_f1_parameters(	n,
								E_MeV_u,
								particle_no,
								material_name,
								parameter,
								// return:
								LET_MeV_cm2_g,
								r_min_m,
								r_max_m,
								d_min_Gy,
								d_max_Gy,
								k_Gy,
								single_impact_fluence_cm2,
								single_impact_dose_Gy);

	// normalize fluence, get total fluence and dose, eff. LET and mean impact parameter u,
	*total_fluence_cm2		= 0.0f;
	*total_dose_Gy			= 0.0f;

	// if fluences < 0 they are supposed to be D.set.Gy, so in that case convert them first
	// only the first entry will be check
	long 	i;
	float*	fluence_cm2_local		=	(float*)calloc(*n, sizeof(float));
	if (fluence_cm2[0] >= 0){
		for (i = 0; i < *n; i++){
			fluence_cm2_local[i]		= fluence_cm2[i];}
	}
	else{
		for (i = 0; i < *n; i++){
			fluence_cm2_local[i]		= -1.0f * fluence_cm2[i] / (LET_MeV_cm2_g[i] * MeV_g_to_J_kg);}
	}

	for (i = 0; i < *n; i++){
		*total_fluence_cm2	+=	fluence_cm2_local[i];
	}


	float u_single;
	*u								=	0.0f;
	*ave_E_MeV						=	0.0f;
	*dw_E_MeV						=	0.0f;
	*ave_LET_MeV_cm2_g				=	0.0f;
	*dw_LET_MeV_cm2_g				=	0.0f;

	for (i = 0; i < *n; i++){
		norm_fluence[i]				=	fluence_cm2_local[i] / *total_fluence_cm2;
		u_single					=	fluence_cm2_local[i] / single_impact_fluence_cm2[i];
		dose_contribution_Gy[i]		=	u_single * single_impact_dose_Gy[i];
		*total_dose_Gy				+=	dose_contribution_Gy[i];
		*ave_E_MeV					+=	norm_fluence[i] * E_MeV_u[i];
		*dw_E_MeV					+=	dose_contribution_Gy[i] * E_MeV_u[i];
		*ave_LET_MeV_cm2_g			+=	norm_fluence[i] * LET_MeV_cm2_g[i];
		*dw_LET_MeV_cm2_g			+=	dose_contribution_Gy[i] * LET_MeV_cm2_g[i];
		*u							+=	norm_fluence[i] * single_impact_dose_Gy[i];
	}

	*dw_E_MeV					/= *total_dose_Gy;
	*dw_LET_MeV_cm2_g			/= *total_dose_Gy;
	*u							= *total_dose_Gy / *u;

	//  create all-over f1-data-frame
	float	d_max			=	d_max_Gy[0];
	float	d_min			=	d_min_Gy[0];

	for (i = 1; i < *n; i++){
		d_max					=	FMAX(d_max_Gy[i], d_max);
		d_min					=	FMIN(d_min_Gy[i], d_min);}

	float	U				=	(float)(log(2.0f) / (float)(*N2));


	float*	d_df_low			=	(float*)calloc(*n_bins_f1, sizeof(float));
	float*	d_df_mid			=	(float*)calloc(*n_bins_f1, sizeof(float));
	float*	d_df_high			=	(float*)calloc(*n_bins_f1, sizeof(float));
	float*	dd_df				=	(float*)calloc(*n_bins_f1, sizeof(float));

	for (i = 0; i < *n_bins_f1; i++){
	// TO DO: check if n.bins sufficient

		d_df_low[i]					= 	d_min * (float)exp((float)i * U);
		d_df_high[i]				= 	d_min * (float)exp(((float)i + 1) * U);
		d_df_mid[i]					=	d_min * (float)exp(((float)i + 0.5f) * U);
		dd_df[i]					=	d_df_high[i] - d_df_low[i];							// OBS: using Kellerer's bin-width = mid.point * U is not entirely correct

		f1[i]						=	0.0f;
	}

	long	n_bins_used				=	1;

	// loop over all particles and energies, compute contribution to f1
	long 	k;
	for (k = 0; k < *n; k++){

		float	d_max_k				=	d_max_Gy[k];
		float	d_min_k				=	d_min_Gy[k];

		// find first and last bin to fit this particle's contribution into the all-over f1-frame
		long	i_low, i_high;
		locate(d_df_low, n_bins_f1, &d_min_k, &i_low);
		locate(d_df_low, n_bins_f1, &d_max_k, &i_high);
		i_low						-=	1;
		i_high						-=	1;

		long	n_bins_df			=	i_high - i_low + 1;  // changed from + 2

		if (n_bins_df > 1){
			float*	d_low				=	(float*)calloc(n_bins_df, sizeof(float));
			float*	d_mid				=	(float*)calloc(n_bins_df, sizeof(float));
			float*	d_high				=	(float*)calloc(n_bins_df, sizeof(float));
			float*	dd					=	(float*)calloc(n_bins_df, sizeof(float));
			float*	r					=	(float*)calloc(n_bins_df, sizeof(float));
			float*	F1_1				=	(float*)calloc(n_bins_df, sizeof(float));
			float*	f1_k				=	(float*)calloc(n_bins_df - 1, sizeof(float));

			// extract the corresponding part from the all-over frame
			for (i = 0; i < n_bins_df; i++){
				d_low[i]					= 	d_df_low[i_low + i];
				d_high[i]					= 	d_df_high[i_low + i];
				d_mid[i]					=	d_df_mid[i_low + i];
				dd[i]						=	d_high[i] - d_low[i];
			}

			// and adjust the edges
			d_low[1-1]					=	d_min_k;
			d_low[n_bins_df-1]			=	d_max_k;

			d_mid[1-1]					=	(float)sqrt(d_low[1 - 1] * d_high[1 - 1]);
			d_mid[n_bins_df-1-1]		=	(float)sqrt(d_low[n_bins_df - 1] * d_high[n_bins_df - 1]);
			d_mid[n_bins_df-1]			=	0;

			d_high[n_bins_df-1-1]		=	d_max_k;
			d_high[n_bins_df-1]			=	0.0f;	//??

			dd[n_bins_df-1]				=	0.0f;


			// now compute r, F1, and f1, this could be any RDD if implemented
			SGP_r_RDD_m	(	&n_bins_df,
								d_low,
								&E_MeV_u[k],
								&particle_no[k],
								material_name,
								parameter,
								r);

			for (i = 0; i < n_bins_df; i++){
				F1_1[i]						= (r[i] / r_max_m[k]) * (r[i] / r_max_m[k]);}				// F1 - 1 instead of F1 to avoid numeric cut-off problems

			F1_1[n_bins_df-1]		=	0.0f;

			// now compute f1 as the derivative of F1
			for (i = 0; i < (n_bins_df - 1); i++){
				f1_k[i]					=	-1.0f * (F1_1[i + 1] - F1_1[i]) / dd[i];}

			// adjust the density in first and last bin, because upper limit is not d.max.Gy and lower not d.min.Gy
			f1_k[1-1]				=	f1_k[1-1] * dd[1-1] / dd_df[i_low];
			f1_k[n_bins_df-1-1]		=	f1_k[n_bins_df-1-1] * dd[n_bins_df-1-1] / dd_df[i_high - 1];

			// and paste f1 for this energy /particle into the over all data frame according to rel. fluence
			for (i = 0; i < (n_bins_df - 1); i++){
				f1[i_low + i]			+=	norm_fluence[k] * f1_k[i];}


			free(d_low);
			free(d_mid);
			free(d_high);
			free(dd);
			free(r);
			free(F1_1);
			free(f1_k);
			}
		else{ // n_bins_df == 1
			f1[i_low ]				+=	norm_fluence[k] * 1.0f / dd_df[i_low];
		}

		// remember highest bin used
		n_bins_used				=	LMAX(n_bins_used, i_high);
	}

	// copy back for the dose axis
	for (i = 0; i < *n_bins_f1; i++){
		f1_d_Gy[i]		=	d_df_mid[i];
		f1_dd_Gy[i]		=	dd_df[i];
	}

	free(d_df_low);
	free(d_df_mid);
	free(d_df_high);
	free(dd_df);

	free(fluence_cm2_local);

	// DEBUG //
	if(*debug){
		float d_min_debug	=	DEBUG_MEAN / (float)pow(2, DEBUG_INTERVALS / 2);
		for (i = 0; i <*n_bins_f1; i++){
			f1_d_Gy[i]		=	d_min_debug * (float)exp(((float)i + 0.5f) * U);
			f1_dd_Gy[i]		=	d_min_debug * (float)(exp(((float)i + 1.0f) * U) - exp(((float)i) * U));
			f1[i]			=	1 / DEBUG_SIGMA * (float)exp(-0.5f * pow((f1_d_Gy[i] - DEBUG_MEAN)/DEBUG_SIGMA, 2));
		}
	}
	// DEBUG //


	// normalize f1 (should be ok anyway but there could be small round-off errors)
	float	f1_norm		=	0.0f;
	for (i = 0; i < *n_bins_f1; i++){
		f1_norm		+=		f1[i] * f1_dd_Gy[i];}
	for (i = 0; i < *n_bins_f1; i++){
		f1[i]		/=		f1_norm;}


	return;
}


void	SGP_SC_get_f_array_size(	float*	u,
									float*	fluence_factor,
									long*	N2,
									long*	n_bins_f1,
									float*	f1_d_Gy,
									float*	f1_dd_Gy,
									float*	f1,
									// from here: return values
									long*	n_bins_f,
									float*	u_start,
									long*	n_convolutions)
{
	// Get expectation value of dose from f1
	float	d_f1_Gy		=	0.0f;

	long 	i;
	for (i = 0; i < *n_bins_f1; i++){
		d_f1_Gy		+=	f1_d_Gy[i] * f1_dd_Gy[i] * f1[i];
	}

	// The dose set by the input data is therefore
	float	d_f_Gy		= (*u) * (*fluence_factor) * d_f1_Gy;

	// How many convolution are necessary starting from a small mean
	// impact number that allows linear approximation (e.g. 0.002)
	*n_convolutions		= 0;

	*u_start			=		d_f_Gy	/ d_f1_Gy;
	while(*u_start	> MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT){
		*u_start			=		0.5f * (*u_start);
		(*n_convolutions)++;
	}

	// Every convolution will add a factor of two, so the array for f has to
	// be expanded from f1 by N2 * n_convolutions
	*n_bins_f			=	(*n_convolutions + 1) * (*N2);
	*n_bins_f			+=	*n_bins_f1;
	return;
}


void	SGP_SC_get_f_start(			float*	u_start,
									long*	n_bins_f1,
									long*	N2,
									float*	f1_d_Gy,
									float*	f1_dd_Gy,
									float*	f1,
									long*	n_bins_f,
									// from here: return values
									float*	f_d_Gy,
									float*	f_dd_Gy,
									float*	f_start)
{
	// temporary arrays
	float*	d_low				=	(float*)calloc(*n_bins_f, sizeof(float));
	float*	d_high				=	(float*)calloc(*n_bins_f, sizeof(float));

	float	U					=	(float)(log(2.0f) / (float)(*N2));

	long 	i;
	for (i = 0; i < *n_bins_f; i++){
		d_low[i]					= 	f1_d_Gy[0] * (float)exp(((float)i - 0.5f)* U);
		d_high[i]					= 	f1_d_Gy[0] * (float)exp(((float)i + 0.5f)* U);
		if (i < *n_bins_f1){
			f_d_Gy[i]					=	f1_d_Gy[i];
			f_dd_Gy[i]					=	f1_dd_Gy[i];
			f_start[i]					=	f1[i];}
		else{
			f_d_Gy[i]					=	(float)sqrt(d_low[i] * d_high[i]);
			f_dd_Gy[i]					=	d_high[i] - d_low[i];
			f_start[i]					=	0.0f;}
	}

	free(d_low);
	free(d_high);

	return;
}

/*******************************************************************************
/ Successive convolutions (Kellerer Algorithm)
*******************************************************************************/

typedef struct{

	long			array_size;									// size of function arrays F..BI
	long			N2;
	float			U;

	float			X;
	float			FINAL;


	float			CN;
	long			N1;


	float			CM1;
	float			CM2;
	float			CM3;
	float			CM4;

	float			D1;
	float			D2;
	float			D3;
	float			D4;

	float			F0;
	float*			F;
	long			MIF;
	long			LEF;

	float			H0;
	float*			H;
	long			MIH;
	long			LEH;

	float			E0;
	float*			E;
	float*			DE;
	long			MIE;

	float*			DI;
	float*			A;
	float*			BI;

	bool			write_output;
	FILE*			output_file;

	bool			shrink_tails;
	float			shrink_tails_under;
	bool			adjust_N2;
}   	aKList;


///////////////////////////////////////////////////////////////////////////////////////
// SGP_SC_NORMAL
///////////////////////////////////////////////////////////////////////////////////////
aKList	SGP_SC_NORMAL(aKList theKList){

	if(theKList.write_output){
		fprintf(theKList.output_file,	"\n\nThis is subroutine SGP_SC_NORMAL\n");
		fprintf(theKList.output_file,	    "=========================\n");
	}

	float	Y					= theKList.CM1 * 2;
	float	Z					= theKList.CM2 * 2;
	float	CM0					= theKList.H0;
			theKList.CM1		= 0;

	long		N					= theKList.MIH - theKList.MIE;
	long 		L;
	for (L = 1; L <= theKList.LEH; L++){
		long		LE					=		L + N;
		float	S					=		theKList.H[L-1] * theKList.DE[LE-1];
				CM0					=		CM0 + S;
				theKList.CM1		=		theKList.CM1 + S * theKList.E[LE-1];
		}

	float	TT					=		(1.0f - theKList.H0) / (CM0 - theKList.H0);
			theKList.CM1		=		theKList.CM1 * TT;
	float	R					=		theKList.CM1 * theKList.CM1;
			theKList.CM2		=		R * theKList.H0;
			theKList.CM3		=		-1.0f * theKList.CM1 * R * theKList.H0;
			theKList.CM4		=		R * R * theKList.H0;

	for (L = 1; L <= theKList.LEH; L++){
		long		LE					=		L + N;
		float	EC					=		theKList.E[LE-1] - theKList.CM1;
		float	E2					=		EC * EC;
				theKList.H[L-1]		=		theKList.H[L-1] * TT;
		float	S					=		theKList.H[L-1] * theKList.DE[LE-1] * E2;
				theKList.CM2		=		theKList.CM2 + S;
				theKList.CM3		=		theKList.CM3 + S * EC;
				theKList.CM4		=		theKList.CM4 + S * E2;
		}

	theKList.X			=	theKList.X * CM0;
	Y					=	theKList.CM1 / Y;
	Z					=	theKList.CM2 / Z;

	if(theKList.write_output){
		if(theKList.N1 > 0){
			fprintf(theKList.output_file,	"\nPrecision control (ratio actual/theoretical):\n");
			fprintf(theKList.output_file,	"Norm:\t\t %4.3e\n", theKList.X);
			fprintf(theKList.output_file,	"Mean:\t\t %4.3e\n", Y);
			fprintf(theKList.output_file,	"Var:\t\t %4.3e\n", Z);
		}else{
			fprintf(theKList.output_file,	"\nSingle hit distribution integral:\t\t %4.3e\n", CM0);
		}
	}

	return(theKList);
}



///////////////////////////////////////////////////////////////////////////////////////
// SGP_SC_OUTPUT
///////////////////////////////////////////////////////////////////////////////////////
aKList	SGP_SC_OUTPUT(aKList theKList){

	if(theKList.write_output){
		fprintf(theKList.output_file,	"\n\nThis is subroutine SGP_SC_OUTPUT\n");
		fprintf(theKList.output_file,	    "=========================\n");
	}

	float*	SD						=	(float*)calloc(theKList.array_size, sizeof(float));

	float	B						=	theKList.CM2 / (theKList.CM1 * theKList.CM1);
	float	C						=	theKList.CM3 / (float)sqrt(theKList.CM2 * theKList.CM2 * theKList.CM2);
	float	D						=	theKList.CM4 / (theKList.CM2 * theKList.CM2);
	float	S1						=	theKList.CN * theKList.D1;
	float	S2						=	theKList.CN * theKList.D2 / (S1 * S1);
	float	S3						=	theKList.D3 / (float)sqrt(theKList.CN * theKList.D2 * theKList.D2 * theKList.D2);
	float	S4						=	theKList.D4 / (theKList.D2 * theKList.D2 * theKList.CN) + 3;

	if(theKList.N1 <= 0){
		S2								=	B;
		S3								=	C;
		S4								=	D;}

	if(theKList.write_output){
		fprintf(	theKList.output_file,
					"\nZero component:\t\t%4.3e\n",
					theKList.H0);

		fprintf(	theKList.output_file,
					"\nMean\t\t\t\t\tActual:\t%4.3e\tTheory:\t%4.3e\n",
					theKList.CM1, S1);
		fprintf(	theKList.output_file,
					"Variance/Mean^2\t\t\tActual:\t%4.3e\tTheory:\t%4.3e\n",
					B, S2);
		fprintf(	theKList.output_file,
					"Central3/Variance^3/2\tActual:\t%4.3e\tTheory:\t%4.3e\n",
					C, S3);
		fprintf(	theKList.output_file,
					"Central4/Variance^2\t\tActual:\t%4.3e\tTheory:\t%4.3e\n",
					D, S4);

		fprintf(	theKList.output_file,
					"\nMIF: %d, LEF: %d, MIH: %d, LEH: %d, MIE: %d\n\n",
					theKList.MIF,
					theKList.LEF,
					theKList.MIH,
					theKList.LEH,
					theKList.MIE);

		fprintf(	theKList.output_file,
					"i\tE\t\t\tDE\t\t\tH\t\t\tF\n");

		long 	L;
		for (L = 1; L <= theKList.array_size; L++){
			fprintf(	theKList.output_file,
						"%d\t%4.3e\t%4.3e\t%4.3e\t%4.3e\n",
						L,
						theKList.E[L-1],
						theKList.DE[L-1],
						theKList.H[L-1],
						theKList.F[L-1]);
		}
	}

	return(theKList);
}



///////////////////////////////////////////////////////////////////////////////////////
// SGP_SC_INTERP
///////////////////////////////////////////////////////////////////////////////////////

aKList	SGP_SC_INTERP(aKList theKList){

	theKList.A[1-1]					=	theKList.F[2-1] - theKList.F[1-1];
	theKList.BI[1-1]				=	0.0f;
	theKList.F[theKList.LEF + 1 - 1]=	0.0f;

	long 	K;
	for (K = 1; K <= theKList.N2; K++){
		long L					=	theKList.LEF + K;
		theKList.A[L-1]			=	0.0f;
		theKList.BI[L-1]		=	0.0f;
	}

	long 	L;
	for (L = 2; L <= theKList.LEF; L++){
		theKList.A[L -1]		=	0.5f * (theKList.F[L + 1 -1] - theKList.F[L - 1 -1]);
		theKList.BI[L -1]		=	theKList.A[L-1] + theKList.F[L - 1 -1] - theKList.F[L -1];
	}

	return(theKList);
}

///////////////////////////////////////////////////////////////////////////////////////
// SGP_SC_RESET
///////////////////////////////////////////////////////////////////////////////////////

aKList	SGP_SC_RESET(aKList theKList){

	if (theKList.N2 <= 256){
		if(theKList.LEF <= 64){

/*			/////////////////////////////////////////////////////////////////////////////
			FILE* checkResetFile	=	fopen("CheckReset.log","w");

			fprintf(	checkResetFile,
						"BEFORE RESET:\nN2: %d, MIF: %d, LEF: %d, MIH: %d, LEH: %d, MIE: %d, E0: %4.3g\n\n",
						theKList.N2,
						theKList.MIF,
						theKList.LEF,
						theKList.MIH,
						theKList.LEH,
						theKList.MIE,
						theKList.E0);

			fprintf(	checkResetFile,
						"i\tE\t\t\tDE\t\t\tH\t\t\tF\n");

			for (long L = 1; L <= theKList.array_size; L++){
				fprintf(	checkResetFile,
							"%d\t%4.3e\t%4.3e\t%4.3e\t%4.3e\n",
							L,
							theKList.E[L-1],
							theKList.DE[L-1],
							theKList.H[L-1],
							theKList.F[L-1]);}

			///////////////////////////////////////////////////////////////////////////// */

			float S							=	(float)log(2.0f);
			float TT						=	(float)theKList.N2;
//			theKList.N2						=	theKList.N2 * 2;
			theKList.N2						+=	(long)(0.1 + exp((float)((long)(log(TT) / S - 0.99f)) * S));
			TT								=	TT / (float)theKList.N2;
			theKList.U						=	S / (float)theKList.N2;

			if(theKList.write_output){
				fprintf(theKList.output_file,			"\n\nThis is subroutine SGP_SC_RESET\n");
				fprintf(theKList.output_file,			"========================\n");
				fprintf(theKList.output_file,			"\nCoordinate change with new N2: %d\n", theKList.N2);
			}

			theKList						=	SGP_SC_INTERP(theKList);
			theKList.F[theKList.LEF + 1 -1]	=	0;
			long N							=	theKList.MIF;
			theKList.MIF					=	(long)((float)theKList.MIF / TT) + 1;		/////////////////////
			theKList.LEF					=	(long)((float)theKList.LEF / TT) - 1;		// added (SG) : -1 //
																							/////////////////////
			long 	K;
			for (K = 1; K <= theKList.LEF; K++){
				long   L						=	theKList.LEF - K + 1;
				float  FLF						=	(float)(L + theKList.MIF) * TT - (float)N;
				long   LFF						=	(long)(FLF + 0.5f);
				float S							=	FLF - (float)LFF;

				////////////////////////////////////////////////////////////////////////////////
				// Replaced Kellerer's original quadratic interpolation by a slower           //
				// but more correct approach. The original produced wrong interpolation &     //
				// negative values of F on the new N2's grid at very steep, irregular parts   //
				// of F	and eventually systematic deviation of moments.                       //
				////////////////////////////////////////////////////////////////////////////////

// original:	theKList.F[L -1]				=	theKList.F[LFF -1] + S * (theKList.A[LFF - 1] + S * theKList.BI[LFF - 1]);

				theKList.F[L -1]				=	theKList.F[LFF -1];
				if((S < 0 ) && (LFF >= 2)){
					theKList.F[L -1]			=	(float)(pow(theKList.F[LFF - 1 -1], -1.0f * S) * pow(theKList.F[LFF -1], 1.0f + S));
				}
				if((S > 0 ) && (LFF <= theKList.LEF - 1)){
					theKList.F[L -1]			=	(float)(pow(theKList.F[LFF -1], 1.0f - S) * pow(theKList.F[LFF + 1 -1], S));
				}
			}

			long 	L;
			for (L = theKList.N2; L <= theKList.array_size; L++){;
				float S					= (float)(L - theKList.N2) * theKList.U;
				float tmp				= (float)(-1.0f * log(1.0f - 0.5f * exp(-S)) / theKList.U);
				theKList.DI[L -1]		= tmp - (float)theKList.N2;		// type casts necessary to prevent round of errors (that will eventually yield negative H-values in SGP_SC_FOLD
			}

			theKList.MIE					=	theKList.MIF;

			long 	J;
			for (J = 1; J <= theKList.array_size; J++){
				float S						=	(float)(J + theKList.MIE);
				theKList.E[J -1]			=	(float)exp(S * theKList.U) * theKList.E0;
				///////////////////////////////////////////////////////////////////////////
				// addition SG: not to use Kellerer's formula for new DE's, as it is     //
				// not exact (but deviation are small)	                                 //
				///////////////////////////////////////////////////////////////////////////
				float* high_E				=	(float*)calloc(theKList.array_size, sizeof(float));
				S							=	(float)(J + theKList.MIE + 1);
				high_E[J - 1]				=	(float)exp(S * theKList.U) * theKList.E0;
				theKList.DE[J -1]			=	high_E[J -1] - theKList.E[J -1];
				free(high_E);
			}

/*			/////////////////////////////////////////////////////////////////////////////
			fprintf(	checkResetFile,
						"###################################################################\nAFTER RESET:\nN2: %d, MIF: %d, LEF: %d, MIH: %d, LEH: %d, MIE: %d, E0: %4.3g\n\n",
						theKList.N2,
						theKList.MIF,
						theKList.LEF,
						theKList.MIH,
						theKList.LEH,
						theKList.MIE,
						theKList.E0);

			fprintf(	checkResetFile,
						"i\tE\t\t\tDE\t\t\tH\t\t\tF\n");

			for ( L = 1; L <= theKList.array_size; L++){
				fprintf(	checkResetFile,
							"%d\t%4.3e\t%4.3e\t%4.3e\t%4.3e\n",
							L,
							theKList.E[L-1],
							theKList.DE[L-1],
							theKList.H[L-1],
							theKList.F[L-1]);}

			fclose(checkResetFile);
			///////////////////////////////////////////////////////////////////////////// */
		}else{
			return(theKList);
		}
	}else{
		theKList.MIE					=	theKList.MIF;

		long 	J;
		for (J = 1; J <= theKList.array_size; J++){
				float S						=	(float)(J + theKList.MIE);
				theKList.E[J -1]			=	(float)exp(S * theKList.U) * theKList.E0;
				///////////////////////////////////////////////////////////////////////////
				// addition SG: not to use Kellerer's formula for new DE's, as it is     //
				// not exact (but deviation are small)	                                 //
				///////////////////////////////////////////////////////////////////////////
				float* high_E				=	(float*)calloc(theKList.array_size, sizeof(float));
				S							=	(float)(J + theKList.MIE + 1);
				high_E[J - 1]				=	(float)exp(S * theKList.U) * theKList.E0;
				theKList.DE[J -1]			=	high_E[J -1] - theKList.E[J -1];
				free(high_E);
		}
	}

	return(theKList);
}



///////////////////////////////////////////////////////////////////////////////////////
// SGP_SC_ZERO
///////////////////////////////////////////////////////////////////////////////////////

aKList	SGP_SC_ZERO(aKList theKList){

	if(theKList.write_output){
		fprintf(theKList.output_file,			"\n\nThis is subroutine SGP_SC_ZERO\n");
		fprintf(theKList.output_file,			"========================\n");
	}

	theKList.X					=	0;
	long N						=	theKList.MIH - theKList.MIE;

	long 	L;
	for (L = 1; L <= theKList.LEH; L++){
		long K					=	L + N;
		theKList.X				=	theKList.X + theKList.H[L -1] * theKList.DE[K -1];
	}

	float S					=	(1.0f - theKList.F0) * (1.0f - theKList.F0) / theKList.X;
	theKList.X				=	2.0f / S;

	for (L = 1; L <= theKList.LEH; L++){;
		theKList.H[L -1]		=	theKList.H[L -1] * S;
	}

	N							=	theKList.MIH - theKList.MIF;
	theKList.MIH				=	theKList.MIF;
	theKList.LEH				=	theKList.LEH + N;

	long 	LL;
	for (LL = 1; LL <= theKList.LEH; LL++){
		long L						=	theKList.LEH + 1 - LL;
		long K						=	L + N;
		theKList.H[K -1]			=	theKList.H[L -1];
	}

	for (L = 1; L <= N; L++){
		theKList.H[L -1]			=	0.0f;
	}

	S							=	theKList.F0 * 2.0f;

	for (L = 1; L <= theKList.LEF; L++){
		theKList.H[L -1]			=	theKList.H[L -1] + theKList.F[L -1] * S;
	}

	return(theKList);

}



///////////////////////////////////////////////////////////////////////////////////////
// SGP_SC_SHRINK
///////////////////////////////////////////////////////////////////////////////////////

aKList	SGP_SC_SHRINK(aKList theKList){

	float	EX						=	theKList.shrink_tails_under;
	float	S						=	0.0f;
	long	N						=	theKList.MIH - theKList.MIE;

	long 	L;
	for (L = 1; L <= theKList.LEH; L++){
		long K							=	L + N;
		S								=	S + theKList.H[L -1] * theKList.DE[K -1];
		if(S > 1000.0 * EX){
			theKList.MIH 				=	theKList.MIH + L - 1;
			break;}
		}

	long		M						=	L - 1;
	S								=	0;

	long 	K;
	for (K = 1; K <= theKList.LEH; K++){
		L								=	theKList.LEH + 1 - K;
		long KK							=	L + N;
		S								=	S + theKList.H[L - 1] * theKList.DE[KK - 1];
		if(S > EX){
			break;}
	}

	theKList.LEH					=	L - M;
	for (L = 1; L <= theKList.LEH; L++){
		K								=	L + M;
		theKList.H[L -1]				=	theKList.H[K -1];
	}

	K								=	theKList.LEH + 1;
	long	KK							=	theKList.LEH + M;
	for (L = K; L <= KK; L++){
		theKList.H[L -1]				=	0;
	}

	return(theKList);

}


///////////////////////////////////////////////////////////////////////////////////////
// SGP_SC_FOLD
///////////////////////////////////////////////////////////////////////////////////////

aKList	SGP_SC_FOLD(aKList theKList){

	float*	FDE					=	(float*)calloc(theKList.array_size, sizeof(float));

	if((theKList.CN >= 10.0) && (theKList.adjust_N2 == true)){
		theKList					=	SGP_SC_RESET(theKList);}

	theKList.H0					=	theKList.F0 * theKList.F0;
	theKList.MIH				=	theKList.MIF + theKList.N2;
	theKList.LEH				=	theKList.LEF;
	long	K					=	theKList.LEF + 1;
	long KK						=	K + theKList.N2;

	long 	L;
	for (L = K; L <= KK; L++){
		theKList.F[L -1]			=	0;
	}

	theKList					=	SGP_SC_INTERP(theKList);
	long N						=	theKList.MIF - theKList.MIE;

	for (L = 1; L <= theKList.LEH; L++){
		K							=	L + N;
		FDE[L -1]					=	theKList.F[L -1] * theKList.DE[K -1];
	}

	long 	LH;
	for (LH = 1; LH <= theKList.LEH; LH++){
		float 	HLH					=	0;
		long 	LL					=	LH + theKList.N2;
		long 	LF;
		for (LF = 1; LF <= LH; LF++){
			K							=	LL - LF;
			float FLF					=	(float)LH - theKList.DI[K -1];
			long LFF					=	(long)(FLF + 0.5f);
			float S						=	FLF - (float)LFF;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Modification SG: if Kellerer's quadratic interpolation fails, use simple estimate
			float tmp					=	theKList.F[LFF -1] + S * (theKList.A[LFF -1] + S * theKList.BI[LFF -1]);
			if (tmp <0){
				tmp = 0.0f;				// Very crude - better to replace by interpolation as done in RESET
			}							// which is time-consuming, however.
			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			HLH							=	HLH + FDE[LF -1] * tmp;
		}
		theKList.H[LH -1]			=	HLH - FDE[LH -1] * theKList.F[LH -1] * 0.5f;
	}

	if (theKList.F0 < 1e-10){
		theKList.X						=	2.0f;
	}else{
		theKList						=	SGP_SC_ZERO(theKList);
	}

	return(theKList);

}



void	 SGP_SuccessiveConvolutions(				float*	u,
													long*	n_bins_f,
													long*	N2,
													long*	n_bins_f_used,
													float*	f_d_Gy,
													float*	f_dd_Gy,
													float*	f,
													float*	f0,
													float*	fdd,									// frequence:			H * DE			(f * dd)
													float*	dfdd,									// dose contribution:	H * E * DE		(f * d * dd)
													float*	d,										// first moment:						(<d>)
													bool*	write_output,
													bool*	shrink_tails,
													float*	shrink_tails_under,
													bool*	adjust_N2)
{
	// index variables
	long		i;


	//////////////////////////////////////////
	// Init KList structure
	// (Constructor)
	//////////////////////////////////////////
	aKList				KList;



	KList.array_size	= *n_bins_f;

	KList.N2			= *N2;
	KList.U				= (float)log(2.0f) / KList.N2;

	//////////////////////////////////

	KList.write_output		=	*write_output;
	if(KList.write_output){
		KList.output_file		=	fopen("SuccessiveConvolutions.log","w");
		if (KList.output_file == NULL) return;											// File error

		fprintf(KList.output_file, "##############################################################\n");
		fprintf(KList.output_file, "##############################################################\n");
		fprintf(KList.output_file, "This is LGC2.2 core - successive convolution mode (2008/08/12).\n");
	}

	//////////////////////////////////

	KList.shrink_tails		=	*shrink_tails;
	KList.shrink_tails_under=	*shrink_tails_under;
	KList.adjust_N2			=	*adjust_N2;

	//////////////////////////////////

	KList.F				= (float*)calloc(KList.array_size, sizeof(float));
	KList.H				= (float*)calloc(KList.array_size, sizeof(float));
	KList.E				= (float*)calloc(KList.array_size, sizeof(float));
	KList.DE			= (float*)calloc(KList.array_size, sizeof(float));
	KList.DI			= (float*)calloc(KList.array_size, sizeof(float));
	KList.A				= (float*)calloc(KList.array_size, sizeof(float));
	KList.BI			= (float*)calloc(KList.array_size, sizeof(float));

	// Some other initializations
	KList.MIH			= 0;
	KList.MIE			= 0;
	KList.N1			= 0;
	KList.H0			= 0;
	KList.X				= 1;
	KList.CN			= 1;
	KList.CM1			= 1;
	KList.CM2			= 1;


	if(KList.write_output){
		fprintf(KList.output_file,	"\n\nThis is main\n");
		fprintf(KList.output_file,	    "============\n");
	}

	// Copy input data
	KList.E0			= f_d_Gy[1 -1] * (float)exp(-1.0f * KList.U);

	long 	L;
	for (L = 1; L <= KList.array_size; L++){
		KList.E[L -1]			= f_d_Gy[L -1];
		KList.DE[L -1]			= f_dd_Gy[L -1];
		KList.H[L -1]			= f[L -1];}

	KList.LEH				= *n_bins_f_used;

	///////////////////////////////////////
	// Fill array for auxilary function that enables easy index operations
	for	(L = KList.N2; L <= KList.array_size; L++){
		float S				= (float)(L - KList.N2) * KList.U;
		float tmp			= (float)(-1.0f * log(1.0f - 0.5f * exp(-S)) / KList.U);
		KList.DI[L -1]		= tmp - (float)KList.N2;}		// type casts necessary to prevent round of errors (that will eventually yield negative H-values in SGP_SC_FOLD

	///////////////////////////////////////
	// Normalize distribution
	///////////////////////////////////////
	KList	= SGP_SC_NORMAL(KList);

	if(KList.write_output){
		fprintf(KList.output_file,	"\n\nThis is main\n");
		fprintf(KList.output_file,	    "============\n");

		fprintf(KList.output_file, "\nNormalized single hit distribution in KList:\n");
		for (i = 0; i < KList.array_size; i++){
			fprintf(	KList.output_file,
						"i: %d\t\tKList.E: %4.3e Gy\t\tKList.DE: %4.3e\t\tKList.H: %4.3e\t\t\n",
						i,
						KList.E[i],
						KList.DE[i],
						KList.H[i]);
		}

		fprintf(KList.output_file,	"\n\nThis is main\n");
		fprintf(KList.output_file,	    "============\n");

		fprintf(KList.output_file, "\nMoments of the single hit distribution:\n");
	}

	///////////////////////////////////////
	// Get moments of single impact f1
	///////////////////////////////////////

			KList.D1		=		KList.CM1;
	float	S				=		KList.D1 * KList.D1;
			KList.D2		=		KList.CM2 + S;
			KList.D3		=		KList.CM3 + 3.0f * KList.CM2 * KList.D1 + S * KList.D1;
			KList.D4		=		KList.CM4 + 4.0f * KList.CM3 * KList.D1 + 6.0f * S * KList.CM2 + S * S;

	float	S2				=		KList.D2 / KList.D1;
	float	S3				=		KList.D3 / KList.D1;
	float	S4				=		KList.D4 / KList.D1;
			S				=		S3 / (float)sqrt(S2 * S2 * S2);
	float	TT				=		S4	/ (S2 * S2);

	if(KList.write_output){
		fprintf(	KList.output_file,
					"\nInitial distribution:\n");
		fprintf(	KList.output_file,
					"Delta 1:\t\t%4.3e\n",
					KList.D1);
		fprintf(	KList.output_file,
					"Delta 2:\t\t%4.3e\n",
					KList.D2);

		fprintf(	KList.output_file,
					"\nCharacteristics of the solution to the mean value E\n");
		fprintf(	KList.output_file,
					"Rel. variance:\t%4.3e / E\n",
					S2);
		fprintf(	KList.output_file,
					"Skewness:\t\t%4.3e / sqrt(E)\n",
					S);
		fprintf(	KList.output_file,
					"Kurtosis:\t\t%4.3e / E + 3\n",
					TT);
	}

	///////////////////////////////////////
	// SGP_SC_SHRINK distribution
	///////////////////////////////////////
	if(KList.shrink_tails){
		KList	= SGP_SC_SHRINK(KList);}

	///////////////////////////////////////
	// SGP_SC_OUTPUT
	///////////////////////////////////////
	KList	= SGP_SC_OUTPUT(KList);

	///////////////////////////////////////
	// Get approximation for small hit
	// numbers
	///////////////////////////////////////

	KList.FINAL				= *u * KList.D1;													// Final mean impact number

	long	n_convolutions		= 0;
	KList.CN				=		KList.FINAL	/ KList.D1;
	while(KList.CN > MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT){
		KList.CN				=		0.5f * KList.CN;
		n_convolutions++;
	}

	if(KList.write_output){
		fprintf(KList.output_file,	"\n\nThis is main\n");
		fprintf(KList.output_file,	    "============\n");

		fprintf(	KList.output_file,	"\nSmall hit number approximation:\n");
		fprintf(	KList.output_file,
					"\nTarget hit value:\t%4.3e\t\tStart hit value:\t%4.3e\t\tNumber of convolutions:\t%d\n",
					*u,
					KList.CN,
					n_convolutions);
	}

	KList.H0				=		1.0f - KList.CN;

	for (L = 1; L <= KList.LEH; L++){
		KList.H[L -1]			=	KList.H[L -1] * KList.CN;
	}


	///////////////////////////////////////
	// Convolution loop
	///////////////////////////////////////
	long 	j;
	for(j = 0; j < n_convolutions; j++){
		KList.N1				=	KList.N1 + 1;
		KList.CN				=	KList.CN * 2;

		if(KList.write_output){
			fprintf(	KList.output_file,
						"\n\n##############################################################\n");
			fprintf(	KList.output_file,	"\n\nThis is main\n");
			fprintf(	KList.output_file,	    "============\n");

			fprintf(	KList.output_file,
						"\nConvolution number:\t\t%d\nMean hit number:\t\t%4.3e\n",
						KList.N1,
						KList.CN);
		}

		for (L = 1; L <= KList.LEH; L++){
			KList.F[L -1]			=	KList.H[L -1];}

		KList.F0				=	KList.H0;
		KList.LEF				=	KList.LEH;
		KList.MIF				=	KList.MIH;
		KList					=	SGP_SC_FOLD(KList);
		if(KList.shrink_tails){
			KList					=	SGP_SC_SHRINK(KList);}
		KList					=	SGP_SC_NORMAL(KList);
		KList					=	SGP_SC_OUTPUT(KList);
	}


	//////////////////////////////////////////
	// Copy results back to input structure
	// and adjust according to MIH, MIE
	//////////////////////////////////////////


	*d		= 0.0f;

	for (L = 1; L <= KList.array_size; L++){
		f_d_Gy[L -1]			=	0.0f;
		f_dd_Gy[L -1]			=	0.0f;
		f[L -1]					=	0.0f;
		fdd[L -1]				=	0.0f;
		dfdd[L -1]				=	0.0f;}

	long	N				= KList.MIH - KList.MIE;
	for (L = 1; L <= KList.LEH; L++){
			long LE					=	L + N;
			f_d_Gy[L -1]			=	KList.E[LE -1];
			f_dd_Gy[L -1]			=	KList.DE[LE -1];
			f[L -1]					=	KList.H[L-1];
			fdd[L -1]				=	f[L -1] * f_dd_Gy[L -1];
			dfdd[L -1]				=	fdd[L -1] * f_d_Gy[L -1];
			*d						+=	dfdd[L -1];}

	*n_bins_f_used = KList.LEH;

	*f0						= KList.H0;

	*N2						= KList.N2;			// could have been changed by RESET --> report back

	//////////////////////////////////////////
	// Free allocated KList structures
	// (Destructor)
	//////////////////////////////////////////
	free(KList.F);
	free(KList.H);
	free(KList.E);
	free(KList.DE);
	free(KList.DI);
	free(KList.A);
	free(KList.BI);

	if(KList.write_output){
		fprintf(KList.output_file, "\n\nThis is main\n");
		fprintf(KList.output_file, "============\n");
		fprintf(KList.output_file, "\nDone.\n");
		fprintf(KList.output_file, "##############################################################\n");
		fprintf(KList.output_file, "##############################################################\n");

		// Close file
		fclose(KList.output_file);
	}

	return;
}


void SGP_SC_Loop(	long*	n,
					float*	E_MeV_u,
					long*	particle_no,
					float*	fluence_cm2,
					long*	slab_no,
					char*	material_name,
					float*	parameter,
					long*	N2,
					long*	n_slabs,
					long*	n_gamma_parameter,
					long*	gamma_model,
					float*	gamma_parameter,
					long*	verbose_level,
					char*	output_fileName,
					// return values
					float*	u,
					float*	total_d_Gy,
					float*	d,
					float*	S_HCP,
					float*	S_gamma,
					float*	efficiency,
					float*	S_HCP_total,
					float*	S_gamma_total,
					float*	efficiency_total)

{
	FILE* output_file	=	fopen(output_fileName,"w");
	if (output_file == NULL) return;											// File error

	*S_HCP_total		=	0.0f;
	*S_gamma_total		=	0.0f;
	*efficiency_total	=	0.0f;

	long 	i, j;
	if(*verbose_level == 0){
		fprintf(output_file,	"LCG3 run\n");

		fprintf(output_file,	"N2, material, r.min.m");

		for (i = 0; i < *n_gamma_parameter; i++){
			fprintf(output_file,	", parameter.%d", i);}
		fprintf(output_file,	"\n%d, \"%s\", %4.3g",
								*N2, material_name, parameter[0]);
		fprintf(output_file,	", %4.3g", gamma_parameter[0]);
		for (i = 1; i < *n_gamma_parameter; i++){
			fprintf(output_file,	", %4.3g", gamma_parameter[i]);}
		fprintf(output_file,	"\n");
		fprintf(output_file,	"\nslab.no, E.ave.MeV, E.dw.MeV, LET.ave.MeV.cm2.g, LET.dw.MeV.cm2.g, ");
		fprintf(output_file,	"total.fluence.cm2, mean.hit.number, fluence.factor, d.spectrum.Gy, ");
		fprintf(output_file,	"d.set.Gy, dose.SC.Gy, n.bins.used, S.HCP, S.gamma, efficiency\n");
	}
	else{
		fprintf(output_file,	"###################################################################\n");
		fprintf(output_file,	"# This is LGC3 - SC loop\n");
		fprintf(output_file,	"###################################################################\n");
		fprintf(output_file,	"\nLooping over %d uniform slabs (in beam direction) of same size of a extended detector.\n",
								*n_slabs);
		fprintf(output_file,	"\n%d steps within a factor of two (N2).\n",
								*N2);
		fprintf(output_file,	"\nMaterial: %s.\n",
								material_name);
		fprintf(output_file,	"\nRadial dose distribution according to Geiß (1997), r.min = %3.2g m.\n",
								parameter[0]);
		fprintf(output_file,	"\nUsing gamma response model %d with parameters (",
								*gamma_model);

		for (i = 0; i < *n_gamma_parameter; i++){
			if(i != 0)	{	fprintf(output_file,	", %4.3f",	gamma_parameter[i]);}
			else		{	fprintf(output_file,	"%4.3f",	gamma_parameter[i]);}
		}

		fprintf(output_file,	")\n\n");
	}

	for (i = 0; i < *n_slabs; i++){

		// How many particles for this slab?
		long	nLines						=	0;
		for (j = 0; j < *n; j++){
			if (slab_no[j]	== i+1){
				nLines++;
			}
		}


		if(nLines > 0){
			float*	E_MeV_u_slab;
			if (E_MeV_u[0] > 0){																			// If E < 0 are passed, non PSTAR-LET are provided, see SGP_RDD_Parameters
				E_MeV_u_slab				=	(float*)calloc(nLines, sizeof(float));}
			else{
				E_MeV_u_slab				=	(float*)calloc(2 * nLines, sizeof(float));}
			long*	particle_no_slab			=	(long*)calloc(nLines, sizeof(long));
			float*	fluence_cm2_slab			=	(float*)calloc(nLines, sizeof(float));

			// extract E, F, p for current slab
			long	m	=	0;
			for (j = 0; j < *n; j++){
				if (slab_no[j]	== i+1){
					if(E_MeV_u[0] > 0){																	// If E < 0 are passed, non PSTAR-LET are provided, see SGP_RDD_Parameters
						E_MeV_u_slab[m]			=	E_MeV_u[j];}
					else{
						E_MeV_u_slab[m]			=	E_MeV_u[j];
						E_MeV_u_slab[m+nLines]	=	E_MeV_u[j+(*n)];
					}
					particle_no_slab[m]		=	particle_no[j];
					fluence_cm2_slab[m++]	=	fluence_cm2[j];
				}
			}

			// Get array size for f1
			long	n_bins_f1;
			bool	debug	= false;
			SGP_SC_get_f1_array_size(	&nLines,
										E_MeV_u_slab,
										particle_no_slab,
										material_name,
										parameter,
										N2,
										// from here: return values
										&n_bins_f1,
										&debug);

			// Allocate memory for f1 and f1-particle information return data
			float*	f1_d_Gy						=	(float*)calloc(n_bins_f1, sizeof(float));
			float*	f1_dd_Gy					=	(float*)calloc(n_bins_f1, sizeof(float));
			float*	f1							=	(float*)calloc(n_bins_f1, sizeof(float));

			float*	norm_fluence				=	(float*)calloc(nLines, sizeof(float));
			float*	LET_MeV_cm2_g				=	(float*)calloc(nLines, sizeof(float));
			float*	r_min_m						=	(float*)calloc(nLines, sizeof(float));
			float*	r_max_m						=	(float*)calloc(nLines, sizeof(float));
			float*	d_min_Gy					=	(float*)calloc(nLines, sizeof(float));
			float*	d_max_Gy					=	(float*)calloc(nLines, sizeof(float));
			float*	k_Gy						=	(float*)calloc(nLines, sizeof(float));
			float*	single_impact_fluence_cm2	=	(float*)calloc(nLines, sizeof(float));
			float*	single_impact_dose_Gy		=	(float*)calloc(nLines, sizeof(float));
			float*	dose_contribution_Gy		=	(float*)calloc(nLines, sizeof(float));

			float	ave_E_MeV					=	0.0f;
			float	dw_E_MeV					=	0.0f;

			float	ave_LET_MeV_cm2_g			=	0.0f;
			float	dw_LET_MeV_cm2_g			=	0.0f;

			float	total_fluence_cm2			=	0.0f;

			// ...and get f1
			SGP_SC_get_f1(	&nLines,
							E_MeV_u_slab,
							particle_no_slab,
							fluence_cm2_slab,
							material_name,
							parameter,
							N2,
							&n_bins_f1,
							// from here: return values
							norm_fluence,
							LET_MeV_cm2_g,
							r_min_m,
							r_max_m,
							d_min_Gy,
							d_max_Gy,
							k_Gy,
							single_impact_fluence_cm2,
							single_impact_dose_Gy,
							dose_contribution_Gy,
							&total_fluence_cm2,
							&total_d_Gy[i],
							&ave_E_MeV,
							&dw_E_MeV,
							&ave_LET_MeV_cm2_g,
							&dw_LET_MeV_cm2_g,
							&u[i],
							f1_d_Gy,
							f1_dd_Gy,
							f1,
							&debug);


			float	fluence_factor	=	1.0f;

			// Get array size for f, u_start and number
			// of convolutions
			long	n_bins_f;
			float	u_start;
			long	n_convolutions;
			SGP_SC_get_f_array_size(	u,
										&fluence_factor,
										N2,
										&n_bins_f1,
										f1_d_Gy,
										f1_dd_Gy,
										f1,
										// from here: return values
										&n_bins_f,
										&u_start,
										&n_convolutions);

			// Allocate memory for f(_start)
			float*	f_d_Gy						=	(float*)calloc(n_bins_f, sizeof(float));
			float*	f_dd_Gy						=	(float*)calloc(n_bins_f, sizeof(float));
			float*	f							=	(float*)calloc(n_bins_f, sizeof(float));

			// Get f_start
			SGP_SC_get_f_start(	&u_start,
								&n_bins_f1,
								N2,
								f1_d_Gy,
								f1_dd_Gy,
								f1,
								&n_bins_f,
								// from here: return values
								f_d_Gy,
								f_dd_Gy,
								f);

			// Allocate memory for convolutions
			float*	fdd							=	(float*)calloc(n_bins_f, sizeof(float));
			float*	dfdd						=	(float*)calloc(n_bins_f, sizeof(float));
			float	f0							=	0.0f;
			bool	write_output				=	false;

			bool	shrink_tails				=	true;
			float	shrink_tails_under			=	1e-40f;

			bool	adjust_N2					=	true;

			long	n_bins_f_used				=	n_bins_f1;

			SGP_SuccessiveConvolutions(	&u[i],
										&n_bins_f,
										N2,
										&n_bins_f_used,
										f_d_Gy,
										f_dd_Gy,
										f,
										&f0,
										fdd,
										dfdd,
										&d[i],
										&write_output,
										&shrink_tails,
										&shrink_tails_under,
										&adjust_N2);


			// Allocate memory for gamma response
			float*	S							=	(float*)calloc(n_bins_f, sizeof(float));

			// Apply gamma response
			SGP_get_gamma_response(	&n_bins_f_used,				// only the first XXX bins are used
									f_d_Gy,
									f_dd_Gy,
									f,
									&f0,
									n_gamma_parameter,
									gamma_model,
									gamma_parameter,
									S,
									// return
									&S_HCP[i],
									&S_gamma[i],
									&efficiency[i]);

			if(!_isnan(S_HCP[i])){
				*S_HCP_total		+=	S_HCP[i];}
			if(!_isnan(S_gamma[i])){
				*S_gamma_total		+=	S_gamma[i];}
			if(!_isnan(S_HCP[i]) && S_HCP[i] > 0){
				*efficiency_total	+=	S_HCP[i];}

			// Output results

			if(*verbose_level == 0){
				fprintf(output_file,	"%4d, %4.3g, %4.3g, %4.3g, %4.3g, %4.3g, %4.3g, %4.3g, %4.3g, %4.3g, %4.3g, %4d, %4.3g, %4.3g, %4.3g\n",
										i+1, ave_E_MeV, dw_E_MeV, ave_LET_MeV_cm2_g, dw_LET_MeV_cm2_g,
										total_fluence_cm2, u[i], fluence_factor, total_d_Gy[i], total_d_Gy[i] * fluence_factor, d[i],
										n_bins_f_used, S_HCP[i], S_gamma[i], efficiency[i]);
			}
			else{
				fprintf(output_file,	"\n###################################################################\n");
				fprintf(output_file,	"Slab no. %d, %d entries in particle spectrum\n",
										i + 1,
										nLines);
				fprintf(output_file,	"\n dose from spectrum/ Gy:\t\t\t%4.3g\n",			total_d_Gy[i]);
				fprintf(output_file,	" average E / MeV:\t\t\t\t\t%4.3g\n",				ave_E_MeV);
				fprintf(output_file,	" dose-weigthed E / MeV):\t\t\t%4.3g\n",			dw_E_MeV);
				fprintf(output_file,	" average LET / (MeV*cm2/g):\t\t\t%4.3g\n",			ave_LET_MeV_cm2_g);
				fprintf(output_file,	" dose-weigthed LET / (MeV*cm2/g):\t%4.3g\n",		dw_LET_MeV_cm2_g);
				fprintf(output_file,	" total fluence / cm-2:\t\t\t\t%4.3g\n",			total_fluence_cm2);
				fprintf(output_file,	" mean hit number µ:\t\t\t\t\t%4.3g\n",				u[i]);
				fprintf(output_file,	" fluence factor:\t\t\t\t\t%4.3g\n",				fluence_factor);
				fprintf(output_file,	" dose set/ Gy:\t\t\t\t\t\t%4.3g\n",				total_d_Gy[i] * fluence_factor);
				fprintf(output_file,	"\n average dose from SC / Gy:\t\t\t%4.3g\n",				d[i]);
				fprintf(output_file,	"\n HCP response:\t\t\t\t\t\t%4.3g\n",				S_HCP[i]);
				fprintf(output_file,	" gamma response:\t\t\t\t\t%4.3g\n",				S_gamma[i]);
				fprintf(output_file,	" efficiency:\t\t\t\t\t\t%4.3g\n",					efficiency[i]);
			}

			// free slab related memory again
			free(E_MeV_u_slab);
			free(particle_no_slab);
			free(fluence_cm2_slab);

			free(norm_fluence);
			free(LET_MeV_cm2_g);
			free(r_min_m);
			free(r_max_m);
			free(d_min_Gy);
			free(d_max_Gy);
			free(k_Gy);
			free(single_impact_fluence_cm2);
			free(single_impact_dose_Gy);

			free(f1_d_Gy);
			free(f1_dd_Gy);
			free(f1);

			free(f_d_Gy);
			free(f_dd_Gy);
			free(f);

			free(fdd);
			free(dfdd);

			free(S);
		}
		else{	// nLines == 0

			// Write: nix in slab, eff = 0
		}
	}

	*efficiency_total	/=	(*n_slabs * S_HCP[0]) ;

	if(*verbose_level > 0){
		fprintf(output_file,	"\n###################################################################\n");
		fprintf(output_file,	"\n total detector efficiency:\t\t\t%4.3g\n",				*efficiency_total);

		fprintf(output_file,	"\n###################################################################\n");
		fprintf(output_file,	"# Done.\n");
		fprintf(output_file,	"###################################################################\n");
	}

	fclose(output_file);
}


#endif // SGP_SUCCESSIVECONVOLUTIONS_H_
