#include <stdio.h>
#include <stdlib.h>

#define _LINUX // [_LINUX or _WINDOWS] : in Linux we have isnan function while in Windows we have _isnan
#define _R // [_S or _R] in S we can pass long type to the function via as.single, but in R we pass int type
//#define _DEBUG // debugging printouts
#define _SOLVER	// use SOLVER instead of analytical inversion in r_RDD_m

#include "SGP_Constants.h"
#include "SGP_Data.h"
#include "SGP_Utils.h"
#include "SGP_Functions.h"
#include "SGP_RDD.h"
#include "SGP_SuccessiveConvolutions.h"
#include "SGP_GammaResponse.h"
#include "SGP_FileOperations.h"
#include "SGP_ParabolicCylinderFunction.h"
#include "SGP_Transport.h"

#include <gsl/gsl_rng.h>

void SGP_efficiency(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						float*	fluence_cm2,
						long*	material_no,
						long*	RDD_model,
						float*	RDD_parameters,
						long*	ER_model,
						float*	ER_parameters,
						long*	gamma_model,
						float*	gamma_parameters,
//						method					= "SC",
						long*	N2,
						float*	fluence_factor,
						int*	write_output,
						int*	shrink_tails,
						float*	shrink_tails_under,
						int*	adjust_N2,
						float*	results);

void SGP_efficiency_grid(	long*	n,
							float*	E_MeV_u,
							long*	particle_no,
							float*	fluence_cm2,
							long*	material_no,
							long*	RDD_model,
							float*	RDD_parameters,
							long*	ER_model,
							float*	ER_parameters,
							long*	gamma_model,
							float*	gamma_parameters,
	//						method					= "grid",
							long*	N_runs,
							float*	fluence_factor,
							bool*	write_output,
							long*	nX,
							float*	grid_size_m,
							float*	results);

void SGP_efficiency(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						float*	fluence_cm2,
						long*	material_no,
						long*	RDD_model,
						float*	RDD_parameters,
						long*	ER_model,
						float*	ER_parameters,
						long*	gamma_model,
						float*	gamma_parameters,
//						method					= "SC",
						long*	N2,
						float*	fluence_factor,
						int*	write_output,
						int*	shrink_tails,
						float*	shrink_tails_under,
						int*	adjust_N2,
						float*	results)
{

#ifdef _DEBUG
	indnt_init();
	fprintf(debf,"%sSGP_efficiency\n",isp);
	fprintf(debf,"%swrite_output = %d\n",isp,*write_output);
	fprintf(debf,"%sshrink_tails = %d\n",isp,*shrink_tails);
	fprintf(debf,"%sadjust_N2 = %d\n",isp,*adjust_N2);
		//		fprintf(debf,"%sf1 = %g\n",isp,(1.0f/ (M_PI * (*a0_m)*(*a0_m))));
		//		fprintf(debf,"%sf2 = %g\n",isp,SGP_RDD_Katz_point_Gy(t_m,alpha,r_max_m,Katz_point_coeff_Gy));
		//		fprintf(debf,"%sf3 = %g\n",isp,geometryFunctionPhi(r_m,a0_m,t_m));
#endif


#ifdef _R
	bool write_output_b = false;
	if( *write_output == 1)
		write_output_b = true;

	bool shrink_tails_b = false;
	if( *shrink_tails == 1)
		shrink_tails_b = true;

	bool adjust_N2_b = false;
	if( *adjust_N2 == 1)
		adjust_N2_b = true;

	int n_int = (int)(*n);
	*n = (long)n_int;

	int N2_int = (int)(*N2);
	*N2 = (long)N2_int;

	int RDD_model_int = (int)(*RDD_model);
	*RDD_model = (long)RDD_model_int;

	int ER_model_int = (int)(*ER_model);
	*ER_model = (long)ER_model_int;

	int gamma_model_int = (int)(*gamma_model);
	*gamma_model = (long)gamma_model_int;

	int material_no_int	= (int)(*material_no);
	*material_no = (long)material_no_int;

	int particle_no_int = (int)(*particle_no);
	*particle_no = (long)particle_no_int;
#endif

#ifdef _DEBUG
	fprintf(debf,"%sSGP_efficiency 2\n",isp);
	fprintf(debf,"%swrite_output = %d\n",isp,*write_output);
	fprintf(debf,"%sshrink_tails = %d\n",isp,*shrink_tails);
	fprintf(debf,"%sadjust_N2 = %d\n",isp,*adjust_N2);
		//		fprintf(debf,"%sf1 = %g\n",isp,(1.0f/ (M_PI * (*a0_m)*(*a0_m))));
		//		fprintf(debf,"%sf2 = %g\n",isp,SGP_RDD_Katz_point_Gy(t_m,alpha,r_max_m,Katz_point_coeff_Gy));
		//		fprintf(debf,"%sf3 = %g\n",isp,geometryFunctionPhi(r_m,a0_m,t_m));
#endif


	long	n_bins_f1;
	float*	f1_parameters			=	(float*)calloc(9 * (*n), sizeof(float));

	SGP_SC_get_f1_array_size(	n,
								E_MeV_u,
								particle_no,
								material_no,
								RDD_model,
								RDD_parameters,
								ER_model,
								ER_parameters,
								N2,
								&n_bins_f1,
								f1_parameters);

	float*	f_parameters			=	(float*)calloc(7, sizeof(float));

	float*	norm_fluence			=	(float*)calloc(*n, sizeof(float));
	float*	dose_contribution_Gy	=	(float*)calloc(*n, sizeof(float));

	float*	f1_d_Gy					=	(float*)calloc(n_bins_f1, sizeof(float));
	float*	f1_dd_Gy				=	(float*)calloc(n_bins_f1, sizeof(float));
	float*	f1						=	(float*)calloc(n_bins_f1, sizeof(float));

	SGP_SC_get_f1(	n,
					E_MeV_u,
					particle_no,
					fluence_cm2,
					/* detector parameters */
					material_no,
					RDD_model,
					RDD_parameters,
					/* electron range model */
					ER_model,
					ER_parameters,
					/* algorith parameters*/
					N2,
					&n_bins_f1,
					/* f1 parameters*/
					f1_parameters,
					// from here: return values
					norm_fluence,
					dose_contribution_Gy,
					f_parameters,
					f1_d_Gy,
					f1_dd_Gy,
					f1);

	long			n_bins_f;
	float			u_start;
	long			n_convolutions;


	SGP_SC_get_f_array_size(	&f_parameters[0],			// = u
			fluence_factor,
			N2,
			&n_bins_f1,
			f1_d_Gy,
			f1_dd_Gy,
			f1,
			// from here: return values
			&n_bins_f,
			&u_start,
			&n_convolutions);

	float*	f_d_Gy					=	(float*)calloc(n_bins_f, sizeof(float));
	float*	f_dd_Gy					=	(float*)calloc(n_bins_f, sizeof(float));
	float*	f						=	(float*)calloc(n_bins_f, sizeof(float));
	float*	fdd						=	(float*)calloc(n_bins_f, sizeof(float));
	float*	dfdd					=	(float*)calloc(n_bins_f, sizeof(float));
	float	f0						=	0.0f;
	float	d_check					=	0.0f;

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

	SGP_SuccessiveConvolutions(	&f_parameters[0],		// u
			&n_bins_f,
			N2,
			// input + return values
			&n_bins_f1,
			f_d_Gy,
			f_dd_Gy,
			f,
			// return values
			&f0,
			fdd,
			dfdd,
			&d_check,
			&write_output_b,
			&shrink_tails_b,
			shrink_tails_under,
			&adjust_N2_b);

	long			n_bins_f_used	= n_bins_f1;

	float*	S						=	(float*)calloc(n_bins_f_used, sizeof(float));
	float	S_HCP, S_gamma, efficiency;

	SGP_get_gamma_response(	&n_bins_f_used,
							f_d_Gy,
							f_dd_Gy,
							f,
							&f0,
							gamma_model,
							gamma_parameters,
							// return
							S,
							&S_HCP,
							&S_gamma,
							&efficiency);

	results[0]			=	efficiency;				// 0 - 4: algo independent results
	results[1]			=	d_check;
	results[2]			=	S_HCP;
	results[3]			=	S_gamma;
	results[5]			=	f_parameters[0];		// 5 - 9: algo specific: u
	results[6]			=	u_start;
	results[7]			=	n_convolutions;

	free(f1_parameters);
	free(f_parameters);
	free(norm_fluence);
	free(dose_contribution_Gy);
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

void SGP_efficiency_grid(	long*	n,
							float*	E_MeV_u,
							long*	particle_no,
							float*	fluence_cm2,
							long*	material_no,
							long*	RDD_model,
							float*	RDD_parameters,
							long*	ER_model,
							float*	ER_parameters,
							long*	gamma_model,
							float*	gamma_parameters,
	//						method					= "grid",
							long*	N_runs,
							float*	fluence_factor,
							bool*	write_output,
							long*	nX,
							float*	grid_size_m,
							float*	results)
{
	FILE*		output_file;
	long 		i, j, k, m;
	long		n_grid					= (*nX) * (*nX);
	float		calc_grid_size_m		= (*grid_size_m) * (*nX);
	float		calc_grid_area_cm2		= calc_grid_size_m * calc_grid_size_m * 10000;

	output_file		=	fopen("GridSummation.log","w");
	if (output_file == NULL) return;											// File error

	fprintf(output_file, "##############################################################\n");
	fprintf(output_file, "##############################################################\n");
	fprintf(output_file, "This is SGP efficiency grid, version(2009/06/30).\n");
	fprintf(output_file, "##############################################################\n");
	fprintf(output_file, "\n\n\n");
	fprintf(output_file, "\nSystem time:");
	fprintf(output_file, "\n");
	fprintf(output_file, "calc grid:			%d*%d = %d pixels\n", *nX, *nX, n_grid);
	fprintf(output_file, "calc grid size/m:	 	%e\n", calc_grid_size_m);
	fprintf(output_file, "calc grid area/cm2:	%e\n", calc_grid_area_cm2);

	// Alloc checkerboard arrays
	float*		grid_d_Gy		= (float*)calloc(n_grid, sizeof(float));
	float*		grid_S			= (float*)calloc(n_grid, sizeof(float));

	// Clear results
	for (i = 0; i < 10; i++){
		results[i]		= 0.0;
	}

	// Get f1, f parameters
	float*		f1_parameters		= (float*)calloc(*n * 9, sizeof(float));
	float*		f_parameters		= (float*)calloc(7, sizeof(float));
	float*		norm_fluence		= (float*)calloc(*n, sizeof(float));
	float*		dose_contribution_Gy= (float*)calloc(*n, sizeof(float));
	float		max_r_max_m			= 0.0f;

	fprintf(output_file, "f1 parameters for %d particles\n", *n);
	for (i = 0; i < *n; i++){
		SGP_RDD_f1_parameters(	&E_MeV_u[i],
								&particle_no[i],
								material_no,
								RDD_model,
								RDD_parameters,
								ER_model,
								ER_parameters,
								&f1_parameters[i*9]);
		fprintf(output_file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 	f1_parameters[i*9 + 0],
																		f1_parameters[i*9 + 1],
																		f1_parameters[i*9 + 2],
																		f1_parameters[i*9 + 3],
																		f1_parameters[i*9 + 4],
																		f1_parameters[i*9 + 5],
																		f1_parameters[i*9 + 6],
																		f1_parameters[i*9 + 7],
																		f1_parameters[i*9 + 8]);
		max_r_max_m				= 	FMAX(max_r_max_m, f1_parameters[i*9 + 2]);
	}

	fprintf(output_file, "\nOverall r.max/m = %e\n\n", 	max_r_max_m);

	long		n_bins_f1			= 0;
	long		N2					= 0;
	SGP_SC_get_f1(	n,						// for f parameters only
					E_MeV_u,
					particle_no,
					fluence_cm2,
					material_no,
					RDD_model,
					RDD_parameters,
					ER_model,
					ER_parameters,
					&N2,
					&n_bins_f1,
					f1_parameters,
					norm_fluence,
					dose_contribution_Gy,
					f_parameters,
					NULL,
					NULL,
					NULL);

	fprintf(output_file, "f parameters\n");
	fprintf(output_file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 	f_parameters[0],
															f_parameters[1],
															f_parameters[2],
															f_parameters[3],
															f_parameters[4],
															f_parameters[5],
															f_parameters[6]);

	// Largest r.max --> calculate size of sample area
	float sample_grid_size_m	= calc_grid_size_m + 2.01f * max_r_max_m;
	float sample_grid_area_cm2	= sample_grid_size_m * sample_grid_size_m * 10000;
	fprintf(output_file, "sample grid size/m   = %e\n", sample_grid_size_m);
	fprintf(output_file, "sample grid area/cm2 = %e\n", sample_grid_area_cm2);

	// mean and actual number of particles on sample_area
	float*	mean_number_particles	= (float*)calloc(*n, sizeof(float));
	long*	act_number_particles	= (long*)calloc(*n, sizeof(float));
	for (i = 0; i < *n; i++){
		mean_number_particles[i]	= sample_grid_area_cm2 * f_parameters[1] * norm_fluence[i];				// Area * Total_fluence (particle i)
	}

	// create and initialise RNGs
	gsl_rng * rng1 	= gsl_rng_alloc (gsl_rng_taus);
	gsl_rng * rng2 	= gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(rng1, 12345678);
	gsl_rng_set(rng2, 2345678);


	// run loop
	long	n_particles;
	for (m = 0; m < *N_runs; m++){
		float*	run_results			= (float*)calloc(10, sizeof(float));

		// sample particles numbers
		fprintf(output_file, "Actual number of particles (mean)\n");
		n_particles	= 0;
		for (i = 0; i < *n; i++){
			act_number_particles[i]	= 	(long)gsl_ran_poisson(rng1, mean_number_particles[i]);
			n_particles				+=	act_number_particles[i];
			fprintf(output_file, "particle %d: %d (%e)\n", i, act_number_particles[i], mean_number_particles[i]);
		}
		fprintf(output_file, "\nIn total: %d\n", n_particles);

		// alloc particle array
		float*	x_pos				= (float*)calloc(n_particles, sizeof(float));
		float*	y_pos				= (float*)calloc(n_particles, sizeof(float));
		long*	particle_index		= (long*)calloc(n_particles, sizeof(long));
		float*	r_m					= (float*)calloc(n_particles, sizeof(float));
		float*	r_max_m				= (float*)calloc(n_particles, sizeof(float));
		long	n_tmp				= 1;
		float	d_tmp_Gy			= 0.0;

		// fill in index / r_max_m
		j		= 0;
		k		= 0;
		for (i = 0; i < n_particles; i++){
			if(k >= act_number_particles[j]){
				k 		= 0;
				j++;
			}
			particle_index[i]	=	j;
			r_max_m[i]			=	f1_parameters[j*9 + 2];
		}

		// sample particle positions
		for (i = 0; i < n_particles; i++){
			x_pos[i]					= (float)gsl_rng_uniform_pos(rng2) * sample_grid_size_m;
			y_pos[i]					= (float)gsl_rng_uniform_pos(rng2) * sample_grid_size_m;
		}

		// grid loop
		for (j = 0; j < *nX; j++){					// y
			float cur_y_pos			=	max_r_max_m + ((float)j + 0.5f)*(*grid_size_m);
			for (i = 0; i < *nX; i++){				// x
				float cur_x_pos			=	max_r_max_m + ((float)i + 0.5f)*(*grid_size_m);
				grid_d_Gy[j * (*nX) + i]=	0.0f;
				for (k = 0; k < n_particles; k++){	// particles
					r_m[k]					=	sqrt( (x_pos[k] - cur_x_pos) * (x_pos[k] - cur_x_pos) +
													  (y_pos[k] - cur_y_pos) * (y_pos[k] - cur_y_pos));
					if(r_m[k] <= r_max_m[k]){		// does particle contribute?
						SGP_D_RDD_Gy(	&n_tmp,
										&r_m[k],
										&E_MeV_u[particle_index[k]],
										&particle_no[particle_index[k]],
										material_no,
										RDD_model,
										RDD_parameters,
										ER_model,
										ER_parameters,
										&d_tmp_Gy);
						grid_d_Gy[j * (*nX) + i]	+=	d_tmp_Gy;
					} // particle contribution
				}// particle loop
			} // x loop
		} // y loop

		// get gamma response for local dose
		SGP_gamma_response(	&n_grid,
							grid_d_Gy,
							gamma_model,
							gamma_parameters,
							grid_S);
		// averaging
		float d_total_Gy 	= 0.0f;
		float S_HCP			= 0.0f;
		for (i = 0; i < n_grid; i++){
			d_total_Gy		+=	grid_d_Gy[i];
			S_HCP			+=	grid_S[i];
		}

		d_total_Gy		/= n_grid;
		float S_gamma	= 0.0f;
		SGP_gamma_response(	&n_tmp,
							grid_d_Gy,
							gamma_model,
							gamma_parameters,
							&S_gamma);
		S_gamma				*= n_grid;

		float efficiency	= 0.0f;
		if(S_gamma > 0){
			efficiency = S_HCP / S_gamma;
		}

		// write graph (first run)
		bool write_graph = true;
			if(write_graph & m == 1){
				FILE*		graph_file;
				graph_file		=	fopen("GridGraph.csv","w");
				if (graph_file == NULL) return;		// File error

				fprintf(graph_file, "x.m; y.m; d.Gy; S\n");

				for (j = 0; j < *nX; j++){
					for (i = 0; i < *nX; i++){
						fprintf(graph_file, "%e; %e; %e; %e\n",	max_r_max_m + ((float)i + 0.5f)*(*grid_size_m),
																max_r_max_m + ((float)j + 0.5f)*(*grid_size_m),
																grid_d_Gy[j * (*nX) + i],
																grid_S[j * (*nX) + i]);
					}
				}
			}

		run_results[0]		= efficiency;
		run_results[1]		= d_total_Gy;
		run_results[2]		= S_HCP;
		run_results[3]		= S_gamma;
		run_results[4]		= n_particles;

		// copy to results
		results[0]			+= run_results[0];
		results[1]			+= run_results[1];
		results[2]			+= run_results[2];
		results[3]			+= run_results[3];
		results[4]			+= run_results[4];

		results[5]			+= run_results[0]*run_results[0];
		results[6]			+= run_results[1]*run_results[1];
		results[7]			+= run_results[2]*run_results[2];
		results[8]			+= run_results[3]*run_results[3];
		results[9]			+= n_particles * n_particles;

		fprintf(output_file, "\n\nRun %d results\n", m + 1);
		fprintf(output_file, "efficiency 		= %e\n", run_results[0]);
		fprintf(output_file, "d.check.Gy 		= %e\n", run_results[1]);
		fprintf(output_file, "S (HCP) 			= %e\n", run_results[2]);
		fprintf(output_file, "S (gamma) 		= %e\n", run_results[3]);
		fprintf(output_file, "no. particles		= %e\n", run_results[4]);

		free(x_pos);
		free(y_pos);
		free(particle_index);
		free(r_max_m);
	}// end run loop

	results[0]	/= *N_runs;
	results[1]	/= *N_runs;
	results[2]	/= *N_runs;
	results[3]	/= *N_runs;
	results[4]	/= *N_runs;

	results[5]	/= *N_runs;
	results[6]	/= *N_runs;
	results[7]	/= *N_runs;
	results[8]	/= *N_runs;
	results[9]	/= *N_runs;

	results[5]	-= results[0]*results[0];
	results[6]	-= results[1]*results[1];
	results[7]	-= results[2]*results[2];
	results[8]	-= results[3]*results[3];
	results[9]	-= results[4]*results[4];

	results[5]	= FMAX(0, results[5]);
	results[6]	= FMAX(0, results[6]);
	results[7]	= FMAX(0, results[7]);
	results[8]	= FMAX(0, results[8]);
	results[9]	= FMAX(0, results[9]);

    results[5]	= sqrt(results[5] / (*N_runs - 1));
	results[6]	= sqrt(results[6] / (*N_runs - 1));
	results[7]	= sqrt(results[7] / (*N_runs - 1));
	results[8]	= sqrt(results[8] / (*N_runs - 1));
	results[9]	= sqrt(results[9] / (*N_runs - 1));

	fprintf(output_file, "\n###############################################\nResults\n");
	fprintf(output_file, "efficiency 		= %e +/- %e\n", results[0], results[5]);
	fprintf(output_file, "d.check.Gy 		= %e +/- %e\n", results[1], results[6]);
	fprintf(output_file, "S (HCP) 			= %e +/- %e\n", results[2], results[7]);
	fprintf(output_file, "S (gamma) 		= %e +/- %e\n", results[3], results[8]);
	fprintf(output_file, "no. particles		= %e +/- %e\n", results[4], results[9]);
	fprintf(output_file, "###############################################\n");

	gsl_rng_free(rng1);
	gsl_rng_free(rng2);

	free(f1_parameters);
	free(f_parameters);
	free(norm_fluence);
	free(dose_contribution_Gy);

	free(mean_number_particles);
	free(act_number_particles);

	free(grid_d_Gy);
	free(grid_S);

	close(output_file);
}
/*
BOOL APIENTRY DllMain( HANDLE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
    return TRUE;
}
*/

