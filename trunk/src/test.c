/*
 * test.c
 *
 *  Created on: 2009-06-08
 *      Author: grzanka
 */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

//#include "SGP_Constants.h"
//#include "SGP_SuccessiveConvolutions.h"

void SGP_SC_get_f1_array_sizeS(	long*	n,
		float*	E_MeV_u,
		long*	particle_no,
		char**	material_name,
		float*	parameter,
		long*	N2,
		// from here: return values
		long*	n_bins_f1,
		bool*	debug);

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
		long*	n_convolutions);

void SGP_SC_get_f1S(	long*	n,
		float*	E_MeV_u,
		long*	particle_no,
		float*	fluence_cm2,
		char**	material_name,
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

void testRDD(){

    long n = 5;
    float r_m[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
    float E_MeV_u[] = {60, 60, 60, 60, 60};
    long particle_no[] = {1, 1, 1, 1, 1};
    char material_name[50] = "Water, Liquid";
    char ** mn = (char**)calloc(1, sizeof(char*));
    *mn = material_name;
    float parameter = 1e-8;
	float D_RDD_Gy[] = { 0, 0, 0, 0, 0 };
	int i;

	printf("begin %s\n", *mn);

//	SGP_D_RDD_Gy( &n, r_m, &E_MeV_u, &particle_no, material_name, &parameter, &D_RDD_Gy);
	SGP_D_RDD_GyS( &n, r_m, &E_MeV_u, &particle_no, mn, &parameter, &D_RDD_Gy);

	free(mn);

	for( i = 0 ; i < n ; i++){
		printf("end, D_RRD_Gy[%g] = %g\n", r_m[i], D_RDD_Gy[i]);
	}


}


void test_SGP_SC_1(){

	/****************** STEP 1 **************/

	// INPUT :
	long n = 1;
	float E_MeV_u[] = {60.};
	float fluence_cm2[] = {1e8};
	long particle_index[] = {1};
	char material_name[50] = "Water, Liquid";
	char ** mn = (char**)calloc(1,sizeof(char*));
	*mn = material_name;
	float rdd_parameter[] = {1e-8};
	long N2 = 10;
	bool debug = false;

	// OUTPUT :
	long n_bins_f1;

	SGP_SC_get_f1_array_sizeS(&n,E_MeV_u,particle_index,mn,rdd_parameter,&N2,\
			&n_bins_f1,\
			&debug);

	printf("SGP_SC_get_f1_array_sizeS finished:\n");
	printf(" - n_bins_f1 = %ld\n", n_bins_f1);

	/****************** STEP 2 **************/

	// OUTPUT :
	float*	norm_fluence = (float*)calloc(n,sizeof(float));
	float*	LET_MeV_cm2_g = (float*)calloc(n,sizeof(float));
	float*	r_min_m = (float*)calloc(n,sizeof(float));
	float*	r_max_m = (float*)calloc(n,sizeof(float));
	float*	d_min_Gy = (float*)calloc(n,sizeof(float));
	float*	d_max_Gy = (float*)calloc(n,sizeof(float));
	float*	k_Gy = (float*)calloc(n,sizeof(float));
	float*	single_impact_fluence_cm2 = (float*)calloc(n,sizeof(float));
	float*	single_impact_dose_Gy = (float*)calloc(n,sizeof(float));
	float*	dose_contribution_Gy = (float*)calloc(n,sizeof(float));
	float	total_fluence_cm2;
	float	total_dose_Gy;
	float	ave_E_MeV;
	float	dw_E_MeV;
	float	ave_LET_MeV_cm2_g;
	float	dw_LET_MeV_cm2_g;
	float	u;
	float*	f1_d_Gy = (float*)calloc(n_bins_f1,sizeof(float));
	float*	f1_dd_Gy = (float*)calloc(n_bins_f1,sizeof(float));
	float*	f1 = (float*)calloc(n_bins_f1,sizeof(float));


	SGP_SC_get_f1S(&n,E_MeV_u,particle_index,fluence_cm2,mn,rdd_parameter,&N2,&n_bins_f1,\
			norm_fluence,LET_MeV_cm2_g,r_min_m,r_max_m,d_min_Gy,d_max_Gy,k_Gy,\
			single_impact_fluence_cm2, single_impact_dose_Gy,dose_contribution_Gy,\
			&total_fluence_cm2,&total_dose_Gy,&ave_E_MeV,&dw_E_MeV,&ave_LET_MeV_cm2_g,\
			&dw_LET_MeV_cm2_g,&u,f1_d_Gy,f1_dd_Gy,f1,\
			&debug);

	printf("SGP_SC_get_f1S finished:\n");
	printf(" - norm_fluence[0] = %g\n", norm_fluence[0]);
	printf(" - LET_MeV_cm2_g[0] = %g\n", LET_MeV_cm2_g[0]);
	printf(" - r_min_m[0] = %g\n", r_min_m[0]);
	printf(" - r_max_m[0] = %g\n", r_max_m[0]);
	printf(" - d_min_Gy[0] = %g\n", d_min_Gy[0]);
	printf(" - d_max_Gy[0] = %g\n", d_max_Gy[0]);
	printf(" - k_Gy[0] = %g\n", k_Gy[0]);
	printf(" - single_impact_fluence_cm2[0] = %g\n", single_impact_fluence_cm2[0]);
	printf(" - single_impact_dose_Gy[0] = %g\n", single_impact_dose_Gy[0]);
	printf(" - dose_contribution_Gy[0] = %g\n", dose_contribution_Gy[0]);
	printf(" - total_fluence_cm2 = %g\n", total_fluence_cm2);
	printf(" - total_dose_Gy = %g\n", total_dose_Gy);
	printf(" - ave_E_MeV = %g\n", ave_E_MeV);
	printf(" - dw_E_MeV = %g\n", dw_E_MeV);
	printf(" - ave_LET_MeV_cm2_g = %g\n", ave_LET_MeV_cm2_g);
	printf(" - dw_LET_MeV_cm2_g = %g\n", dw_LET_MeV_cm2_g);
	printf(" - u = %g\n", u);
	printf(" - f1_d_Gy[] = %g,%g,%g,...\n", f1_d_Gy[0], f1_d_Gy[1], f1_d_Gy[2]);
	printf(" - f1_dd_Gy[] = %g,%g,%g,...\n", f1_dd_Gy[0], f1_dd_Gy[1], f1_dd_Gy[2]);
	printf(" - f1[] = %g,%g,%g,...\n", f1[0], f1[1], f1[2]);

	/****************** STEP 3 **************/

	// INPUT :
	float fluence_factor = 1.0;
	float u_1 = u;
	// OUTPUT :
	long n_bins_f;
	float u_start;
	long n_convolutions;

	SGP_SC_get_f_array_size(\
			&u_1,&fluence_factor,&N2,&n_bins_f1,\
			f1_d_Gy,f1_dd_Gy,f1,\
			&n_bins_f,&u_start,&n_convolutions);

	printf("SGP_SC_get_f_array_size: finished\n");
	printf(" - n_bins_f = %ld\n", n_bins_f);
	printf(" - u_start = %g\n", u_start);
	printf(" - n_convolutions = %ld\n", n_convolutions);


	/****************** STEP 4 **************/

	// OUTPUT:
	float * f_d_Gy = (float*)calloc(n_bins_f,sizeof(float));
	float * f_dd_Gy = (float*)calloc(n_bins_f,sizeof(float));
	float * f_start = (float*)calloc(n_bins_f,sizeof(float));

	SGP_SC_get_f_start(\
			&u_start,&n_bins_f1,&N2,f1_d_Gy,f1_dd_Gy,f1,&n_bins_f,\
			f_d_Gy,f_dd_Gy,f_start);

	printf("SGP_SC_get_f_start: finished\n");
	printf(" - f_d_Gy[] = %g,%g,%g,...\n", f_d_Gy[0], f_d_Gy[1], f_d_Gy[2]);
	printf(" - f_dd_Gy[] = %g,%g,%g,...\n", f_dd_Gy[0], f_dd_Gy[1], f_dd_Gy[2]);
	printf(" - f_start[] = %g,%g,%g,...\n", f_start[0], f_start[1], f_start[2]);

	/****************** STEP 5 **************/

	// INPUT:
	bool write_output = true;
	bool shrink_tails = false;
	float shrink_tails_under = 1e-30;
	bool adjust_N2 = true;

	long n_bins_f_used = n_bins_f1;
	// OUTPUT:
	float f0 = 0.0;
	float * fdd = (float*)calloc(n_bins_f,sizeof(float));
	float * dfdd = (float*)calloc(n_bins_f,sizeof(float));
	float d = 0.0;

	SGP_SuccessiveConvolutions(&u,&n_bins_f,&N2,&n_bins_f_used,\
			f_d_Gy,f_dd_Gy,f_start,\
			&f0,fdd,dfdd,&d,\
			&write_output,&shrink_tails,&shrink_tails_under,&adjust_N2);

	printf("SGP_SuccessiveConvolutions: finished\n");
	printf(" - f0 = %g\n", f0);
	printf(" - fdd[] = %g,%g,%g,...\n", fdd[0], fdd[1], fdd[2]);
	printf(" - dfdd[] = %g,%g,%g,...\n", dfdd[0], dfdd[1], dfdd[2]);
	printf(" * n_bins_f_used = %ld\n", n_bins_f_used);
	printf(" * f_d_Gy[] = %g,%g,%g,...\n", f_d_Gy[0], f_d_Gy[1], f_d_Gy[2]);
	printf(" * f_dd_Gy[] = %g,%g,%g,...\n", f_dd_Gy[0], f_dd_Gy[1], f_dd_Gy[2]);
	printf(" * f_start[] = %g,%g,%g,...\n", f_start[0], f_start[1], f_start[2]);

	free(norm_fluence);
	free(LET_MeV_cm2_g);
	free(r_min_m);
	free(r_max_m);
	free(d_min_Gy);
	free(d_max_Gy);
	free(k_Gy);
	free(single_impact_fluence_cm2);
	free(single_impact_dose_Gy);
	free(dose_contribution_Gy);
	free(f1_d_Gy);
	free(f1_dd_Gy);
	free(f1);
	free(f_d_Gy);
	free(f_dd_Gy);
	free(f_start);
	free(fdd);
	free(dfdd);
}

void test_SGP_SC_2(){

	// INPUT :
	long n = 1;
	float E_MeV_u[] = {60.};
	long particle_no[] = {1};
	float fluence_cm2[] = {1e8};
	long slab_no[] = {1};
	char material_name[50] = "Water, Liquid";
	//char ** mn = (char**)calloc(1,sizeof(char*));
	//*mn = material_name;
	float rdd_parameter[] = {1e-8};
	long N2 = 1;
	long n_slabs = 1;
	long n_gamma_parameter = 2;
	long gamma_model = 0;
	float gamma_parameter[] = {1., 1.};
	long verbose_level = 2;
	char output_fileName[50] = "SC.dat";
	float*	u = (float*)calloc(n_slabs,sizeof(float));
	float*	total_d_Gy = (float*)calloc(n_slabs,sizeof(float));
	float*	d = (float*)calloc(n_slabs,sizeof(float));
	float*	S_HCP = (float*)calloc(n_slabs,sizeof(float));
	float*	S_gamma = (float*)calloc(n_slabs,sizeof(float));
	float*	efficiency = (float*)calloc(n_slabs,sizeof(float));
	float*	S_HCP_total = (float*)calloc(n_slabs,sizeof(float));
	float*	S_gamma_total = (float*)calloc(n_slabs,sizeof(float));
	float	efficiency_total;

	SGP_SC_Loop(&n, E_MeV_u, particle_no, fluence_cm2, slab_no, material_name,\
			rdd_parameter,&N2,&n_slabs,&n_gamma_parameter,\
			&gamma_model,gamma_parameter,&verbose_level,output_fileName,\
			u,total_d_Gy,d,S_HCP,S_gamma,efficiency,S_HCP_total,\
			S_gamma_total,&efficiency_total);

	free(u);
	free(total_d_Gy);
	free(d);
	free(S_HCP);
	free(S_gamma);
	free(efficiency);
	free(S_HCP_total);
	free(S_gamma_total);
}

int main(){

	//testRDD();
	test_SGP_SC_1();
	//test_SGP_SC_2();

	return 1;
};
