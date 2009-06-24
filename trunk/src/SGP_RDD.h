#ifndef SGP_RDD_H_
#define SGP_RDD_H_

#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>


void SGP_RDD_f1_parameters(	/* radiation field parameters */
							float*	E_MeV_u,
							long*	particle_no,
							/* detector parameters */
							long*	material_no,
							/* radial dose distribution model */
							long*	rdd_model,
							long*	n_rdd_parameter,
							float*	rdd_parameter,
							/* electron range model */
							long*	er_model,
							long*	n_er_parameter,
							float*	er_parameter,
							/* calculated parameters */
							long* 	n_f1_parameters,
							float* 	f1_parameters);

void SGP_D_RDD_Gy	(	long*	n,
		float*	r_m,
		/* radiation field parameters */
		float*	E_MeV_u,
		long*	particle_no,
		/* detector parameters */
		long*	material_no,
		/* radial dose distribution model */
		long*	rdd_model,       /* */
		long*	n_rdd_parameter, /* number of rdd parameters */
		float*	rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		float*	D_RDD_Gy);

void SGP_r_RDD_m	(	long*	n,
		float*	D_RDD_Gy,
		/* radiation field parameters */
		float*	E_MeV_u,
		long*	particle_no,
		/* detector parameters */
		long*	material_no,
		/* radial dose distribution model */
		long*	rdd_model,       /* */
		long*	n_rdd_parameter, /* number of rdd parameters */
		float*	rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		float*	r_m);


inline float SGP_RDD_Katz_point_kernel(float* x, float* alpha){
	return (1.0f/((*x)*(*x)) )*pow(1.0f - (*x), 1.0f / (*alpha));
}

void SGP_RDD_Katz_point_kernelS(int *n, float* x, float* alpha, float* f){
	int i;
	for( i = 0 ; i < *n ; i++){
		f[i] = SGP_RDD_Katz_point_kernel(&(x[i]),alpha);
	};
}

inline float SGP_RDD_Katz_point_coeff_Gy(float* C_J_m,float* Z_eff, float* beta, float* alpha, float* density_kg_m3, float* r_max_m){
	return (*C_J_m) * (*Z_eff)*(*Z_eff) / (2.0f * M_PI * (*beta)*(*beta) * (*alpha) * (*density_kg_m3) * (*r_max_m));
}

inline float SGP_RDD_Katz_point_Gy(float* r_m, float* alpha, float* r_max_m, float* Katz_point_coeff_Gy){
	float x = (*r_m)/(*r_max_m);
	return (*Katz_point_coeff_Gy) * SGP_RDD_Katz_point_kernel(&x,alpha);
}

void SGP_RDD_Katz_point_GyS(int *n, float* r_m, float* alpha, float* r_max_m,float* Katz_point_coeff_Gy, float * D){
	int i;
	for( i = 0 ; i < *n ; i++){
		D[i] = SGP_RDD_Katz_point_Gy(&(r_m[i]),alpha,r_max_m,Katz_point_coeff_Gy);
	};
}

inline float SGP_RDD_Katz_dEdx_kernel(float* x, float* alpha){
	return (1.0f/(*x) )*pow(1.0f - (*x), 1.0f / (*alpha));
}

void SGP_RDD_Katz_dEdx_kernelS(int *n, float* x, float* alpha, float* f){
	int i;
	for( i = 0 ; i < *n ; i++){
		f[i] = SGP_RDD_Katz_dEdx_kernel(&(x[i]),alpha);
	};
}

double SGP_RDD_Katz_dEdx_integrand(double x, void * params){
	float alpha = ((float*)params)[0];
	float f_x = (float)(x);
	return (double)SGP_RDD_Katz_dEdx_kernel( &f_x, &alpha);
}

inline float SGP_RDD_Katz_dEdx_coeff_J_m(float* r_max_m, float* density_kg_m3, float* Katz_point_coeff_Gy){
	return 2 * M_PI * (*density_kg_m3) * (*r_max_m)*(*r_max_m) * (*Katz_point_coeff_Gy);
}

float SGP_RDD_Katz_dEdx_J_m(float* alpha, float* r_min_m, float* r_max_m, float* Katz_dEdx_coeff_J_m){
// integration ...
	double dEdx_integral = 1.0;
	double error;
	gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &SGP_RDD_Katz_dEdx_integrand;
	float params[] = {*alpha};
	F.params = params;
	int status = gsl_integration_qags (&F, (*r_min_m)/(*r_max_m), 1.0, 0, 1e-5, 1000, w1, &dEdx_integral, &error);
	if (status == GSL_EROUND || status == GSL_ESING){
#ifdef _DEBUG
		indnt_inc();
		fprintf(debf,"%s status == %d\n",isp,status);
#endif
		dEdx_integral = -1.0f;
	}

	gsl_integration_workspace_free (w1);

	return (*Katz_dEdx_coeff_J_m)*(float)dEdx_integral;
}

inline float SGP_RDD_Katz_site_Gy(float* r_m, float* alpha, float* r_min_m, float* r_max_m, float* LET_J_m, float* density_kg_m3, float* Katz_dEdx_J_m, float* Katz_point_coeff_Gy){
	if( (*r_m) < (*r_min_m) ){
		return (1.0f / ((*density_kg_m3) * M_PI * (*r_min_m)*(*r_min_m)))*((*LET_J_m) - (*Katz_dEdx_J_m));
	} else {
		return SGP_RDD_Katz_point_Gy(r_m, alpha, r_max_m, Katz_point_coeff_Gy);
	}
	return 0.0f;
}

void SGP_RDD_Katz_site_GyS(int *n, float* r_m, float* alpha, float* r_min_m, float* r_max_m, float* LET_J_m, float* density_kg_m3, float* Katz_dEdx_J_m, float* Katz_point_coeff_Gy, float * D_Gy){
	int i;
	for( i = 0 ; i < *n ; i++){
		D_Gy[i] = SGP_RDD_Katz_site_Gy(&(r_m[i]),alpha,r_min_m,r_max_m,LET_J_m,density_kg_m3,Katz_dEdx_J_m,Katz_point_coeff_Gy);
	};
}

float
geometryFunctionPhi (float* r0_m, float* a0_m, float* r_m)
{
//#ifdef _DEBUG
//		indnt_init();
//		fprintf(debf,"%sr0_m = %g , a0_m = %g , r_m = %g\n",isp,*r0_m,*a0_m, *r_m);
//#endif
	float res = 0.;
	double factor = 0.;
	gsl_complex carg, cres;
	if ((*r_m) <= fabs ((*r0_m) - (*a0_m)))
	{
		if ((*r0_m) >= (*a0_m))
			res = 0.0f;
		else
			res = (float)M_PI;
//#ifdef _DEBUG
//		indnt_init();
//		fprintf(debf,"%sA res = %g\n",isp,res);
//#endif
	}
	else
	{
		factor = gsl_pow_2 (*a0_m) - gsl_pow_2 ((*r0_m) - (*r_m));
		factor /= gsl_pow_2 ((*r_m) + (*r0_m)) - gsl_pow_2 (*a0_m);
//#ifdef _DEBUG
//		indnt_init();
//		fprintf(debf,"%sB factor = %g\n",isp,factor);
//#endif
		GSL_SET_COMPLEX (&carg, sqrt (factor), 0.);
		cres = gsl_complex_arctan (carg);
		res = 2.0f * (float)(GSL_REAL (cres));
//#ifdef _DEBUG
//		indnt_init();
//		fprintf(debf,"%sB res = %g\n",isp,res);
//#endif
	}
	return res;
}

inline float SGP_RDD_Katz_ext_kernel_Gy(float* t_m, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy){
//#ifdef _DEBUG
//		indnt_init();
//		fprintf(debf,"%st_m = %g , a0_m = %g , r_m = %g\n",isp,*t_m,*a0_m, *r_m);
//		fprintf(debf,"%sf1 = %g\n",isp,(1.0f/ (M_PI * (*a0_m)*(*a0_m))));
//		fprintf(debf,"%sf2 = %g\n",isp,SGP_RDD_Katz_point_Gy(t_m,alpha,r_max_m,Katz_point_coeff_Gy));
//		fprintf(debf,"%sf3 = %g\n",isp,geometryFunctionPhi(r_m,a0_m,t_m));
//#endif
	if( (*t_m) < (*r_min_m ) )
		return 0.0;
	if( (*t_m) >= (*a0_m) + (*r_m) )
		return 0.0;
	if( ((*r_m) >= (*a0_m)) && ((*t_m) <= (*r_m) - (*a0_m)))
		return 0.0;
	else
		return  (1.0f/ (M_PI * (*a0_m)*(*a0_m))) * SGP_RDD_Katz_point_Gy(t_m,alpha,r_max_m,Katz_point_coeff_Gy) *  geometryFunctionPhi(r_m,a0_m,t_m) * (*t_m);
}

void SGP_RDD_Katz_ext_kernel_GyS(int *n, float* t_m, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy){
	int i;
	for( i = 0 ; i < *n ; i++){
//#ifdef _DEBUG
//		indnt_init();
//		fprintf(debf,"%st_m[%d] = %g\n",isp,i,t_m[i]);
//#endif
		D_Gy[i] = SGP_RDD_Katz_ext_kernel_Gy(&(t_m[i]),r_m,a0_m,alpha,r_min_m,r_max_m,Katz_point_coeff_Gy);
	};
}

double SGP_RDD_Katz_ext_integrand_Gy(double t_m, void * params){
	float r_m = ((float*)params)[0];
	float a0_m = ((float*)params)[1];
	float alpha = ((float*)params)[2];
	float r_min_m = ((float*)params)[3];
    float r_max_m = ((float*)params)[4];
    float Katz_point_coeff_Gy = ((float*)params)[5];
	float f_t_m = (float)(t_m);
	return (double)SGP_RDD_Katz_ext_kernel_Gy(&f_t_m, &r_m, &a0_m, &alpha, &r_min_m, &r_max_m, &Katz_point_coeff_Gy);
}

inline float SGP_RDD_Katz_ext_Gy(float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy){
	double int_lim_m = 0.0f;
	if( (*r_m) > (*a0_m) ){
		int_lim_m = GSL_MAX((*r_m) - (*a0_m),*r_min_m);
	}
    gsl_set_error_handler_off();

	double ext_integral_Gy;
	double error;
	gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &SGP_RDD_Katz_ext_integrand_Gy;
	float params[] = {*r_m,*a0_m,*alpha,*r_min_m,*r_max_m,*Katz_point_coeff_Gy};
	F.params = params;
	int status = gsl_integration_qags (&F, int_lim_m, (*r_m)+(*a0_m), 0, 1e-5, 1000, w1, &ext_integral_Gy, &error);
	if (status == GSL_EROUND || status == GSL_ESING){
#ifdef _DEBUG
		indnt_init();
		fprintf(debf,"%s status == %d\n",isp,status);
#endif
		ext_integral_Gy = -1.0f;
	}
	gsl_integration_workspace_free (w1);

	return ext_integral_Gy;
}

void SGP_RDD_Katz_ext_GyS(int *n, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy){
	int i;
	for( i = 0 ; i < *n ; i++){
		D_Gy[i] = SGP_RDD_Katz_ext_Gy(&(r_m[i]), a0_m, alpha, r_min_m, r_max_m, Katz_point_coeff_Gy);
	};
}

/*	f1_parameters:
 * 		0 - LET_MeV_cm2_g
 * 		1 - r_min_m
 * 		2 - r_max_m
 * 		3 - d_min_Gy
 * 		4 - d_max_Gy
 * 		5 - k 						(norm. constant)
 * 		6 - single_impact_fluence_cm2
 * 		7 - single_impact_dose_Gy
 * 		8 - dEdx_MeV_cm2_g
 */
void SGP_RDD_f1_parameters(	/* radiation field parameters */
							float*	E_MeV_u,
							long*	particle_no,
							/* detector parameters */
							long*	material_no,
							/* radial dose distribution model */
							long*	rdd_model,
							long*	n_rdd_parameter,
							float*	rdd_parameter,
							/* electron range model */
							long*	er_model,
							long*	n_er_parameter,
							float*	er_parameter,
							/* calculated parameters */
							long * 	n_f1_parameters,
							float * f1_parameters)
{
#ifdef _DEBUG
	indnt_inc();
	fprintf(debf,"%sbegin SGP_RDD_f1_parameters\n",isp);
	fprintf(debf,"%sbegin RDD model = %ld\n",isp,*rdd_model);
	fprintf(debf,"%sMaterial = %ld\n", isp, *material_no);
#endif

	// Get beta, Z and Zeff
	long	n_tmp		= 1;
	float	beta		= 0.0f;
	long	Z			= 0;
	float	Z_eff		= 0.0f;
	SGP_beta_from_particle_no(	&n_tmp,
									E_MeV_u,
									particle_no,
									&beta);
	SGP_Z_from_particle_no(	&n_tmp,
								particle_no,
								&Z);
	SGP_effective_charge_from_beta(	&n_tmp,
									&beta,
									&Z,
									&Z_eff);

	// Get energy constant C == k
	float 	density_g_cm3, density_kg_m3, electron_density_m3, I_eV, alpha_g_cm2_MeV, p_MeV, m_g_cm2;
	SGP_getMaterialData(	&n_tmp,
							material_no,
							&density_g_cm3,
							&electron_density_m3,
							&I_eV,
							&alpha_g_cm2_MeV,
							&p_MeV,
							&m_g_cm2);
	density_kg_m3			=	density_g_cm3 * 1000.0f;

	///////////////////////////////////////////////////////////////
	// PARAMETER 0: Get the LET (same for all models)
	SGP_LET_MeV_cm2_g(	&n_tmp,
						E_MeV_u,
						particle_no,
						material_no,
						&f1_parameters[0]);

	/////////////////////////////////////////////////////////////////////////////
	// PARAMETER 1: Get the LET maximum electron range (same for all RDD models)
	SGP_max_electron_range_m(	&n_tmp,
								E_MeV_u,
								particle_no,
								material_no,
								er_model,
								&f1_parameters[2]);


	/////////////////////////////////////////////////////////////////////////////
	// MODEL SPEZIFIC PARAMETERS
	if( *rdd_model == RDD_Test){
		f1_parameters[1] 	= 0.0f;																								// r_min_m
		f1_parameters[6]	= 1.0f / (pi * f1_parameters[2]*f1_parameters[2] * m_to_cm * m_to_cm);								// pi * r_max_m^2 = Track area -> single_impact_fluence [1/cm2]
		f1_parameters[5]	= f1_parameters[6] * f1_parameters[0] * MeV_to_J * 1000.0f;											// LET  / track area = Norm.constant k
		f1_parameters[3]	= f1_parameters[5];																					// d_min_Gy = k
		f1_parameters[4]	= f1_parameters[5];																					// d_max_Gy = k
		f1_parameters[7]	= f1_parameters[0] * f1_parameters[5] * MeV_g_to_J_kg;												// single_impact_dose = LET * fluence;
		f1_parameters[8]	= f1_parameters[0];																					// dEdx = LET
	}

	if( *rdd_model == RDD_KatzPoint || *rdd_model == RDD_Site || *rdd_model == RDD_ExtTarget){
		//////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////
		float alpha				= 	1.667f;
		float w_el_keV			=	2.0f * electron_mass_MeV_c2 * 1000.0f * beta*beta / (1 - beta*beta);
		if(w_el_keV <= 1.0f){
			alpha					= 1.079f;
		}
		//////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////
		f1_parameters[1]		=	rdd_parameter[0];								// r_min_m

		float	N_el_cm3		=	electron_density_m3 / (100*100*100);
		float	C_J_cm			=	2.0f * pi * N_el_cm3 * (e_esu*e_esu*e_esu*e_esu) / (electron_mass_g * c_cm_s *c_cm_s) * 1e-7;   // energy constant [J/cm] not [erg/cm] hence 10^-7
		float	C_J_m			=	C_J_cm * 100.0f;
		f1_parameters[5]		= 	C_J_m;											// Norm.constant k


		// Get dEdx by simple integration from r_min_m (in case of RDD_Site = a0) to r_max_m
		float	dEdx_J_m		=	0.0f;
		float	LET_J_m 		=	f1_parameters[0] * 1.602e-13 * density_g_cm3 * 100.0f;

		float Katz_point_coeff_Gy = SGP_RDD_Katz_point_coeff_Gy(&C_J_m,&Z_eff,&beta,&alpha,&density_kg_m3,&(f1_parameters[2]));
		float Katz_dEdx_coeff_J_m = SGP_RDD_Katz_dEdx_coeff_J_m(&(f1_parameters[2]),&density_kg_m3,&Katz_point_coeff_Gy);
		dEdx_J_m = SGP_RDD_Katz_dEdx_J_m(&alpha,&(f1_parameters[1]),&(f1_parameters[2]),&Katz_dEdx_coeff_J_m);

		float	dEdx_MeV_g_cm2		=	dEdx_J_m / 100.0f / density_g_cm3 / MeV_to_J;

#ifdef _DEBUG
	fprintf(debf,"%sLET_J_m = %g\n",isp,LET_J_m);
	fprintf(debf,"%sKatz_point_coeff_Gy = %g\n",isp,Katz_point_coeff_Gy);
	fprintf(debf,"%sKatz_dEdx_coeff_J_m = %g\n",isp,Katz_dEdx_coeff_J_m);
#endif

		f1_parameters[8]		= 	dEdx_MeV_g_cm2;

		// d_min_Gy = 0 --> causes problem in SPIFF algorithm. Replaced therefore by minimal dose
		if( *rdd_model != RDD_ExtTarget )
			f1_parameters[3]		=	rdd_parameter[1];

		// d_max_Gy
		if(*rdd_model == RDD_Site){
			float tmp = 0.0f;
//			f1_parameters[4]		=	1.0f / (density_kg_m3 * pi * f1_parameters[1]*f1_parameters[1]) * (LET_J_m - dEdx_J_m);
			f1_parameters[4]		=	SGP_RDD_Katz_site_Gy(&tmp,&alpha,&(f1_parameters[1]),&(f1_parameters[2]),&LET_J_m,&density_kg_m3,&Katz_dEdx_coeff_J_m,&Katz_point_coeff_Gy);
		}else if(*rdd_model == RDD_KatzPoint){
//			f1_parameters[4]		=	C_J_m * Z_eff*Z_eff / (2.0f * pi * f1_parameters[1]*f1_parameters[1] * beta*beta * alpha * density_kg_m3) * pow(1.0f - f1_parameters[1] / f1_parameters[2], 1.0f / alpha);
			f1_parameters[4]		=	SGP_RDD_Katz_point_Gy(&(f1_parameters[1]),&alpha,&(f1_parameters[2]),&Katz_point_coeff_Gy);
		} else {
			f1_parameters[4]		=	0.;
		}

		// single impact fluence
		f1_parameters[6]	= 1.0f / (pi * (f1_parameters[2] * m_to_cm) * (f1_parameters[2] * m_to_cm));	// single_impact_fluence [1/cm2]
		// single_impact_dose
		if(*rdd_model == RDD_Site){
			f1_parameters[7]		=	f1_parameters[0] * MeV_g_to_J_kg * f1_parameters[6];				// LET * fluence
		}else{
			f1_parameters[7]		=	f1_parameters[8] * MeV_g_to_J_kg * f1_parameters[6];				// dEdx * fluence
		}
	}

	if( *rdd_model == RDD_Geiss){
		f1_parameters[1]	=	rdd_parameter[0];															// "r_min_m" = a0
		// Normalization to match with LET
		float	tmp			= (float)(0.5f + log(f1_parameters[2] / f1_parameters[1]));
		tmp					*= 2.0f * pi * (f1_parameters[1] * m_to_cm) * (f1_parameters[1] * m_to_cm);
		f1_parameters[5]	= f1_parameters[0] * MeV_g_to_J_kg / tmp;										// k = LET / tmp
		f1_parameters[6]	= 1.0f / (pi * (f1_parameters[2] * m_to_cm) * (f1_parameters[2] * m_to_cm));	// single_impact_fluence [1/cm2]
		f1_parameters[7]	= f1_parameters[0] * MeV_g_to_J_kg * f1_parameters[6];							// single_impact_dose = LET * single_impact_fluence
		f1_parameters[4]	= f1_parameters[5];																// d_max_Gy = k
		f1_parameters[3]	= f1_parameters[5] * f1_parameters[1]*f1_parameters[1] / (f1_parameters[2]*f1_parameters[2]);				// d_min_Gy
		f1_parameters[8]	= f1_parameters[0];																// dEdx = LET
	}

#ifdef _DEBUG
	fprintf(debf,"%send SGP_RDD_f1_parameters\n",isp);
	indnt_dec();
#endif

}

void SGP_D_RDD_Gy	(	long*	n,
		float*	r_m,
		/* radiation field parameters */
		float*	E_MeV_u,
		long*	particle_no,
		/* detector parameters */
		long*	material_no,
		/* radial dose distribution model */
		long*	rdd_model,
		long*	n_rdd_parameter,
		float*	rdd_parameter,
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		float*	D_RDD_Gy)
{
#ifdef _DEBUG
	indnt_init();
	indnt_inc();
	fprintf(debf,"%sbegin SGP_D_RDD_Gy\n",isp);
#endif

	// conversion through int
#ifdef _R
	int n_int = (int)(*n);
	*n = (long)n_int;

	int n_rdd_parameter_int = (int)(*n_rdd_parameter);
	*n_rdd_parameter = (long)n_rdd_parameter_int;

	int n_er_parameter_int = (int)(*n_er_parameter);
	*n_er_parameter = (long)n_er_parameter_int;

	int rdd_model_int = (int)(*rdd_model);
	*rdd_model = (long)rdd_model_int;

	int er_model_int = (int)(*er_model);
	*er_model = (long)er_model_int;

	int material_no_int	= (int)(*material_no);
	*material_no = (long)material_no_int;

	long ii;
	int particle_no_int;
	for(ii = 0 ; ii < *n ; ii++){
		particle_no_int = (int)(particle_no[ii]);
		particle_no[ii] = (long)particle_no_int;
//		printf("output particle_no[%ld]=%ld\n", i , particle_no_long[i]);
	}
#endif

#ifdef _DEBUG
	fprintf(debf,"%sn = %ld\n", isp, *n);
	fprintf(debf,"%sModel = %ld (no of parameters : %ld) \n", isp, *rdd_model, *n_rdd_parameter);
	fprintf(debf,"%sMaterial = %ld\n", isp, *material_no);
#endif

	long 		i;
	long		n_tmp			= 1;
	// Get f1 parameters
	long 		n_f1_parameters = 9;
	float*		f1_parameters 	= (float*)calloc(n_f1_parameters, sizeof(float));
	SGP_RDD_f1_parameters(			/* radiation field parameters */
			E_MeV_u,
			particle_no,
			/* detector parameters */
			material_no,
			/* radial dose distribution model */
			rdd_model,
			n_rdd_parameter,
			rdd_parameter,
			/* electron range model */
			er_model,
			n_er_parameter,
			er_parameter,
			/* calculated parameters */
			&n_f1_parameters,
			f1_parameters);

	// Get material data
	float 		density_g_cm3, electron_density_m3, I_eV, alpha_g_cm2_MeV, p_MeV, m_g_cm2;
	SGP_getMaterialData(	&n_tmp,
				material_no,
				&density_g_cm3,
				&electron_density_m3,
				&I_eV,
				&alpha_g_cm2_MeV,
				&p_MeV,
				&m_g_cm2);
	float		density_kg_m3	=	density_g_cm3 * 1000.0f;


	if( *rdd_model == RDD_Test){
		// Loop over all r_m given
		for (i = 0; i < *n; i++){
			D_RDD_Gy[i]		=	0.0f;
			if (r_m[i] < f1_parameters[2]){
				D_RDD_Gy[i]		=	f1_parameters[5];
			} else {
				D_RDD_Gy[i]		=	0;
			}
		}
	}// end RDD_Test


	if( *rdd_model == RDD_KatzPoint || *rdd_model == RDD_Site || *rdd_model == RDD_ExtTarget){
		// Get beta, Z and Zeff
		float	beta		= 0.0f;
		long	Z			= 0;
		float	Z_eff		= 0.0f;
		SGP_beta_from_particle_no(	&n_tmp,
									E_MeV_u,
									particle_no,
									&beta);
		SGP_Z_from_particle_no(	&n_tmp,
								particle_no,
								&Z);
		SGP_effective_charge_from_beta(	&n_tmp,
										&beta,
										&Z,
										&Z_eff);

		//////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////
		float alpha				= 	1.667f;
		float w_el_keV			=	2.0f * electron_mass_MeV_c2 * 1000.0f * beta*beta / (1 - beta*beta);
		if(w_el_keV <= 1.0f){
			alpha					= 1.079f;
		}
		//////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////

		// Loop over all r_m given
		float Katz_point_coeff_Gy = SGP_RDD_Katz_point_coeff_Gy(&(f1_parameters[5]),&Z_eff,&beta,&alpha,&density_kg_m3,&(f1_parameters[2]));
		float a0 			  = rdd_parameter[1];
#ifdef _DEBUG
	fprintf(debf,"%sKatzPoint_coeff = %g\n", isp, Katz_point_coeff_Gy);
#endif
		for (i = 0; i < *n; i++){
			D_RDD_Gy[i]				=	0.0f;												// r < r_min_m (for RDD_KatzPoint) or r > r_max_m --> D = 0
			if (r_m[i] >= f1_parameters[1] && r_m[i] <= f1_parameters[2] && *rdd_model != RDD_ExtTarget){					// in between r_min and r_max --> D = KatzPoint
				//D_RDD_Gy[i]				= f1_parameters[5] * Z_eff*Z_eff / (2.0f * pi * r_m[i]*r_m[i] * beta*beta * alpha * density_kg_m3) * pow(1.0f - r_m[i] / f1_parameters[2], 1.0f / alpha);
				D_RDD_Gy[i]				= SGP_RDD_Katz_point_Gy(&(r_m[i]),&alpha,&(f1_parameters[2]),&Katz_point_coeff_Gy);
				D_RDD_Gy[i]				= FMAX(D_RDD_Gy[i], f1_parameters[3]);				// Cut-off low doses
			}
			if (r_m[i] <= f1_parameters[1] && *rdd_model == RDD_Site){						// r < r_min_m (for RDD_Site) --> D = d_max_Gy
				D_RDD_Gy[i]				= f1_parameters[4];
			}
			if (r_m[i] <= f1_parameters[2] && *rdd_model == RDD_ExtTarget){
				D_RDD_Gy[i]				= SGP_RDD_Katz_ext_Gy(&(r_m[i]),&a0,&alpha,&(f1_parameters[1]),&(f1_parameters[2]),&Katz_point_coeff_Gy);
			}
		}
	}// end RDD_KatzPoint & RDD_Site

	if( *rdd_model == RDD_Geiss){
		// Loop over all r_m given
		for (i = 0; i < *n; i++){
			D_RDD_Gy[i]		=	0.0f;
			if (r_m[i] < f1_parameters[1]){
				D_RDD_Gy[i]		=	f1_parameters[5];
			}
			if ((f1_parameters[1] <= r_m[i]) && (r_m[i] <= f1_parameters[2])){
				D_RDD_Gy[i]		=	f1_parameters[5] * (f1_parameters[1] / r_m[i]) * (f1_parameters[1] / r_m[i]);
			}
		}
	}// end RDD_Geiss

	free(f1_parameters);

#ifdef _DEBUG
	fprintf(debf,"%send SGP_D_RDD_Gy\n", isp);
	indnt_dec();
#endif
}

typedef struct {
	long*	n;
	float*	r_m;
	/* radiation field parameters */
	float*	E_MeV_u;
	long*	particle_no;
	/* detector parameters */
	long*	material_no;
	/* radial dose distribution model */
	long*	rdd_model;
	long*	n_rdd_parameter;
	float*	rdd_parameter;
	/* electron range model */
	long*	er_model;
	long*	n_er_parameter;
	float*	er_parameter;
	/* calculated parameters */
	float*	D_RDD_Gy;
	float   D0;
} SGP_D_RDD_Gy_parameters;

float SGP_D_RDD_Gy_solver( float r , void * params ){
	SGP_D_RDD_Gy_parameters* params_struct = (SGP_D_RDD_Gy_parameters*)(params);
	*((*params_struct).n) = 1;
	params_struct->r_m = (float*)calloc(1,sizeof(float));
	(params_struct->r_m)[0] = r;
	SGP_D_RDD_Gy(	params_struct->n,
					params_struct->r_m,
					params_struct->E_MeV_u,
					params_struct->particle_no,
					params_struct->material_no,
					params_struct->rdd_model,
					params_struct->n_rdd_parameter,
					params_struct->rdd_parameter,
					params_struct->er_model,
					params_struct->n_er_parameter,
					params_struct->er_parameter,
					params_struct->D_RDD_Gy);
	free(params_struct->r_m);
	return (params_struct->D_RDD_Gy)[0]-(params_struct->D0);
}


void SGP_r_RDD_m	(	long*	n,
		float*	D_RDD_Gy,
		/* radiation field parameters */
		float*	E_MeV_u,
		long*	particle_no,
		/* detector parameters */
		long*	material_no,
		/* radial dose distribution model */
		long*	rdd_model,       /* */
		long*	n_rdd_parameter, /* number of rdd parameters */
		float*	rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		float*	r_RDD_m)
{
	#ifdef _DEBUG
		indnt_init();
		indnt_inc();
		fprintf(debf,"%sbegin SGP_r_RDD_m\n",isp);
	#endif

		// conversion through int
	#ifdef _R
		int n_int = (int)(*n);
		*n = (long)n_int;

		int n_rdd_parameter_int = (int)(*n_rdd_parameter);
		*n_rdd_parameter = (long)n_rdd_parameter_int;

		int n_er_parameter_int = (int)(*n_er_parameter);
		*n_er_parameter = (long)n_er_parameter_int;

		int rdd_model_int = (int)(*rdd_model);
		*rdd_model = (long)rdd_model_int;

		int er_model_int = (int)(*er_model);
		*er_model = (long)er_model_int;
	#endif

	#ifdef _DEBUG
		fprintf(debf,"%sn = %ld\n", isp, *n);
		fprintf(debf,"%sModel = %ld (no of parameters : %ld) \n", isp, *rdd_model, *n_rdd_parameter);
	#endif

	long 		i;
	long		n_tmp			= 1;
	// Get f1 parameters
	long 		n_f1_parameters = 9;
	float*		f1_parameters 	= (float*)calloc(n_f1_parameters, sizeof(float));
	SGP_RDD_f1_parameters(			/* radiation field parameters */
		E_MeV_u,
		particle_no,
		/* detector parameters */
		material_no,
		/* radial dose distribution model */
		rdd_model,
		n_rdd_parameter,
		rdd_parameter,
		/* electron range model */
		er_model,
		n_er_parameter,
		er_parameter,
		/* calculated parameters */
		&n_f1_parameters,
		f1_parameters);
	// Get material data
	float 		density_g_cm3, electron_density_m3, I_eV, alpha_g_cm2_MeV, p_MeV, m_g_cm2;
	SGP_getMaterialData(	&n_tmp,
				material_no,
				&density_g_cm3,
				&electron_density_m3,
				&I_eV,
				&alpha_g_cm2_MeV,
				&p_MeV,
				&m_g_cm2);
//		float		density_kg_m3	=	density_g_cm3 * 1000.0f;

	if( *rdd_model == RDD_Test){
		// Loop over all doses given
		for (i = 0; i < *n; i++){
			r_RDD_m[i]		=	0.0f;
			if (D_RDD_Gy[i] > 0.0f){
				r_RDD_m[i]		=	f1_parameters[2];
			}
		}
	}// end RDD_Test

	if( *rdd_model == RDD_Geiss){
		// Loop over all doses given
		for (i = 0; i < *n; i++){
			// if D is outside the definition, return -1
			r_RDD_m[i]		=	-1.0f;
			if ((f1_parameters[3] <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= f1_parameters[4])){
				r_RDD_m[i]		=	f1_parameters[1] * (float)sqrt(f1_parameters[5] / D_RDD_Gy[i]);}
		}
	}// end RDD_Geiss

	if( (*rdd_model == RDD_KatzPoint) || (*rdd_model == RDD_Site)){
		SGP_D_RDD_Gy_parameters* params = (SGP_D_RDD_Gy_parameters*)calloc(1,sizeof(SGP_D_RDD_Gy_parameters));
		params->n 				= (long*)calloc(1,sizeof(long));
		params->E_MeV_u 		= E_MeV_u;
		params->particle_no		= particle_no;
		params->material_no		= material_no;
		params->rdd_model 		= rdd_model;
		params->n_rdd_parameter = n_rdd_parameter;
		params->rdd_parameter 	= rdd_parameter;
		params->er_model 		= er_model;
		params->n_er_parameter 	= n_er_parameter;
		params->er_parameter 	= er_parameter;
		params->D_RDD_Gy 		= (float*)calloc(1,sizeof(float));
		// Loop over all doses given

		float critical_d_Gy		= 0.0f;
		float critical_r_m		= rdd_parameter[0] * (1.0f + 1e-6f);
		if(*rdd_model == RDD_Site){
			SGP_D_RDD_Gy	(	&n_tmp,										// Use D(r) to find dose at jump of D_Site
								&critical_r_m,
								/* radiation field parameters */
								E_MeV_u,
								particle_no,
								/* detector parameters */
								material_no,
								/* radial dose distribution model */
								rdd_model,
								n_rdd_parameter,
								rdd_parameter,
								/* electron range model */
								er_model,
								n_er_parameter,
								er_parameter,
								&critical_d_Gy);
		}
		printf("Critical radius = %e m, critical dose = %e Gy\n", critical_r_m, critical_d_Gy);


		for (i = 0; i < *n; i++){
			params->D0 = D_RDD_Gy[i];
			// if D is outside the definition, return -1
			r_RDD_m[i]				=	-1.0f;
			float	solver_accuracy	=	1e-13f;
			if ((f1_parameters[3] <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= f1_parameters[4])){
				if(*rdd_model == RDD_Site && D_RDD_Gy[i] >= critical_d_Gy){
					r_RDD_m[i] = rdd_parameter[0];
				}else{
					r_RDD_m[i] = zriddr(SGP_D_RDD_Gy_solver,
										(void*)params,
										f1_parameters[1],
										f1_parameters[2],
										solver_accuracy);
				}
			}
		}
		free(params->n);
		free(params->D_RDD_Gy);
		free(params);
	}// end RDD_KatzPoint & RDD_Site

	free(f1_parameters);

#ifdef _DEBUG
		fprintf(debf,"%send SGP_r_RDD_m\n", isp);
		indnt_dec();
#endif
}

#endif // SGP_RDD_H_
