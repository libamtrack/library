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

inline float SGP_RDD_Katz_point_coeff_Gy(float* C_J_m,float* Z_eff, float* beta, float* alpha, float* density_kg_m3, float* r_max_m){
	return (*C_J_m) * (*Z_eff)*(*Z_eff) / (2.0f * M_PI * (*beta)*(*beta) * (*alpha) * (*density_kg_m3) * (*r_max_m));
}

inline float SGP_RDD_Katz_point_Gy(float* r_m, float* alpha, float* r_max_m, float* Katz_point_coeff_Gy){
	float x = (*r_m)/(*r_max_m);
	return (*Katz_point_coeff_Gy) * SGP_RDD_Katz_point_kernel(&x,alpha);
}

inline float SGP_RDD_Katz_dEdx_kernel(float* x, float* alpha){
	return (1.0f/(*x) )*pow(1.0f - (*x), 1.0f / (*alpha));
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
	gsl_integration_qags (&F, (*r_min_m)/(*r_max_m), 1.0, 0, 1e-5, 1000, w1, &dEdx_integral, &error);
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


//double SGP_RDD_f1_dEdx_integrand(double r, void * params){
////	float LET = ((float*)params)[0];
//	float density_g_cm3 = ((float*)params)[1];
////	float r_min_m = ((float*)params)[2];
//	float r_max_m = ((float*)params)[3];
//	float C_J_m = ((float*)params)[4];
//	float Z_eff = ((float*)params)[5];
//	float beta = ((float*)params)[6];
//	float alpha = ((float*)params)[7];
//	float	res		=	0.0f;
//	float 	density_kg_m3 	=   1e3 * density_g_cm3;
//	res = r * SGP_RDD_Katz_coeff(C_J_m,Z_eff, beta, alpha, density_kg_m3) * SGP_RDD_Katz_kernel(r_m, r_max_m, alpha);
//	//res	+=	r * C_J_m * Z_eff*Z_eff / (2.0f * M_PI * r * r * beta * beta * alpha * density_kg_m3) * pow(1.0f - r / r_max_m, 1.0f / alpha);
//	res	*=	2.0f * M_PI * density_kg_m3;
//	res	/= 100.0f / density_g_cm3 / MeV_to_J;
//	return	res;
//}
//
//
//
float SGP_RDD_f1_getdEdx(float LET, float density_g_cm3, float r_min_m, float r_max_m, float C_J_m, float Z_eff, float beta, float alpha)
{
	double result, error;

//	gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
//
//	gsl_function F;
//	F.function = &SGP_RDD_f1_dEdx_integrand;
//	float params[8] = {LET, density_g_cm3, r_min_m, r_max_m, C_J_m, Z_eff, beta, alpha};
//	F.params = params;
//
//	gsl_integration_qags (&F, r_min_m, r_max_m, 0, 1e-5, 1000, w1, &result, &error);
//
//	gsl_integration_workspace_free (w1);

	return	result;
}

//double
//geometryFunctionPhi (double r0, double a0, double r)
//{
//	double res = 0.;
//	double factor = 0.;
//	gsl_complex carg, cres;
//	if (r <= fabs (r0 - a0))
//	{
//		if (r0 >= a0)
//			res = 0.;
//		else
//			res = M_PI;
//	}
//	else
//	{
//		factor = gsl_pow_2 (a0) - gsl_pow_2 (r0 - r);
//		factor /= gsl_pow_2 (r + r0) - gsl_pow_2 (a0);
//		GSL_SET_COMPLEX (&carg, sqrt (factor), 0.);
//		cres = gsl_complex_arctan (carg);
//		res = 2. * GSL_REAL (cres);
//	}
//	return res;
//}
//
//double SGP_RDD_f1_extTarget_integrand(double r, void * params){
//
//	float LET = ((float*)params)[0];
//	float density_g_cm3 = ((float*)params)[1];
//	float r_min_m = ((float*)params)[2];
//	float r_max_m = ((float*)params)[3];
//	float C_J_m = ((float*)params)[4];
//	float Z_eff = ((float*)params)[5];
//	float beta = ((float*)params)[6];
//	float alpha = ((float*)params)[7];
//
//	float dEdx_MeV_g_cm2 = SGP_RDD_f1_getdEdx(LET,density_g_cm3,r_min_m,r_max_m,C_J_m,Z_eff,beta,alpha);
//
//	float doseSite = 0.0f;
//
////	f1_parameters[3]		=	C_J_m * Z_eff*Z_eff / (2.0f * pi * f1_parameters[2]*f1_parameters[2] * beta*beta * alpha * density_kg_m3) * pow(1.0f - f1_parameters[2] / f1_parameters[2], 1.0f / alpha);
////
////	// single impact fluence
////	f1_parameters[6]	= 1.0f / (pi * (f1_parameters[2] * m_to_cm) * (f1_parameters[2] * m_to_cm));	// single_impact_fluence [1/cm2]
////	// single_impact_dose
//
////	f1_parameters[7]		=	f1_parameters[0] * MeV_g_to_J_kg * f1_parameters[6];				// LET * fluence
//
//	float	dEdx_J_m		=	dEdx_MeV_g_cm2 * 100.0f * density_g_cm3 * MeV_to_J;
//	float 	density_kg_m3	= 	density_g_cm3 * 1e3;
//	float 	LET_J_m			=	LET * 1.602e-13 * density_g_cm3 * 100.0f;
//
//	if (r >= r_min_m && r <= r_max_m){					// in between r_min and r_max --> D = KatzPoint
//		doseSite = C_J_m * Z_eff*Z_eff / (2.0f * pi * r*r * beta*beta * alpha * density_kg_m3) * pow(1.0f - r / r_max_m, 1.0f / alpha);
//	}
//	if (r <= r_min_m){						// r < r_min_m (for RDD_Site) --> D = d_max_Gy
//		doseSite = 1.0f / (density_kg_m3 * pi * r_min_m*r_min_m) * (LET_J_m - dEdx_J_m);
//	}
//
//	float	res		=	0.0f;
//	float 	a0		= 	1e-8;
//
//	res = (1 / (pi* a0*a0)) * doseSite * geometryFunctionPhi(r_min_m,a0,r);
//
//	return	res;
//}
//
//
//float SGP_RDD_f1_getExtTargetDose(float LET, float density_g_cm3, float r_min_m, float r_max_m, float C_J_m, float Z_eff, float beta, float alpha)
//{
//	double result, error;
//
//	gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
//
//	gsl_function F;
//	F.function = &SGP_RDD_f1_extTarget_integrand;
//	float params[8] = {LET, density_g_cm3, r_min_m, r_max_m, C_J_m, Z_eff, beta, alpha};
//	F.params = params;
//
//	gsl_integration_qags (&F, r_min_m, r_max_m, 0, 1e-5, 1000, w1, &result, &error);
//
//	gsl_integration_workspace_free (w1);
//
//	return	result;
//}

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
		f1_parameters[1] 	= 0.0f;													// r_min_m
		f1_parameters[6]	= 1.0f / (pi * f1_parameters[2]*f1_parameters[2]);		// pi * r_max_m^2 = Track area -> single_impact_fluence [1/m2]
		f1_parameters[5]	= f1_parameters[0] * f1_parameters[6];					// LET / track area = Norm.constant k
		f1_parameters[3]	= f1_parameters[5];										// d_min_Gy = k
		f1_parameters[4]	= f1_parameters[5];										// d_max_Gy = k
		f1_parameters[7]	= f1_parameters[0] * f1_parameters[5] * MeV_g_to_J_kg;	// single_impact_dose = LET * fluence;
		f1_parameters[6]	/= 10000;												// single_impact_fluence [1/cm2]
		f1_parameters[8]	= f1_parameters[0];										// dEdx = LET
	}

	if( *rdd_model == RDD_KatzPoint || *rdd_model == RDD_Site){
		//////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////
		float alpha				= 	1.667f;
		float w_el_keV			=	2.0f * electron_mass_MeV_c2 * 1000.0f * beta*beta / (1 - beta*beta);
		if(w_el_keV <= 1.0f){
			alpha					= 1.079f;
		}
		//////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////
		f1_parameters[1]		=	rdd_parameter[0];								// r_min_m

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
		float	N_el_cm3		=	electron_density_m3 / (100*100*100);
		float	C_J_cm			=	2.0f * pi * N_el_cm3 * (e_esu*e_esu*e_esu*e_esu) / (electron_mass_g * c_cm_s *c_cm_s) * 1e-7;   // energy constant [J/cm] not [erg/cm] hence 10^-7
		float	C_J_m			=	C_J_cm * 100.0f;
		f1_parameters[5]		= 	C_J_m;											// Norm.constant k


		// Get dEdx by simple integration from r_min_m (in case of RDD_Site = a0) to r_max_m
		float	dEdx_J_m		=	0.0f;
		float	LET_J_m 		=	f1_parameters[0] * 1.602e-13 * density_g_cm3 * 100.0f;
		long 	j;
		long 	n_steps = 100;
		float 	log_r_min_m		=	log10(f1_parameters[1]);
		float 	log_r_max_m		=	log10(f1_parameters[2]);
		float	log_step_size	=	(log_r_max_m - log_r_min_m) / n_steps;
		for (j = 0; j < n_steps; j++){
			float	r_int_m			=	pow(10, log_r_min_m + (j + 0.5f) * log_step_size);
			float	dr_int_m		=	pow(10, log_r_min_m + (j + 1.0f) * log_step_size) - pow(10, log_r_min_m + (j + 0.0f) * log_step_size);
			dEdx_J_m				+=	r_int_m * dr_int_m *		\
										C_J_m * Z_eff*Z_eff / (2.0f * pi * r_int_m*r_int_m * beta*beta * alpha * density_kg_m3) * pow(1.0f - r_int_m / f1_parameters[2], 1.0f / alpha);
		}
		dEdx_J_m					*=	2.0f * pi * density_kg_m3;
		float	dEdx_MeV_g_cm2		=	dEdx_J_m / 100.0f / density_g_cm3 / MeV_to_J;

		f1_parameters[8]		= 	dEdx_MeV_g_cm2;

		// d_min_Gy = 0 --> causes problem in SPIFF algorithm. Replaced therefore by minimal dose
		f1_parameters[3]		=	1e-30f;
//		f1_parameters[3]		=	C_J_m * Z_eff*Z_eff / (2.0f * pi * f1_parameters[2]*f1_parameters[2] * beta*beta * alpha * density_kg_m3) * pow(1.0f - f1_parameters[2] / f1_parameters[2], 1.0f / alpha);
		// d_max_Gy
		if(*rdd_model == RDD_Site){
			f1_parameters[4]		=	1.0f / (density_kg_m3 * pi * f1_parameters[1]*f1_parameters[1]) * (LET_J_m - dEdx_J_m);
		}else{
			f1_parameters[4]		=	C_J_m * Z_eff*Z_eff / (2.0f * pi * f1_parameters[1]*f1_parameters[1] * beta*beta * alpha * density_kg_m3) * pow(1.0f - f1_parameters[1] / f1_parameters[2], 1.0f / alpha);
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


	if( *rdd_model == RDD_KatzPoint || *rdd_model == RDD_Site){
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
		for (i = 0; i < *n; i++){
			D_RDD_Gy[i]				=	0.0f;												// r < r_min_m (for RDD_KatzPoint) or r > r_max_m --> D = 0
			if (r_m[i] >= f1_parameters[1] && r_m[i] <= f1_parameters[2]){					// in between r_min and r_max --> D = KatzPoint
				D_RDD_Gy[i]				= f1_parameters[5] * Z_eff*Z_eff / (2.0f * pi * r_m[i]*r_m[i] * beta*beta * alpha * density_kg_m3) * pow(1.0f - r_m[i] / f1_parameters[2], 1.0f / alpha);
			}
			if (r_m[i] <= f1_parameters[1] && *rdd_model == RDD_Site){						// r < r_min_m (for RDD_Site) --> D = d_max_Gy
				D_RDD_Gy[i]				= f1_parameters[4];
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
		for (i = 0; i < *n; i++){
			params->D0 = D_RDD_Gy[i];
			// if D is outside the definition, return -1
			r_RDD_m[i]				=	-1.0f;
			float	solver_accuracy	=	1e-13f;
			if ((f1_parameters[3] <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= f1_parameters[4])){
					r_RDD_m[i] = zriddr(	SGP_D_RDD_Gy_solver,
											(void*)params,
											f1_parameters[1],
											f1_parameters[2],
											solver_accuracy);
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
