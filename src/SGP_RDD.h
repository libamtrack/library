#ifndef SGP_RDD_H_
#define SGP_RDD_H_

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
