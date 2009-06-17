#ifndef SGP_RDD_H_
#define SGP_RDD_H_


//void SGP_RDD_f1_parameters(	long*			n,
//							float*			E_MeV_u,
//							long*			particle_no,
//							char*			material_name,
//							float*			parameter,
//							float*			LET_MeV_cm2_g,
//							float*			r_min_m,
//							float*			r_max_m,
//							float*			d_min_Gy,
//							float*			d_max_Gy,
//							float*			k_Gy,
//							float*			single_impact_fluence_cm2,
//							float*			single_impact_dose_Gy);

void SGP_RDD_f1_parameters(	long*			n,
		/* radial dose distribution model */
		long*	rdd_model,       /* model: 0 - test, 1 - Katz, 2 - LEM*/
		long*	n_rdd_parameter, /* number of rdd parameters */
		float*	rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		/* calculated parameters */
		long * n_f1_parameters,
		float * f1_parameters);

void SGP_RDD_f1_parameters_Geiss(	long*			n,
									float*			E_MeV_u,
									long*			particle_no,
									char*			material_name,
									float*			parameter,
									float*			LET_MeV_cm2_g,
									float*			r_min_m,
									float*			r_max_m,
									float*			d_min_Gy,
									float*			d_max_Gy,
									float*			k_Gy,
									float*			single_impact_fluence_cm2,
									float*			single_impact_dose_Gy);

void SGP_RDD_f1_parametersS(	long*			n,
		float*			E_MeV_u,
		long*			particle_no,
		char**			material_name,
		float*			parameter,
		float*			LET_MeV_cm2_g,
		float*			r_min_m,
		float*			r_max_m,
		float*			d_min_Gy,
		float*			d_max_Gy,
		float*			k_Gy,
		float*			single_impact_fluence_cm2,
		float*			single_impact_dose_Gy);

void SGP_D_RDD_Gy	(	long*	n,
		float*	r_m,
		/* radial dose distribution model */
		long*	rdd_model,       /* */
		long*	n_rdd_parameter, /* number of rdd parameters */
		float*	rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		float*	D_RDD_Gy);

/*
void SGP_D_RDD_Gy	(	long*	n,
						float*	r_m,
						float*	E_MeV_u,
						long*	particle_no,
						char*	material_name,
						float*	parameter,
						float*	D_RDD_Gy);*/

/* obsolete: to be deleted*/
void SGP_D_RDD_GyS(		long*	n,
		float*	r_m,
		float*	E_MeV_u,
		int*	particle_no,
		char**	material_name,
		float*	parameter,
		float*	D_RDD_Gy){};

void SGP_r_RDD_m	(	long*	n,
		float*	D_RDD_Gy,
		/* radial dose distribution model */
		long*	rdd_model,       /* */
		long*	n_rdd_parameter, /* number of rdd parameters */
		float*	rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		float*	r_m);

void SGP_r_RDD_m_old	(	long*	n,
				float*	D_Gy,
				float*	E_MeV_u,
				long*	particle_no,
				char*	material_name,
				float*	parameter,
				float*	r_RDD_m);

//void SGP_RDD_f1_parameters(	long*			n,
//							float*			E_MeV_u,
//							long*			particle_no,
//							char*			material_name,
//							float*			parameter,					// general formalism for later implementation of other RDDs
//							float*			LET_MeV_cm2_g,				// here: n_parameter = 1, parameter[0] = r_min_m
//							float*			r_min_m,
//							float*			r_max_m,
//							float*			d_min_Gy,
//							float*			d_max_Gy,
//							float*			k_Gy,
//							float*			single_impact_fluence_cm2,
//							float*			single_impact_dose_Gy)

void SGP_RDD_f1_parameters(	long*			n,
		/* radial dose distribution model */
		long*	rdd_model,       /* */
		long*	n_rdd_parameter, /* number of rdd parameters */
		float*	rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		/* calculated parameters */
		long * n_f1_parameters,
		float * f1_parameters)
{
#ifdef _DEBUG
	indnt_inc();
	fprintf(debf,"%sbegin SGP_RDD_f1_parameters\n",isp);
#endif

	if( *rdd_model == RDD_Test){
		f1_parameters[0] = 137.;
	}

	if( *rdd_model == RDD_KatzPoint){
		f1_parameters[0] = 137.;
	}

	if( *rdd_model == RDD_Geiss){

		float*	LET_MeV_cm2_g				=	(float*)calloc(*n, sizeof(float));
		float*	r_min_m						=	(float*)calloc(*n, sizeof(float));
		float*	r_max_m						=	(float*)calloc(*n, sizeof(float));
		float*	d_min_Gy					=	(float*)calloc(*n, sizeof(float));
		float*	d_max_Gy					=	(float*)calloc(*n, sizeof(float));
		float*	k_Gy						=	(float*)calloc(*n, sizeof(float));
		float*	single_impact_fluence_cm2	=	(float*)calloc(*n, sizeof(float));
		float*	single_impact_dose_Gy		=	(float*)calloc(*n, sizeof(float));

		float*	E_MeV_u						=	(float*)calloc(*n, sizeof(float));
		long*	particle_no					=	(long*)calloc(*n, sizeof(long));
		char	material_name[100];
		float*	parameter					=	(float*)calloc(*n, sizeof(float));

		long i;
		for( i = 0 ; i < *n ; i++){
			E_MeV_u[i] = rdd_parameter[0];
			particle_no[i] = (long)(rdd_parameter[1]);
			parameter[i] = rdd_parameter[3];
		}

		switch( (int)(rdd_parameter[2]) ){
		case 0:
			strcpy(material_name,"Water, Liquid");
			break;
		case 1:
			strcpy(material_name,"Aluminum Oxide");
			break;
		case 2:
			strcpy(material_name,"Aluminum");
			break;
		case 3:
			strcpy(material_name,"PMMA");
			break;
		default:
			strcpy(material_name,"");
			break;
		}

		// Get maximum electron range
		SGP_max_electron_range_m(	n,
										E_MeV_u,
										particle_no,
										material_name,
										r_max_m);

		f1_parameters[2] = r_max_m[0];

		// Get the LET
		// Usually the routine looks up the LET value corresponding to
		// the E values given. This can, however, be incorrect in case of
		// dose-weighted mean values. There is therefore the possibility
		// to pass on LET values different from the PSTAR tables.
		// In that case, n NEGATIVE E values are given, followed by
		// n LET value, coerced in the E_MeV_u array

		if(E_MeV_u[0] > 0){									// Usual case: get LET from PSTAR table
			long nn = 1;
			SGP_LET_MeV_cm2_g(	&nn,
							E_MeV_u,
							particle_no,
							material_name,
							LET_MeV_cm2_g);
		}else{												// Alternatively: use provided LET values
			long	i;
			for (i = 0; i < *n; i++){
				LET_MeV_cm2_g[i]	=	E_MeV_u[i + (*n)];
				E_MeV_u[i]		=	-1.0f * E_MeV_u[i];
			}
		}
		f1_parameters[0] = LET_MeV_cm2_g[0];

		// Loop over all particles and energies given
		for (i = 0; i < *n; i++){
		//		fprintf(debf,"%sr_max_m[%ld] = %e\n", i, r_max_m[i]);
			r_min_m[i]	= parameter[0];
			//fprintf(debf,"%sr_min_m[%ld] = %e\n", i, r_min_m[i]);
			if(r_min_m[i] >= r_max_m[i]){
				r_min_m[i] = r_max_m[i];
			}

			// Normalization to match with LET
			float	tmp						=	(float)(0.5f + log(r_max_m[i] / r_min_m[i]));

				tmp								*=	2.0f * pi * (r_min_m[i] * m_to_cm) * (r_min_m[i] * m_to_cm);
				//fprintf(debf,"%stmp = %e\n", tmp);

				k_Gy[i]							=	LET_MeV_cm2_g[0] * MeV_g_to_J_kg / tmp;
#ifdef _DEBUG
				fprintf(debf,"%sLET[%ld] = %e\n", isp, i, LET_MeV_cm2_g[0]);
#endif

				tmp								=	1 / (pi * (r_max_m[i] * m_to_cm) * (r_max_m[i] * m_to_cm));
				single_impact_fluence_cm2[i]	=	tmp;
				single_impact_dose_Gy[i]		=	LET_MeV_cm2_g[0] * MeV_g_to_J_kg * tmp;

				d_max_Gy[i]						=	k_Gy[i];
				d_min_Gy[i]						=	k_Gy[i] * r_min_m[i] * r_min_m[i] / (r_max_m[i] * r_max_m[i]);

			}

		f1_parameters[1] = r_min_m[0];
		f1_parameters[3] = d_min_Gy[0];
		f1_parameters[4] = d_max_Gy[0];
		f1_parameters[5] = k_Gy[0];
		f1_parameters[6] = single_impact_fluence_cm2[0];
		f1_parameters[7] = single_impact_dose_Gy[0];

		free(LET_MeV_cm2_g);
		free(r_min_m);
		free(r_max_m);
		free(d_min_Gy);
		free(d_max_Gy);
		free(k_Gy);
		free(single_impact_fluence_cm2);
		free(single_impact_dose_Gy);

		free(E_MeV_u);
		free(particle_no);
		free(parameter);

	}


#ifdef _DEBUG
	fprintf(debf,"%send SGP_RDD_f1_parameters\n",isp);
	indnt_dec();
#endif

}

void SGP_RDD_f1_parameters_Geiss(	long*			n,
									float*			E_MeV_u,
									long*			particle_no,
									char*			material_name,
									float*			parameter,
									float*			LET_MeV_cm2_g,
									float*			r_min_m,
									float*			r_max_m,
									float*			d_min_Gy,
									float*			d_max_Gy,
									float*			k_Gy,
									float*			single_impact_fluence_cm2,
									float*			single_impact_dose_Gy){
		/* radial dose distribution model */
		long	rdd_model = 2;       /* model: 0 - test, 1 - Katz, 2 - LEM*/
		long	n_rdd_parameter = 4; /* number of rdd parameters */
		float*	rdd_parameter;   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long	er_model = 0;
		long	n_er_parameter = 1;
		float*	er_parameter;
		/* calculated parameters */
		long  n_f1_parameters = 8;
		float * f1_parameters;

		rdd_parameter = (float*)calloc(n_rdd_parameter, sizeof(float));
		er_parameter = (float*)calloc(n_er_parameter, sizeof(float));
		f1_parameters = (float*)calloc(n_f1_parameters, sizeof(float));

		if( strcmp(material_name,"Water, Liquid") == 0)
			rdd_parameter[2] = 0.;
		if( strcmp(material_name,"Aluminum Oxide") == 0)
			rdd_parameter[2] = 1.;
		if( strcmp(material_name,"Aluminum") == 0)
			rdd_parameter[2] = 2.;
		if( strcmp(material_name,"PMMA") == 0)
			rdd_parameter[2] = 3.;
		rdd_parameter[0] = E_MeV_u[0];
		rdd_parameter[1] = (float)(particle_no[0]);
		rdd_parameter[3] = parameter[0];

		SGP_RDD_f1_parameters(n,&rdd_model,&n_rdd_parameter,rdd_parameter,&er_model, &n_er_parameter, er_parameter, &n_f1_parameters, f1_parameters);

		long i;
		for( i = 0 ; i < *n ; i++){
			LET_MeV_cm2_g[i] = f1_parameters[0];
			r_min_m[i] = f1_parameters[1];
			r_max_m[i] = f1_parameters[2];
			d_min_Gy[i] = f1_parameters[3];
			d_max_Gy[i] = f1_parameters[4];
			k_Gy[i] = f1_parameters[5];
			single_impact_fluence_cm2[i] = f1_parameters[6];
			single_impact_dose_Gy[i] = f1_parameters[7];
		}

		free(rdd_parameter);
		free(er_parameter);
		free(f1_parameters);
}




void SGP_RDD_f1_parametersS(	long*			n,
		float*			E_MeV_u,
		long*			particle_no,
		char**			material_name,
		float*			parameter,
		float*			LET_MeV_cm2_g,
		float*			r_min_m,
		float*			r_max_m,
		float*			d_min_Gy,
		float*			d_max_Gy,
		float*			k_Gy,
		float*			single_impact_fluence_cm2,
		float*			single_impact_dose_Gy)
{
	//	SGP_RDD_f1_parameters(	n,
	//							E_MeV_u,
	//							particle_no,
	//							*material_name,
	//							parameter,
	//							LET_MeV_cm2_g,
	//							r_min_m,
	//							r_max_m,
	//							d_min_Gy,
	//							d_max_Gy,
	//							k_Gy,
	//							single_impact_fluence_cm2,
	//							single_impact_dose_Gy);
}



//void SGP_D_RDD_Gy	(	long*	n,
//						float*	r_m,
//						/* radial dose distribution model */
//						long*	rdd_model,
//						long*	n_rdd_parameter,
//						float*	rdd_parameter,
//						/* electron range model */
//						long*	er_model,
//						long*	n_er_parameter,
//						float*	er_parameter,
//						float*	D_RDD_Gy);

/*
float RDD_TEST_Gy(	float	r_m,
					long*	rdd_model,
					long*	n_rdd_parameter,
					float*	rdd_parameter,
					long*	er
*/
void SGP_D_RDD_Gy	(	long*	n,
		float*	r_m,
		/* radial dose distribution model */
		long*	rdd_model,
		long*	n_rdd_parameter,
		float*	rdd_parameter,
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		float*	D_RDD_Gy)

//
//void SGP_D_RDD_Gy	(	long*	n,
//						float*	r_m,
//						float*	E_MeV_u,
//						long*	particle_no,
//						char*	material_name,
//						float*	parameter,
//						float*	D_RDD_Gy)
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

	long i;

	if( *rdd_model == RDD_KatzPoint){

#ifdef _DEBUG
		fprintf(debf,"%sRDD_KatzPoint, rdd_parameter[0] = %g, rdd_parameter[1] = %g, rdd_parameter[2] = %g\n",
				isp, rdd_parameter[0], rdd_parameter[1], rdd_parameter[2]);
#endif

		long	n_tmp		= 1;
		float	beta		= 0.0f;
		long	Z			= 0;
		if(rdd_parameter[0]	< 0){
			beta		= -1.0f * rdd_parameter[0];
			Z			= (long)rdd_parameter[1];
		}else{
			float	E_MeV_u		= rdd_parameter[0];
			long	particle_no	= (long)rdd_parameter[1];
			beta		= 0.0f;
			SGP_beta_from_particle_index(	&n_tmp,
										&E_MeV_u,
										&particle_no,
										&beta);
			// PRELIMINARY!!!
			Z			= 1;
		}
		long	material_no	= (long)rdd_parameter[2];
		float	Z_eff		= 0.0f;
		SGP_effective_charge_from_beta(	&n_tmp,
										&beta,
										&Z,
										&Z_eff);

		char	material_name[256];
		char*	pmn[1];
		pmn[0]	= (char*)&material_name;
		getMaterialName(&material_no, material_name);

		float 	density_g_cm3;
		float	electron_density_m3;
		float	I_eV;
		float	alpha_g_cm2_MeV;
		float	p_MeV;
		float	m_g_cm2;

		SGP_getMaterialData(	&n_tmp,
						pmn,
						&density_g_cm3,
						&electron_density_m3,
						&I_eV,
						&alpha_g_cm2_MeV,
						&p_MeV,
						&m_g_cm2);

		float	N_el_cm3		=	electron_density_m3 / (100*100*100);
		float	p_g_cm3			=	density_g_cm3;
		float	I_erg			=	I_eV / 6.24150974e11;

		float	C				=	2.0f * pi * N_el_cm3 * (e_esu*e_esu*e_esu*e_esu) / (electron_mass_g * c_cm_s *c_cm_s) * 1e-7;   // energy constant [J/cm] not [erg/cm] hence 10^-7

		float	w_max_keV		=	2.0f * electron_mass_MeV_c2 * 1000 * beta*beta / (1.0f - beta*beta);

		float	alpha			=	1.667;
		float	k_el			=	6.13e-6f;			// range constant [g*cm^(-2)*keV^(-alpha)]
		if(w_max_keV <= 1){								// low E alpha (< 1 keV)
			alpha			=	1.079;
		}

		float	r_max_cm		=	10e-6f * w_max_keV;

		long i;
		for (i = 0; i < *n; i++){
			float	r_cm			=	r_m[i] * 100.0f;
			if (r_cm <= r_max_cm){
				D_RDD_Gy[i]				=	C * Z_eff*Z_eff / (2.0f * pi * r_cm * beta*beta) * pow(1.0f / r_cm - 1.0f / r_max_cm, 1.0f / alpha);
			}else{
				D_RDD_Gy[i]				=	0.0f;
			}
		}
	}

	if( *rdd_model == RDD_Test){

#ifdef _DEBUG
		fprintf(debf,"%sRDD_Test, rdd_parameter[0] = %g\n", isp, rdd_parameter[0]);
#endif

		long n_f1_parameters = 1;
		float*	f1_parameters = (float*)calloc(n_f1_parameters, sizeof(float));

		SGP_RDD_f1_parameters(	n,
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

		/* rewriting output of SGP_RDD_f1_parameters function
		 * to more convenient format
		 */
		float*	factor = (float*)calloc(*n, sizeof(float));
		for( i = 0 ; i < *n ; i++){
			factor[i] = f1_parameters[0];
		}
		free(f1_parameters);

		/*
		 * calculation of dose
		 */
		for (i = 0; i < *n; i++){
			D_RDD_Gy[i]		=	0.0f;
			if (r_m[i] < rdd_parameter[0]){
				D_RDD_Gy[i]		=	factor[i];
			} else {
				D_RDD_Gy[i]		=	0;
			}
		}

	}


	if( *rdd_model == RDD_Geiss){

#ifdef _DEBUG
		fprintf(debf,"%sRDD_Geiss, rdd_parameter[0] = %g, rdd_parameter[1] = %g, rdd_parameter[2], rdd_parameter[3] = %g\n",
				isp, rdd_parameter[0], rdd_parameter[1], rdd_parameter[2], rdd_parameter[3]);
#endif
		long n_f1_parameters = 8;
		float*	f1_parameters = (float*)calloc(n_f1_parameters, sizeof(float));

		SGP_RDD_f1_parameters(	n,
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

		/* rewriting output of SGP_RDD_f1_parameters function
		 * to more convenient format
		 */
		float*	LET_MeV_cm2_g				=	(float*)calloc(*n, sizeof(float));
		float*	r_min_m						=	(float*)calloc(*n, sizeof(float));
		float*	r_max_m						=	(float*)calloc(*n, sizeof(float));
		float*	d_min_Gy					=	(float*)calloc(*n, sizeof(float));
		float*	d_max_Gy					=	(float*)calloc(*n, sizeof(float));
		float*	k_Gy						=	(float*)calloc(*n, sizeof(float));
		float*	single_impact_fluence_cm2	=	(float*)calloc(*n, sizeof(float));
		float*	single_impact_dose_Gy		=	(float*)calloc(*n, sizeof(float));
		for( i = 0 ; i < *n ; i++){
			LET_MeV_cm2_g[i] = f1_parameters[0];
			r_min_m[i] = f1_parameters[1];
			r_max_m[i] = f1_parameters[2];
			d_min_Gy[i] = f1_parameters[3];
			d_max_Gy[i] = f1_parameters[4];
			k_Gy[i] = f1_parameters[5];
			single_impact_fluence_cm2[i] = f1_parameters[6];
			single_impact_dose_Gy[i] = f1_parameters[7];
		}
		free(f1_parameters);

		//long	i;
		for (i = 0; i < *n; i++){
#ifdef _DEBUG
			fprintf(debf,"%sr_m[%ld] = %g\n", isp, i , r_m[i]);
			fprintf(debf,"%sk_Gy[%ld] = %g\n", isp, i , k_Gy[i]);
#endif

			// Differentiate between the three cases (i) r < r_min (ii) r_min <= r <= r_max (iii) r > r_max
			D_RDD_Gy[i]		=	0.0f;

			if (r_m[i] < r_min_m[i]){
				D_RDD_Gy[i]		=	k_Gy[i];}
			if ((r_min_m[i] <= r_m[i]) && (r_m[i] <= r_max_m[i])){
				D_RDD_Gy[i]		=	k_Gy[i] * (r_min_m[i] / r_m[i]) * (r_min_m[i] / r_m[i]);}
		}

		free(LET_MeV_cm2_g);
		free(r_min_m);
		free(r_max_m);
		free(d_min_Gy);
		free(d_max_Gy);
		free(k_Gy);
		free(single_impact_fluence_cm2);
		free(single_impact_dose_Gy);


	}

	if( *rdd_model == RDD_Site){

#ifdef _DEBUG
		fprintf(debf,"%sRDD_Site, rdd_parameter[0] = %g\n", isp, rdd_parameter[0]);
#endif

		long n_f1_parameters = 1;
		float*	f1_parameters = (float*)calloc(n_f1_parameters, sizeof(float));

		SGP_RDD_f1_parameters(	n,
				rdd_model,       /* model: 0 - test, 1 - Katz, 2 - LEM*/
				n_rdd_parameter, /* number of rdd parameters */
				rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
				/* electron range model */
				er_model,
				n_er_parameter,
				er_parameter,
				/* calculated parameters */
				&n_f1_parameters,
				f1_parameters);

		/* rewriting output of SGP_RDD_f1_parameters function
		 * to more convenient format
		 */
		float*	factor = (float*)calloc(*n, sizeof(float));
		for( i = 0 ; i < *n ; i++){
			factor[i] = f1_parameters[0];
		}
		free(f1_parameters);

		/*
		 * calculation of dose
		 */
		for (i = 0; i < *n; i++){
			D_RDD_Gy[i]		=	0.0f;
			if (r_m[i] < rdd_parameter[0]){
				D_RDD_Gy[i]		=	factor[i];
			} else {
				D_RDD_Gy[i]		=	0;
			}
		}

	}

#ifdef _DEBUG
	fprintf(debf,"%send SGP_D_RDD_Gy\n", isp);
	indnt_dec();
#endif
}

/* obsolete: to be deleted
void SGP_D_RDD_GyS(		long*	n,
		float*	r_m,
		float*	E_MeV_u,
		int*	particle_no,
		char**	material_name,
		float*	parameter,
		float*	D_RDD_Gy){

	////	fprintf(debf,"%sbegin SGP_D_RDD_GyS\n");
	////	fprintf(debf,"%sn = %ld, material_name = %s, parameter = %e\n", *n, *material_name, *parameter);
	////	long i;
	////	for( i = 0 ; i < *n ; i++){
	////		fprintf(debf,"%sr_m[%ld]=%e\n", i , r_m[i]);
	////		fprintf(debf,"%sE_MeV_u[%ld]=%e\n", i , E_MeV_u[i]);
	////		fprintf(debf,"%sparticle_no[%ld]=%d\n", i , particle_no[i]);
	////	}
	//
	//	long particle_no_long = (long)(*particle_no);
	//
	//	SGP_D_RDD_Gy	(	n,
	//						r_m,
	//						E_MeV_u,
	//						&particle_no_long,
	//						*material_name,
	//						parameter,
	//						D_RDD_Gy);
	//
	//	fprintf(debf,"%send SGP_D_RDD_GyS\n");

};
*/



void SGP_r_RDD_m	(	long*	n,
		float*	D_RDD_Gy,
		/* radial dose distribution model */
		long*	rdd_model,       /* */
		long*	n_rdd_parameter, /* number of rdd parameters */
		float*	rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
		/* electron range model */
		long*	er_model,
		long*	n_er_parameter,
		float*	er_parameter,
		float*	r_RDD_m)
/*void SGP_r_RDD_m	(	long*	n,
		float*	D_Gy,
		float*	E_MeV_u,
		long*	particle_no,
		char*	material_name,
		float*	parameter,
		float*	r_RDD_m)*/
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


		long i;

		if( *rdd_model == RDD_Test){

	#ifdef _DEBUG
			fprintf(debf,"%sRDD_Test, rdd_parameter[0] = %g\n", isp, rdd_parameter[0]);
	#endif

			long n_f1_parameters = 1;
			float*	f1_parameters = (float*)calloc(n_f1_parameters, sizeof(float));

			SGP_RDD_f1_parameters(	n,
					rdd_model,       /* model: 0 - test, 1 - Katz, 2 - LEM*/
					n_rdd_parameter, /* number of rdd parameters */
					rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
					/* electron range model */
					er_model,
					n_er_parameter,
					er_parameter,
					/* calculated parameters */
					&n_f1_parameters,
					f1_parameters);

			/* rewriting output of SGP_RDD_f1_parameters function
			 * to more convenient format
			 */
			float*	factor = (float*)calloc(*n, sizeof(float));
			for( i = 0 ; i < *n ; i++){
				factor[i] = f1_parameters[0];
			}
			free(f1_parameters);

			/*
			 * calculation of radius from dose (to be revised!)
			 */
			for (i = 0; i < *n; i++){
				r_RDD_m[i]		=	0.0f;
				if (D_RDD_Gy[i] > 0.0f){
					r_RDD_m[i]		=	666.6f;}
			}

		}

		if( *rdd_model == RDD_Geiss){

			long n_f1_parameters = 8;
			float*	f1_parameters = (float*)calloc(n_f1_parameters, sizeof(float));

			SGP_RDD_f1_parameters(	n,
					rdd_model,       /* model: 0 - test, 1 - Katz, 2 - LEM*/
					n_rdd_parameter, /* number of rdd parameters */
					rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
					/* electron range model */
					er_model,
					n_er_parameter,
					er_parameter,
					/* calculated parameters */
					&n_f1_parameters,
					f1_parameters);

			/* rewriting output of SGP_RDD_f1_parameters function
			 * to more convenient format
			 */
			float*	LET_MeV_cm2_g				=	(float*)calloc(*n, sizeof(float));
			float*	r_min_m						=	(float*)calloc(*n, sizeof(float));
			float*	r_max_m						=	(float*)calloc(*n, sizeof(float));
			float*	d_min_Gy					=	(float*)calloc(*n, sizeof(float));
			float*	d_max_Gy					=	(float*)calloc(*n, sizeof(float));
			float*	k_Gy						=	(float*)calloc(*n, sizeof(float));
			float*	single_impact_fluence_cm2	=	(float*)calloc(*n, sizeof(float));
			float*	single_impact_dose_Gy		=	(float*)calloc(*n, sizeof(float));
			for( i = 0 ; i < *n ; i++){
				LET_MeV_cm2_g[i] = f1_parameters[0];
				r_min_m[i] = f1_parameters[1];
				r_max_m[i] = f1_parameters[2];
				d_min_Gy[i] = f1_parameters[3];
				d_max_Gy[i] = f1_parameters[4];
				k_Gy[i] = f1_parameters[5];
				single_impact_fluence_cm2[i] = f1_parameters[6];
				single_impact_dose_Gy[i] = f1_parameters[7];
			}
			free(f1_parameters);

			for (i = 0; i < *n; i++){
	#ifdef _DEBUG
				fprintf(debf,"%sD_RDD_Gy[%ld] = %g\n", isp, i , D_RDD_Gy[i]);
				fprintf(debf,"%sk_Gy[%ld] = %g\n", isp, i , k_Gy[i]);
	#endif

				// if D is outside the definition, return -1
				r_RDD_m[i]		=	-1.0f;

				if ((d_min_Gy[i] <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= d_max_Gy[i])){
					r_RDD_m[i]		=	r_min_m[i] * (float)sqrt(k_Gy[i] / D_RDD_Gy[i]);}
			}

			free(LET_MeV_cm2_g);
			free(r_min_m);
			free(r_max_m);
			free(d_min_Gy);
			free(d_max_Gy);
			free(k_Gy);
			free(single_impact_fluence_cm2);
			free(single_impact_dose_Gy);


		}

	#ifdef _DEBUG
		fprintf(debf,"%send SGP_r_RDD_m\n", isp);
		indnt_dec();
	#endif
}

void SGP_r_RDD_m_old	(	long*	n,
				float*	D_Gy,
				float*	E_MeV_u,
				long*	particle_no,
				char*	material_name,
				float*	parameter,
				float*	r_RDD_m)
		{
			long		n_particle	=	1;
			float		LET_MeV_cm2_g;
			float		r_min_m;
			float		r_max_m;
			float		d_min_Gy;
			float		d_max_Gy;
			float		k_Gy;
			float		single_impact_fluence_cm2;
			float		single_impact_dose_Gy;

			SGP_RDD_f1_parameters_Geiss(	&n_particle,
					E_MeV_u,
					particle_no,
					material_name,
					parameter,
					&LET_MeV_cm2_g,
					&r_min_m,
					&r_max_m,
					&d_min_Gy,
					&d_max_Gy,
					&k_Gy,
					&single_impact_fluence_cm2,
					&single_impact_dose_Gy);

			long	i;
			for (i = 0; i < *n; i++){
				// if D is outside the definition, return -1
				r_RDD_m[i]		=	-1.0f;

				if ((d_min_Gy <= D_Gy[i]) && (D_Gy[i] <= d_max_Gy)){
					r_RDD_m[i]		=	r_min_m * (float)sqrt(k_Gy/D_Gy[i]);}
			}
		}

#endif // SGP_RDD_H_
