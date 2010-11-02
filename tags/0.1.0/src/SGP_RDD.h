#ifndef SGP_RDD_H_
#define SGP_RDD_H_


void SGP_RDD_f1_parameters(	long*			n,
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
								float*			single_impact_dose_Gy)
{
	SGP_RDD_f1_parameters(	n,
								E_MeV_u,
								particle_no,
								*material_name,
								parameter,
								LET_MeV_cm2_g,
								r_min_m,
								r_max_m,
								d_min_Gy,
								d_max_Gy,
								k_Gy,
								single_impact_fluence_cm2,
								single_impact_dose_Gy);
}

void SGP_D_RDD_Gy	(	long*	n,
						float*	r_m,
						float*	E_MeV_u,
						long*	particle_no,
						char*	material_name,
						float*	parameter,
						float*	D_RDD_Gy);

void SGP_D_RDD_GyS(	long*	n,
						float*	r_m,
						float*	E_MeV_u,
						long*	particle_no,
						char**	material_name,
						float*	parameter,
						float*	D_RDD_Gy){

	SGP_D_RDD_Gy	(	n,
						r_m,
						E_MeV_u,
						particle_no,
						*material_name,
						parameter,
						D_RDD_Gy);

};

void SGP_r_RDD_m	(	long*	n,
						float*	D_Gy,
						float*	E_MeV_u,
						long*	particle_no,
						char*	material_name,
						float*	parameter,
						float*	r_RDD_m);

void SGP_RDD_f1_parameters(	long*			n,
							float*			E_MeV_u,
							long*			particle_no,
							char*			material_name,
							float*			parameter,					// general formalism for later implementation of other RDDs
							float*			LET_MeV_cm2_g,				// here: n_parameter = 1, parameter[0] = r_min_m
							float*			r_min_m,
							float*			r_max_m,
							float*			d_min_Gy,
							float*			d_max_Gy,
							float*			k_Gy,
							float*			single_impact_fluence_cm2,
							float*			single_impact_dose_Gy)
{
	// Get maximum electron range
	SGP_max_electron_range_m(	n,
								E_MeV_u,
								particle_no,
								material_name,
								r_max_m);

	// Get the LET
	// Usually the routine looks up the LET value corresponding to
	// the E values given. This can, however, be incorrect in case of
	// dose-weighted mean values. There is therefore the possibility
	// to pass on LET values different from the PSTAR tables.
	// In that case, n NEGATIVE E values are given, followed by
	// n LET value, coerced in the E_MeV_u array

	if(E_MeV_u[0] > 0){									// Usual case: get LET from PSTAR table
		SGP_LET_MeV_cm2_g(	n,
						E_MeV_u,
						particle_no,
						material_name,
						LET_MeV_cm2_g);}
	else{												// Alternatively: use provided LET values
		long	i;
		for (i = 0; i < *n; i++){
			LET_MeV_cm2_g[i]	=	E_MeV_u[i + (*n)];
//			E_MeV_u[i]		=	-1.0f * E_MeV_u[i];
		}
	}


	// Loop over all particles and energies given
	long	i;
	for (i = 0; i < *n; i++){

		r_min_m[i]						=	parameter[0];
		if(r_min_m[i] >= r_max_m[i]){
			r_min_m[i]						=	r_max_m[i];}

		// Normalization to match with LET
		float	tmp						=	(float)(0.5f + log(r_max_m[i] / r_min_m[i]));
		tmp								*=	2.0f * pi * (r_min_m[i] * m_to_cm) * (r_min_m[i] * m_to_cm);

		k_Gy[i]							=	LET_MeV_cm2_g[i] * MeV_g_to_J_kg / tmp;

		tmp								=	1 / (pi * (r_max_m[i] * m_to_cm) * (r_max_m[i] * m_to_cm));
		single_impact_fluence_cm2[i]	=	tmp;
		single_impact_dose_Gy[i]		=	LET_MeV_cm2_g[i] * MeV_g_to_J_kg * tmp;

		d_max_Gy[i]						=	k_Gy[i];
		d_min_Gy[i]						=	k_Gy[i] * r_min_m[i] * r_min_m[i] / (r_max_m[i] * r_max_m[i]);
	}
}


void SGP_D_RDD_Gy	(	long*	n,
						float*	r_m,
						float*	E_MeV_u,
						long*	particle_no,
						char*	material_name,
						float*	parameter,
						float*	D_RDD_Gy)
{
	float*	LET_MeV_cm2_g				=	(float*)calloc(*n, sizeof(float));
	float*	r_min_m						=	(float*)calloc(*n, sizeof(float));
	float*	r_max_m						=	(float*)calloc(*n, sizeof(float));
	float*	d_min_Gy					=	(float*)calloc(*n, sizeof(float));
	float*	d_max_Gy					=	(float*)calloc(*n, sizeof(float));
	float*	k_Gy						=	(float*)calloc(*n, sizeof(float));
	float*	single_impact_fluence_cm2	=	(float*)calloc(*n, sizeof(float));
	float*	single_impact_dose_Gy		=	(float*)calloc(*n, sizeof(float));

	SGP_RDD_f1_parameters(	n,
								E_MeV_u,
								particle_no,
								material_name,
								parameter,
								LET_MeV_cm2_g,
								r_min_m,
								r_max_m,
								d_min_Gy,
								d_max_Gy,
								k_Gy,
								single_impact_fluence_cm2,
								single_impact_dose_Gy);

	long	i;
	for (i = 0; i < *n; i++){
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


void SGP_r_RDD_m	(	long*	n,
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

	SGP_RDD_f1_parameters(	&n_particle,
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
