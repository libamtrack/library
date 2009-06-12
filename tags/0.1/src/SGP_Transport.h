#ifndef SGP_TRANSPORT_H_
#define SGP_TRANSPORT_H_

#include "SGP_ParabolicCylinderFunction.h"
#include "SGP_Constants.h"

void SGP_BortfeldTransportProton(	float*	E_initial_MeV,
									float*	sE_initial_MeV,
									float*	fluence_initial_cm2,
									char*	plateau_dose_material_name,
									long*	n_shielding_slabs,
									float*	shielding_thickness_m,
									char**	shielding_material_name,
									long*	n_detector_slabs,
									float*	detector_thickness_m,
									char*	detector_material_name,
									/* return values */
									float*  plateau_dose_Gy,
									float*  plateau_dose_noNuc_Gy,
									float*	detector_z_cm,
									float*	detector_E_MeV,
									float*	detector_sE_MeV,
									float*	detector_fluence_cm2,
									float*	detector_dfluencedz_cm,
									float*	detector_LET_MeV_g_cm2,
									float*	detector_dose_Gy,
									float*	detector_dose_noNuc_Gy,
									float*	geom_factor);

void SGP_BortfeldTransportProtonL(	float*	E_initial_MeV,
									float*	sE_initial_MeV,
									float*	fluence_initial_cm2,
									char*	plateau_dose_material_name,
									long*	n_shielding_slabs,
									float*	shielding_thickness_m,
									char*	shielding_material_nameL,
									long*	n_detector_slabs,
									float*	detector_thickness_m,
									char*	detector_material_name,
									/* return values */
									float*  plateau_dose_Gy,
									float*  plateau_dose_noNuc_Gy,
									float*	detector_z_cm,
									float*	detector_E_MeV,
									float*	detector_sE_MeV,
									float*	detector_fluence_cm2,
									float*	detector_dfluencedz_cm,
									float*	detector_LET_MeV_g_cm2,
									float*	detector_dose_Gy,
									float*	detector_dose_noNuc_Gy,
									float*	geom_factor)
{
	// S/R and L cannot pass string arrays, therefore they send a string with ";" delimiters
	// that is taken apart here and fed in the actual routine
	char**	shielding_material_name;
	long	i;
	if(*n_shielding_slabs > 0){
		shielding_material_name		= (char**)calloc(*n_shielding_slabs, sizeof(char*));
		for (i = 0; i < *n_shielding_slabs; i++){
			shielding_material_name[i]	= (char*)calloc(255, sizeof(char));
		}

		char	sep[]	= ";";
		char*	token;
		token	=	strtok( shielding_material_nameL, sep );
		strcpy(shielding_material_name[0], token);
		for (i = 1; i < *n_shielding_slabs; i++){
			token						=	strtok( NULL, sep );
			strcpy(shielding_material_name[i], token);
		}
	}else{
		char*	dummy;
		shielding_material_name = &dummy;
	}

	SGP_BortfeldTransportProton(	E_initial_MeV,
									sE_initial_MeV,
									fluence_initial_cm2,
									plateau_dose_material_name,
									n_shielding_slabs,
									shielding_thickness_m,
									shielding_material_name,
									n_detector_slabs,
									detector_thickness_m,
									detector_material_name,
									/* return values */
									plateau_dose_Gy,
									plateau_dose_noNuc_Gy,
									detector_z_cm,
									detector_E_MeV,
									detector_sE_MeV,
									detector_fluence_cm2,
									detector_dfluencedz_cm,
									detector_LET_MeV_g_cm2,
									detector_dose_Gy,
									detector_dose_noNuc_Gy,
									geom_factor);

	if(*n_shielding_slabs > 0){
		for (i = 0; i < *n_shielding_slabs; i++){
			free(shielding_material_name[i]);}
		free(shielding_material_name);}
};


void SGP_BortfeldTransportProtonS(	float*	E_initial_MeV,
									float*	sE_initial_MeV,
									float*	fluence_initial_cm2,
									char**	plateau_dose_material_name,
									long*	n_shielding_slabs,
									float*	shielding_thickness_m,
									char**	shielding_material_name,
									long*	n_detector_slabs,
									float*	detector_thickness_m,
									char**	detector_material_name,
									/* return values */
									float*  plateau_dose_Gy,
									float*  plateau_dose_noNuc_Gy,
									float*	detector_z_cm,
									float*	detector_E_MeV,
									float*	detector_sE_MeV,
									float*	detector_fluence_cm2,
									float*	detector_dfluencedz_cm,
									float*	detector_LET_MeV_g_cm2,
									float*	detector_dose_Gy,
									float*	detector_dose_noNuc_Gy,
									float*	geom_factor)
{
	// S/R cannot pass string arrays, therefore n_shielding_slabs = 0 or 1
	SGP_BortfeldTransportProton(	E_initial_MeV,
									sE_initial_MeV,
									fluence_initial_cm2,
									*plateau_dose_material_name,
									n_shielding_slabs,
									shielding_thickness_m,
									shielding_material_name,
									n_detector_slabs,
									detector_thickness_m,
									*detector_material_name,
									/* return values */
									plateau_dose_Gy,
									plateau_dose_noNuc_Gy,
									detector_z_cm,
									detector_E_MeV,
									detector_sE_MeV,
									detector_fluence_cm2,
									detector_dfluencedz_cm,
									detector_LET_MeV_g_cm2,
									detector_dose_Gy,
									detector_dose_noNuc_Gy,
									geom_factor);

};

void	SGP_Dose_Gy(	float*	density_g_cm3,
						float*	fluence_cm2,
						float*	LET_MeV_g_cm2,
						float*	dfluencedz_cm,
						float*	E_MeV,
						float*	dose_Gy,
						float*	dose_noNuc_Gy)
{
	// Contribution of inelastic nuclear collision to dose (Berger, "Penetration of proton beam through water, 1993)
	float	gamma	=	0.6f;
	float	tmp		=	(*fluence_cm2) * (*LET_MeV_g_cm2);
	*dose_noNuc_Gy	=	tmp * MeV_to_J * g_cm3_to_kg_m3;
	tmp				+=	gamma * (*dfluencedz_cm) * (*E_MeV) * (*density_g_cm3);
	*dose_Gy		=	tmp * MeV_to_J * g_cm3_to_kg_m3;
}

void SGP_BortfeldTransportProton(	float*	E_initial_MeV,
									float*	sE_initial_MeV,
									float*	fluence_initial_cm2,
									char*	plateau_dose_material_name,
									long*	n_shielding_slabs,
									float*	shielding_thickness_m,
									char**	shielding_material_name,
									long*	n_detector_slabs,
									float*	detector_thickness_m,
									char*	detector_material_name,
									/* return values */
									float*  plateau_dose_Gy,
									float*  plateau_dose_noNuc_Gy,
									float*	detector_z_cm,
									float*	detector_E_MeV,
									float*	detector_sE_MeV,
									float*	detector_fluence_cm2,
									float*	detector_dfluencedz_cm,
									float*	detector_LET_MeV_g_cm2,
									float*	detector_dose_Gy,
									float*	detector_dose_noNuc_Gy,
									float*	geom_factor)
{
	/////////////////////////////////////////////////////////
	// The following code only works for protons
	// R(E) and S(E) are approximated by power laws
	// which is valid from 1 MeV - 200 MeV (Bortfeld, 1997)

	/////////////////////////////////////////////////////////
	// GET DATA FOR MATERIALS INVOLVED

	long	i;

	long	nMaterials				= *n_shielding_slabs + 2;				// max. one material for each shielding slab + plateau material + detector

	float*	density_g_cm3			= (float*)calloc(nMaterials, sizeof(float));
	float*	electron_density_m3		= (float*)calloc(nMaterials, sizeof(float));
	float*	alpha_g_cm2_MeV			= (float*)calloc(nMaterials, sizeof(float));
	float*	alpha_cm_MeV			= (float*)calloc(nMaterials, sizeof(float));
	float*	p_MeV					= (float*)calloc(nMaterials, sizeof(float));
	float*	m_g_cm2					= (float*)calloc(nMaterials, sizeof(float));
	float*	m_cm					= (float*)calloc(nMaterials, sizeof(float));

	float*	dsE2dz					= (float*)calloc(nMaterials, sizeof(float));

	char**	materials				= (char**)calloc(nMaterials, sizeof(char*));
	materials[0]					= plateau_dose_material_name;
	materials[1]					= detector_material_name;
	if (*n_shielding_slabs > 0){
		for (i = 2; i < nMaterials; i++){
			materials[i]					= shielding_material_name[i - 2];}
	}

	SGP_getMaterialData(	&nMaterials,
							materials,
							density_g_cm3,
							electron_density_m3,
							alpha_g_cm2_MeV,
							p_MeV,
							m_g_cm2);

	SGP_Bohr_Energy_Straggling_g_cm2(	&nMaterials,
										materials,
										dsE2dz);

	for (i = 0; i < nMaterials; i++){
		alpha_cm_MeV[i]		=	alpha_g_cm2_MeV[i] / density_g_cm3[i];
		m_cm[i]				=	m_g_cm2[i] / density_g_cm3[i];
	}

	/////////////////////////////////////////////////////////
	// PLATEAU DOSE

	float	R0_cm					= alpha_cm_MeV[0] * pow(*E_initial_MeV, p_MeV[0]);
	float	LET_MeV_g_cm2			= 1.0f / (alpha_cm_MeV[0] * p_MeV[0]) * pow(*E_initial_MeV, 1.0f - p_MeV[0]) / density_g_cm3[0];
	// if fluence < 0 then dose given, so compute fluence
	if(*fluence_initial_cm2 < 0){
		*fluence_initial_cm2			= -1.0f * (*fluence_initial_cm2) / (LET_MeV_g_cm2 * MeV_to_J * g_cm3_to_kg_m3);}
	float	dfluencedz				= (*fluence_initial_cm2) * m_cm[0] / (1 + m_cm[0] * R0_cm);


	SGP_Dose_Gy(	&density_g_cm3[0],
					fluence_initial_cm2,
					&LET_MeV_g_cm2,
					&dfluencedz,
					E_initial_MeV,
					// return
					plateau_dose_Gy,
					plateau_dose_noNuc_Gy);


	/////////////////////////////////////////////////////////
	// TRANSPORT THROUGH SHIELDING SLABS
	float	E_MeV			=	*E_initial_MeV;
	float	sE2_MeV2		=	(*sE_initial_MeV) * (*sE_initial_MeV);
	float	fluence_cm2		=	*fluence_initial_cm2;

	if (*n_shielding_slabs > 0){
		for (i = 0; i < *n_shielding_slabs; i++){
			// some conversions
			long	k				=	i + 2;					// pointer to material data in material list
			float	dz_cm			=	shielding_thickness_m[i] * m_to_cm;
			R0_cm					=	alpha_cm_MeV[k] * pow(E_MeV, p_MeV[k]);
			float	u				=	R0_cm - dz_cm;

			// Add energy straggling
			sE2_MeV2				=	dsE2dz[i] * dz_cm + sE2_MeV2 * sE2_MeV2;
//  Is the beam stopped in the slab?
//	ACCOUNT FOR STRAGGLING!
//  this approach is only valid if slabs are far from Bragg peak!
//	SGP_Funs(&z_cm, &range_cm, &sigma, &ni, &tmp);
			if (u <= 0){
				E_MeV			=	0.0f;
				fluence_cm2		=	0.0f;
				break;}
			//	If not, get energy and fluence after passing the shielding slab
			E_MeV			=	1.0f / pow(alpha_cm_MeV[k], 1.0f / p_MeV[k]) * pow(u, 1.0f / p_MeV[k]);
			fluence_cm2		*=	(1.0f + m_cm[k] * u) / (1.0f + m_cm[0] * R0_cm);
		}
	}

	/////////////////////////////////////////////////////////
	// TRANSPORT THROUGH DETECTOR

	float	dz_cm			=	(*detector_thickness_m / *n_detector_slabs) * m_to_cm;
	R0_cm					=	alpha_cm_MeV[1] * pow(E_MeV, p_MeV[1]);

	float	sz2_cm2			=	dsE2dz[1] *																// Bortfeld's s_mono
								p_MeV[1] * p_MeV[1] * pow(alpha_cm_MeV[1], 2.0 / p_MeV[1]) /
								(3.0 - 2.0 / p_MeV[1]) * pow(R0_cm, 3.0 - 2.0 / p_MeV[1]);

	float	sigma2_cm2		=	sz2_cm2 + sE2_MeV2	* pow(alpha_cm_MeV[1], 2)							// Add straggling in detector and from slabs
													* pow(p_MeV[1], 2)								// see Bortfeld, eq. 19
													* pow(E_MeV, 2 * p_MeV[1] - 2.0f);


	float	sigma_cm		=	sqrt(sigma2_cm2);

	*geom_factor			=	0.0f;

	for (i = 0; i < *n_detector_slabs; i++){

		detector_z_cm[i]			=	dz_cm * (i + 0.5f);

		float	tmp					=	0.0f;

		// E incl. straggling
		float	ni					=	1.0f / p_MeV[1];
		SGP_Funs(&detector_z_cm[i], &R0_cm, &sigma_cm, &ni, &tmp);
		detector_E_MeV[i]			=	1.0f / pow(alpha_cm_MeV[1], ni) * tmp;
		detector_sE_MeV[i]			=	sqrt( sE2_MeV2 + pow(detector_z_cm[i] * dsE2dz[1], 2));

		// L incl. straggling
		ni							=	1.0f / p_MeV[1] - 1.0f;
		SGP_Funs(&detector_z_cm[i], &R0_cm, &sigma_cm, &ni, &tmp);
		detector_LET_MeV_g_cm2[i]	=	1.0f / (p_MeV[1] * pow(alpha_cm_MeV[1], 1.0f / p_MeV[1])) * tmp / density_g_cm3[1];

		// F incl. straggling
		ni							=	1.0f;
		SGP_Funs(&detector_z_cm[i], &R0_cm, &sigma_cm, &ni, &tmp);
		detector_fluence_cm2[i]		=	fluence_cm2 *  (1.0f + m_cm[1] * R0_cm) * (1.0f + m_cm[1] * tmp);

		detector_dfluencedz_cm[i]	=	fluence_cm2 * m_cm[1] / (1.0f + m_cm[1] * R0_cm);

		SGP_Dose_Gy(	&density_g_cm3[1],
						&detector_fluence_cm2[i],
						&detector_LET_MeV_g_cm2[i],
						&detector_dfluencedz_cm[i],
						&detector_E_MeV[i],
						// return
						&detector_dose_Gy[i],
						&detector_dose_noNuc_Gy[i]);

		*geom_factor				+=	detector_dose_Gy[i];
	}

	*geom_factor				/=	(*n_detector_slabs) * detector_dose_Gy[0];
}



#endif // SGP_TRANSPORT_H_
