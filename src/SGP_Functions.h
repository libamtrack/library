#ifndef SGP_FUNCTIONS_H_
#define SGP_FUNCTIONS_H_

#include <math.h>
#include <stdbool.h>
#include "SGP_Constants.h"
#include "SGP_Data.h"
#include "SGP_Utils.h"

///////////////////////////////////////////////////////////////////////
// DATA ACCESS ROUTINE EXPORT (MAINLY FOR DEBUGGING)
void SGP_LET_MeV_cm2_g(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						char*	material_name,
						float*	LET_MeV_cm2_g);

void SGP_LET_MeV_cm2_gS(	long*	n,
							float*	E_MeV_u,
							long*	particle_no,
							char**	material_name,
							float*	LET_MeV_cm2_g)
{
	SGP_LET_MeV_cm2_g(	n,
						E_MeV_u,
						particle_no,
						*material_name,
						LET_MeV_cm2_g);
};


void SGP_LET_keV_um(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						char*	material_name,
						float*	LET_keV_um);

void SGP_CSDA_range_g_cm2(	long*	n,
							float*	E_MeV_u,
							long*	particle_no,
							char*	material_name,
							float*	CSDA_range_g_cm2);

void SGP_E_MeV(	long*	n,
				float*	CSDA_range_g_cm2,
				long*	particle_no,
				char*	material_name,
				float*	E_MeV);

void SGP_Particle_Properties(	long*	n,
								char**	particle_name,
								/* return values*/
								long*	particle_index,
								char**	USRTRACK_name,
								char**	element_name,
								long*	Z,
								long*	A,
								float*	mass);

/* As S-Plus does not allow to pass character arrays
 * to subroutines this wrapping routine handles
 * S-Plus calls of SGP_Particle_Properties
 */
void SGP_Particle_PropertiesS(	char**	particle_name,
								/* return values*/
								long*	particle_index,
								char**	USRTRACK_name,
								char**	element_name,
								long*	Z,
								long*	A,
								float*	mass)
{
	long n	= 1;
	SGP_Particle_Properties(	&n,
								particle_name,
								/* return values*/
								particle_index,
								USRTRACK_name,
								element_name,
								Z,
								A,
								mass);
};

///////////////////////////////////////////////////////////////////////
// MISC CONV. ROUTINES
void SGP_beta_from_particle_index(	long*	n,
									float*	E_MeV_u,
									long*	particle_no,
									float*	beta);

void SGP_beta_from_mass(	long*	n,
							float*	E_MeV_u,
							float*	mass,
							float*	beta);

void SGP_effective_charge_from_particle_index(	long*	n,
												float*	E_MeV_u,
												long*	particle_no,
												float*	effective_charge);
void SGP_effective_charge_from_beta(	long*	n,
										float*	beta,
										long*	Z,
										float*	effective_charge);

void SGP_scaled_energy(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						float*	effective_charge);


void SGP_max_E_transfer_MeV(	long*	n,
								float*	E_MeV_u,
								long*	particle_no,
								float*	max_E_transfer_MeV);

void SGP_max_electron_range_m(	long*	n,
								float*	E_MeV_u,
								long*	particle_no,
								char*	material_name,
								float*	max_electron_range_m);

void SGP_Z_from_particle_index(	long*	n,								// should rather be SGP_ParticleProperties(particle.index)
								long*	particle_no,
								long*	Z);

///////////////////////////////////////////////////////////////////////
// Routines to access PARTICLE data
///////////////////////////////////////////////////////////////////////
void SGP_Particle_Properties(	long*	n,
								char**	particle_name,
								/* return values*/
								long*	particle_index,
								char**	USRTRACK_name,
								char**	element_name,
								long*	Z,
								long*	A,
								float*	mass)
{
	long*	match	=	(long*)calloc(*n, sizeof(long));
	pmatchc(	particle_name,
				n,
				SGP_Particle_Data.particle_name,
				&SGP_Particle_Data.n,
				match);

	long i;
	for(i = 0; i < *n; i++){
		particle_index[i]			= SGP_Particle_Data.particle_index[match[i]];
		strcpy(USRTRACK_name[i], SGP_Particle_Data.USRTRACK_name[match[i]]);
		strcpy(element_name[i], SGP_Particle_Data.element_name[match[i]]);
		Z[i]						= SGP_Particle_Data.Z[match[i]];
		A[i]						= SGP_Particle_Data.A[match[i]];
		mass[i]						= SGP_Particle_Data.mass[match[i]];
	}
}

///////////////////////////////////////////////////////////////////////
// Routines to access MATERIAL data
///////////////////////////////////////////////////////////////////////
void SGP_getMaterialData(		long*	n,
								char**	material_names,
								float*	density_g_cm3,
								float*	electron_density_m3,
								float*	I_eV,
								float*	alpha_g_cm2_MeV,
								float*	p_MeV,
								float*	m_g_cm2)
{
	long*	match	=	(long*)calloc(*n, sizeof(long));
	pmatchc(	material_names,
				n,
				SGP_Material_Data.material_name,
				&SGP_Material_Data.n,
				match);

	long i;
	for(i = 0; i < *n; i++){
		density_g_cm3[i]		= SGP_Material_Data.density_g_cm3[match[i]];
		electron_density_m3[i]	= SGP_Material_Data.electron_density_m3[match[i]];
		I_eV[i]					= SGP_Material_Data.I_eV[match[i]];
		alpha_g_cm2_MeV[i]		= SGP_Material_Data.alpha_g_cm2_MeV[match[i]];
		p_MeV[i]				= SGP_Material_Data.p_MeV[match[i]];
		m_g_cm2[i]				= SGP_Material_Data.m_g_cm2[match[i]];
	}

	free(match);
}

void SGP_getMaterialDataS(		char**	material_name,
								float*	density_g_cm3,
								float*	electron_density_m3,
								float*	I_eV,
								float*	alpha_g_cm2_MeV,
								float*	p_MeV,
								float*	m_g_cm2){
	long	n;
	n		= 1;
	SGP_getMaterialData(		&n,
								material_name,
								density_g_cm3,
								electron_density_m3,
								I_eV,
								alpha_g_cm2_MeV,
								p_MeV,
								m_g_cm2);
}

#define matchIt			long	n_mat	= 1;									\
						pmatchc(	&material_name,								\
									&n_mat,										\
									SGP_Material_Data.material_name,			\
									&SGP_Material_Data.n,						\
									&match);

void SGP_density_g_cm3(			long*	n,
								char*	material_name,
								float*	density_g_cm3)
{	long	match; matchIt;
	long	i;
	for(i = 0; i < *n; i++){
		density_g_cm3[i]		= SGP_Material_Data.density_g_cm3[match];}}

void SGP_density_g_cm3S(		char**	material_name,
								float*	density_g_cm3){
	long	n;
	n		= 1;
	SGP_density_g_cm3(	&n,
						*material_name,
						density_g_cm3);
}
void SGP_electron_density_m3(	long*	n,
								char*	material_name,
								float*	electron_density_m3)
{	long	match; matchIt;
	long	i;
	for(i = 0; i < *n; i++){
		electron_density_m3[i]	= SGP_Material_Data.electron_density_m3[match];}}

void SGP_electron_density_m3S(	char**	material_name,
								float*	electron_density_m3){
	long	n;
	n		= 1;
	SGP_electron_density_m3(	&n,
								*material_name,
								electron_density_m3);
}

void SGP_alpha_g_cm2_MeV(		long*	n,
								char*	material_name,
								float*	alpha_g_cm2_MeV)
{	long	match; matchIt;
	long	i;
	for(i = 0; i < *n; i++){
		alpha_g_cm2_MeV[i]		= SGP_Material_Data.alpha_g_cm2_MeV[match];}}

void SGP_p_MeV(					long*	n,
								char*	material_name,
								float*	p_MeV)
{	long	match; matchIt;
	long	i;
	for(i = 0; i < *n; i++){
		p_MeV[i]				= SGP_Material_Data.p_MeV[match];}}

void SGP_m_g_cm2(				long*	n,
								char*	material_name,
								float*	m_g_cm2)
{	long	match; matchIt;
	long	i;
	for(i = 0; i < *n; i++){
		m_g_cm2[i]				= SGP_Material_Data.m_g_cm2[match];}}



///////////////////////////////////////////////////////////////////////
// Routines to access PSTAR data
///////////////////////////////////////////////////////////////////////

void SGP_LET_MeV_cm2_g(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						char*	material_name,
						float*	LET_MeV_cm2_g)
{

//	printf("begin SGP_LET_MeV_cm2_g\n");
//
//	printf("n = %ld, material_name = %s\n", *n, material_name);
//	long ii;
//	for( ii = 0 ; ii < *n ; ii++){
//		printf("E_MeV_u[%d]=%e\n", ii , E_MeV_u[ii]);
//		printf("particle_no[%d]=%d\n", ii , particle_no[ii]);
//	}


//	printf("E_MeV_u=%e\n", E_MeV_u[0]);
//	printf("particle_no=%ld\n", particle_no[0]);

	// get scaled energies for all given particles and energies
	float*	sE	=	(float*)calloc(*n, sizeof(float));
	SGP_scaled_energy(	n,
						E_MeV_u,
						particle_no,
						sE);

//	printf("sE[0]=%e\n",  sE[0]);

	// get effective charge for all given particles and energies
	float*	eC	=	(float*)calloc(*n, sizeof(float));
	SGP_effective_charge_from_particle_index(	n,
												E_MeV_u,
												particle_no,
												eC);

	// get LET-data for given material:
	// first: find all those PSTAR entries that match the material name
	bool*	matches	=	(bool*)calloc(SGP_PSTAR_Data.n, sizeof(bool));
	matchc(	material_name,
			SGP_PSTAR_Data.material_name,
			&SGP_PSTAR_Data.n,
			matches);

	// how many are this?
	long	n_entries	= 0;
	long	i;
	for (i = 0; i < SGP_PSTAR_Data.n; i++){
		if (matches[i]) {n_entries++;}
	}

	// allocate vectors for extracted LET entries
	float*	E	=	(float*)calloc(n_entries, sizeof(float));
	float*	L	=	(float*)calloc(n_entries, sizeof(float));

	// and get the values
	long j = 0;
	for (i = 0; i < SGP_PSTAR_Data.n; i++){
		if(matches[i]){
			E[j]	= SGP_PSTAR_Data.kin_E_MeV[i];
			L[j]	= SGP_PSTAR_Data.stp_pow_tot_MeV_cm2_g[i];
			j++;
		}
	}

	long	n_pol			= 4 + 1;
	for (i = 0; i < *n; i++){

		// Get proton-LET for scaled energy from table E, L using 4th degree polynomial (n_pol - 1 = 2) interpolation
		float	proton_LET		= 0.0f;
		float	err_proton_LET	= 0.0f;		// dummy
		interp(		E,
					L,
					&n_entries,
					&n_pol,
					&sE[i],
					&proton_LET,
					&err_proton_LET);

		// and get LET by scaling with effective charge^2
		LET_MeV_cm2_g[i] = 	eC[i] * eC[i] * proton_LET;

	}

	free(E);
	free(L);
	free(eC);
	free(sE);
	free(matches);

//	printf("end SGP_LET_MeV_cm2_g\n");

}


void SGP_LET_keV_um(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						char*	material_name,
						float*	LET_keV_um)
{
	// Get mass-norm. LET
	SGP_LET_MeV_cm2_g(	n,
						E_MeV_u,
						particle_no,
						material_name,
						LET_keV_um);

	// Get density matching to material_name (only 1 name therefore n_mat = 1)
	long	n_mat	= 1;
	long	match;
	pmatchc(	&material_name,
				&n_mat,
				SGP_Material_Data.material_name,
				&SGP_Material_Data.n,
				&match);

	long	i;
	for (i = 0; i < *n; i++){
		LET_keV_um[i]	*=	SGP_Material_Data.density_g_cm3[match];
	}

}

void SGP_CSDA_range_g_cm2(	long*	n,
							float*	E_MeV_u,
							long*	particle_no,
							char*	material_name,
							float*	CSDA_range_g_cm2)
{
	// first: find all those PSTAR entries that match the material name
	bool*	matches	=	(bool*)calloc(SGP_PSTAR_Data.n, sizeof(bool));
	matchc(	material_name,
			SGP_PSTAR_Data.material_name,
			&SGP_PSTAR_Data.n,
			matches);

	// how many are this?
	long	n_entries	= 0;
	long	i;
	for (i = 0; i < SGP_PSTAR_Data.n; i++){
		if (matches[i]) {n_entries++;}
	}

	// allocate vectors for extracted CDSA range entries
	float*	E	=	(float*)calloc(n_entries, sizeof(float));
	float*	r	=	(float*)calloc(n_entries, sizeof(float));

	// and get the values
	long j = 0;
	for (i = 0; i < SGP_PSTAR_Data.n; i++){
		if(matches[i]){
			E[j]	= SGP_PSTAR_Data.kin_E_MeV[i];
			r[j]	= SGP_PSTAR_Data.range_cdsa_g_cm2[i];
			j++;
		}
	}

	long	n_pol			= 4 + 1;
	for (i = 0; i < *n; i++){

		// Get CSDA range from table E, r using quadratic (n_pol = 2) interpolation
		float	err_range	= 0.0f;		// dummy
		interp(		E,
					r,
					&n_entries,
					&n_pol,
					&E_MeV_u[i],
					&CSDA_range_g_cm2[i],
					&err_range);
	}

	free(E);
	free(r);
	free(matches);
}

void SGP_E_MeV(	long*	n,
				float*	CSDA_range_g_cm2,
				long*	particle_index,
				char*	material_name,
				float*	E_MeV)
{
	// first: find all those PSTAR entries that match the material name
	bool*	matches	=	(bool*)calloc(SGP_PSTAR_Data.n, sizeof(bool));
	matchc(	material_name,
			SGP_PSTAR_Data.material_name,
			&SGP_PSTAR_Data.n,
			matches);

	// how many are this?
	long	n_entries	= 0;
	long	i;
	for (i = 0; i < SGP_PSTAR_Data.n; i++){
		if (matches[i]) {n_entries++;}
	}

	// allocate vectors for extracted energy entries
	float*	r	=	(float*)calloc(n_entries, sizeof(float));
	float*	E	=	(float*)calloc(n_entries, sizeof(float));

	// and get the values
	long j = 0;
	for (i = 0; i < SGP_PSTAR_Data.n; i++){
		if(matches[i]){
			r[j]	= SGP_PSTAR_Data.range_cdsa_g_cm2[i];
			E[j]	= SGP_PSTAR_Data.kin_E_MeV[i];
			j++;
		}
	}

	long	n_pol			= 2;
	for (i = 0; i < *n; i++){

		// Get energy from table r, E using quadratic (n_pol = 2) interpolation
		float	err_range	= 0.0f;		// dummy
		interp(		r,
					E,
					&n_entries,
					&n_pol,
					&CSDA_range_g_cm2[i],
					&E_MeV[i],
					&err_range);
	}

	free(r);
	free(E);
	free(matches);
}


//////////////////////////////////////////////////
// MISC ROUTINES
//////////////////////////////////////////////////

void SGP_beta_from_particle_index(	long*	n,
									float*	E_MeV_u,
									long*	particle_index,
									float*	beta)
{
	// find look-up indices for A's for particle numbers in particle data
	long*	matches	=	(long*)calloc(*n, sizeof(long));
	float*	mass	=	(float*)calloc(*n, sizeof(float));

	pmatchi(	particle_index,
				n,
				SGP_Particle_Data.particle_index,
				&SGP_Particle_Data.n,
				matches);

	// loop over n to find beta for all given particles and energies
	long	i;
	for(i = 0; i < *n; i++){
		mass[i]	= SGP_Particle_Data.mass[matches[i]];}

	SGP_beta_from_mass(	n,
						E_MeV_u,
						mass,
						beta);
	free(mass);
	free(matches);
}

void SGP_beta_from_mass(	long*	n,
							float*	E_MeV_u,
							float*	mass,
							float*	beta)
{
	// loop over n to find beta for all given particles and energies
	long	i;
	for(i = 0; i < *n; i++){

		float	E_MeV		=	E_MeV_u[i] * mass[i];

		// Get rest energy
		float	E0_MeV		=	(float)proton_mass_MeV_c2 * mass[i];

		// Return relativistic speed
		beta[i]				=	(float)sqrt(1 - 1/((1 + E_MeV/E0_MeV)*(1 + E_MeV/E0_MeV)));
	}
}


void SGP_effective_charge_from_particle_index(	long*	n,
												float*	E_MeV_u,
												long*	particle_index,
												float*	effective_charge)
{
	// get relativistic speeds for all given particles and energies
	float*	beta	=	(float*)calloc(*n, sizeof(float));
	long*	Z		=	(long*)calloc(*n, sizeof(long));

	SGP_beta_from_particle_index(	n,
									E_MeV_u,
									particle_index,
									beta);

	// find look-up indices for Z's for particle numbers in particle data
	long*	matches	=	(long*)calloc(*n, sizeof(long));
	pmatchi(	particle_index,
				n,
				SGP_Particle_Data.particle_index,
				&SGP_Particle_Data.n,
				matches);

	long i;
	for (i = 0; i < *n; i++){
		Z[i]	=	SGP_Particle_Data.Z[matches[i]];
	}

	SGP_effective_charge_from_beta(	n,
									beta,
									Z,
									effective_charge);

	free(beta);
	free(Z);
	free(matches);
}

void SGP_effective_charge_from_beta(	long*	n,
										float*	beta,
										long*	Z,
										float*	effective_charge)
{
	// loop over n
	long	i;
	for (i = 0; i < *n; i++){
		// Return effective charge according to Barkas-Bethe-approximation (but not for protons!)
		if (Z[i]!=1){
			effective_charge[i]	= (float)(Z[i]) * (1 - (float)exp(-125.0f * beta[i] / (pow(Z[i], 2.0f/3.0f))));}
		else{
			effective_charge[i]	= (float)(Z[i]);}
	}
}



void SGP_scaled_energy(	long*	n,
						float*	E_MeV_u,
						long*	particle_index,
						float*	scaled_energy)
{
//	printf("begin SGP_scaled_energy\n");
//	printf("n = %ld\n", *n);
//	long ii;
//	for( ii = 0 ; ii < *n ; ii++){
//		printf("E_MeV_u[%ld]=%e\n", ii , E_MeV_u[ii]);
//		printf("particle_index[%ld]=%ld\n", ii , particle_index[ii]);
//	}

	// find look-up indices for A's for particle numbers in particle data
	long*	matches	=	(long*)calloc(*n, sizeof(long));
	pmatchi(	particle_index,
				n,
				SGP_Particle_Data.particle_index,
				&SGP_Particle_Data.n,
				matches);

	// loop over n to find beta for all given particles and energies
	long	i;
	for(i = 0; i < *n; i++){

//		printf("matches[0] = %d\n", matches[i]);

		float	E_MeV		=	E_MeV_u[i] * SGP_Particle_Data.mass[matches[i]];

		//printf("E_MeV = %d\n", E_MeV);

		// Get rest energy
		float	E0_MeV		=	(float)proton_mass_MeV_c2 * SGP_Particle_Data.mass[matches[i]];

		//printf("E0_MeV = %d\n", E0_MeV);

		// Return mass-scaled energy
		scaled_energy[i]	=	E_MeV / E0_MeV * proton_mass_MeV_c2;
	}

	free(matches);

	//printf("end SGP_scaled_energy\n");

}

void SGP_max_E_transfer_MeV(	long*	n,
								float*	E_MeV_u,
								long*	particle_index,
								float*	max_E_transfer_MeV)
{
	/* if E_MeV_u < 0:		use non-relativistic formula
	 * if E_MeV_u > 0:		use relativistic formula
	 */
	int*	relativistic		=	(int*)calloc(*n, sizeof(int));
	long	i;
	for (i = 0; i < *n; i++){
		if(E_MeV_u[i] >= 0){
			relativistic[i]			= 1;
		}else{
			relativistic[i]			= 0;
			E_MeV_u[i]				*= -1.0f;
		}
	}

	// get relativistic speeds for all given particles and energies
	float*	beta	=	(float*)calloc(*n, sizeof(float));
	SGP_beta_from_particle_index(	n,
									E_MeV_u,
									particle_index,
									beta);

	for (i = 0; i < *n; i++){
		if(relativistic[i] == 0){
			max_E_transfer_MeV[i]	=	4.0f * electron_mass_MeV_c2 / proton_mass_MeV_c2* E_MeV_u[i];
		}else{
			max_E_transfer_MeV[i]	=	2.0f * electron_mass_MeV_c2 * beta[i] * beta[i] / (1.0f - beta[i] * beta[i]);
		}
	}

	free(relativistic);
	free(beta);
}

// ELECTRON RANGE

// Models available:
// 0 - Test
// 1 - Chatterjee & Holley (1993)
// 2 - Butts & Katz (1967), Chunxiang et al., 1985, Walig�rski et al., 1986, Zhang et al., 1994 fitted to data from
//     Kobetich and Katz, 1968; Iskef et al., 1983; Kanter and Sternglass, 1962
// 3 - Kiefer & Straaten (1986)
// 4 - Gei�, 1998
// 5 - Scholz, 2001

void SGP_max_electron_range_m(	long*	n,
								float*	E_MeV_u,
								long*	particle_index,
								char*	material_name,
								float*	max_electron_range_m)
{

//	printf("begin SGP_max_electron_range_m\n");

	// Get density matching to material_name (only 1 name therefore n_mat = 1)
	long	n_mat	= 1;
	long	match;
	pmatchc(	&material_name,
				&n_mat,
				SGP_Material_Data.material_name,
				&SGP_Material_Data.n,
				&match);

	long	i;
	for (i = 0; i < *n; i++){
		float tmpE	= E_MeV_u[i];
		if (tmpE < 0) {tmpE *= -1.0f;}	// E can be set neg. if non-PSTAR are given --> use pos. value
		// Simple power law to get electron range in water (Scholz, Habil 2001)
		max_electron_range_m[i]		=	5.0e-2f * (float)pow(tmpE, 1.7f);
		max_electron_range_m[i]		/=	m_to_um;
/*		# Alternatively, use the one of Geiss (1997) as (wrongly = E is the particle energy!) ref. in Sawakuchi (2007)
		#r.max.cm		<-	4e-5 * E.MeV.u^1.5
		#r.max.m		<-	r.max.cm / m.to.cm
		*/
		// Scale maximum el. range with material density relative to water (1/rho)
		max_electron_range_m[i]		/= SGP_Material_Data.density_g_cm3[match];
//		printf("max_electron_range_m[%d] = %g\n", i , max_electron_range_m[i]);

	}

//	printf("end SGP_max_electron_range_m\n");

}

void SGP_Bohr_Energy_Straggling_g_cm2(	long*	n,
										char**	material_name,
										float*	dsE2dz)
{
	// Get Bohr's energy spread (Wilson, 1947, Phys Rev 71, 385)
	long	i;
	for (i = 0; i < *n; i++){
		float	electron_density_m3;
		long	n_dummy = 1;
		SGP_electron_density_m3(	&n_dummy,
									material_name[i],
									&electron_density_m3);

		double	tmp					=	e_C * e_C * e_C * e_C * electron_density_m3;
		tmp							/=	4.0 * pi * e0_F_m * e0_F_m;
		tmp							/=	MeV_to_J * MeV_to_J * m_to_cm;
		dsE2dz[i]					=	(float)tmp;
	}
}

void SGP_Z_from_particle_index(	long*	n,
								long*	particle_index,
								long*	Z)
{
	// find look-up indices for A's for particle numbers in particle data
	long*	matches	=	(long*)calloc(*n, sizeof(long));

	pmatchi(	particle_index,
				n,
				SGP_Particle_Data.particle_index,
				&SGP_Particle_Data.n,
				matches);

	// loop over n to find Z for all given particles and energies
	long	i;
	for(i = 0; i < *n; i++){
		Z[i]	= SGP_Particle_Data.Z[matches[i]];}

	free(matches);
}

#endif // SGP_FUNCTIONS_H_
