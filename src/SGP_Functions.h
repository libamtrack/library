#ifndef SGP_FUNCTIONS_H_
#define SGP_FUNCTIONS_H_

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "SGP_Constants.h"
#include "SGP_Data.h"
#include "SGP_Utils.h"

extern int indent_counter;
extern char isp[];
extern FILE * debf;

///////////////////////////////////////////////////////////////////////
// DATA ACCESS ROUTINE EXPORT (MAINLY FOR DEBUGGING)
void SGP_LET_MeV_cm2_g(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						long*	material_no,
						float*	LET_MeV_cm2_g);

void SGP_LET_keV_um(	long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						long*	material_no,
						float*	LET_keV_um);

void SGP_CSDA_range_g_cm2(	long*	n,
							float*	E_MeV_u,
							long*	particle_no,
							long*	material_no,
							float*	CSDA_range_g_cm2);

void SGP_E_MeV(	long*	n,
				float*	CSDA_range_g_cm2,
				long*	particle_no,
				long*	material_no,
				float*	E_MeV);

void SGP_Particle_Properties(	long*	particle_no,
								/* return values*/
								char**	particle_name,
								char**	USRTRACK_name,
								char**	element_name,
								long*	Z,
								long*	A,
								float*	mass);

void SGP_getMaterialData(		long*	n,
								long*	material_no,
								float*	density_g_cm3,
								float*	electron_density_m3,
								float*	I_eV,
								float*	alpha_g_cm2_MeV,
								float*	p_MeV,
								float*	m_g_cm2);

///////////////////////////////////////////////////////////////////////
// MISC CONV. ROUTINES
void SGP_beta_from_particle_no(	long*	n,
									float*	E_MeV_u,
									long*	particle_no,
									float*	beta);

void SGP_beta_from_mass(	long*	n,
							float*	E_MeV_u,
							float*	mass,
							float*	beta);

void SGP_Bohr_Energy_Straggling_g_cm2(	long*	n,
										char**	material_name,
										float*	dsE2dz);

void SGP_effective_charge_from_particle_no(	long*	n,
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
								long*	material_no,
								long*	er_model,
								float*	max_electron_range_m);

#ifdef _R
void SGP_max_electron_range_mS(	int*	n,
								float*	E_MeV_u,
								int*	particle_no,
								int*	material_no,
								int*   er_model,
								float*	max_electron_range_m);
#endif
#ifdef _S
void SGP_max_electron_range_mS(	long*	n,
								float*	E_MeV_u,
								long*	particle_no,
								long*	material_no,
								long*	er_model,
								float*	max_electron_range_m);
#endif

void SGP_Z_from_particle_no(	long*	n,								// should rather be SGP_ParticleProperties(particle.index)
								long*	particle_no,
								long*	Z);

#endif // SGP_FUNCTIONS_H_
