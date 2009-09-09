#ifndef SGPARTICLE_H_
#define SGPARTICLE_H_

/*
 * SGParticle.h
 *
 *  Created on: 28.07.2009
 *      Author: greilich
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
#include <gsl/gsl_integration.h>

extern int indent_counter;
extern char isp[];
extern FILE * debf;

/**
 * Computes the efficiency using the SPIFF algorithm
 *
 * @param	n						number of particle types in the mixed particle field (pointer to single variable)
 * @param	E_MeV_u					energy of particles in the mixed particle field (pointer to array of size n)
 * @param	particle_no				type of the particles in the mixed particle field (pointer to array of size n)
 * @see								SGP_Data.h for definition
 * @param	fluence_cm2				fluences for the given particles, doses in Gy if negative (pointer to array of size n)
 * @param	material_no				index number for detector material (pointer to single variable)
 * @see								SGP_Constants.h for definition
 * @param	RDD_model				index number for chosen radial dose distribution (pointer to single variable)
 * @param	RDD_parameters			parameters for chosen radial dose distribution (pointer to array of size depending on chosen model)
 * @see								SGP_Constants.h for definition
 * @param	ER_model				index number for chosen electron-range model (pointer to single variable)
 * @param	ER_parameters			parameters for chosen electron-range model (pointer to array of size depending on chosen model)
 * @see								SGP_Constants.h for definition
 * @param	gamma_model				index number for chosen gamma response (pointer to single variable)
 * @param	gamma_parameters		parameters for chosen gamma response (pointer to array of size depending on chosen model)
 * @see								SGP_Constants.h for definition
 * @param	N2						(algorithm specific) number of bins per factor of two in local dose array (pointer to single variable)
 * @param	fluence_factor			factor to scale the fluences given as "fluence_cm2" with (pointer to single variable)
 * @param	write_output			if true, a protocol is written to "SuccessiveConvolutions.txt" in the working directory (pointer to single variable)
 * @param	shrink_tails			(algorithm specific) if true, tails of the local dose distribution, contributing less than "shrink_tails_under" are cut (pointer to single variable)
 * @param	shrink_tails_under		(algorithm specific) limit for tail cutting in local dose distribution (pointer to single variable)
 * @param	adjust_N2				(algorithm specific) if true, "N2" will be increase if necessary at high fluence to ensure sufficient binning resolution
 * @param	results					pointer to array of size 10 to be allocated by the user which will be used to return the results
 * 			results[0]				efficiency				(algorithm independent) main result: 	particle response at dose D / gamma response at dose D
 *			results[1]				d_check					(algorithm independent) sanity check:	total dose (in Gy) as returned by the alogrithm
 *			results[2]				S_HCP					(algorithm independent) 				absolute particle response
 *	 		results[3]				S_gamma					(algorithm independent)					absolute gamma response
 *	 		results[4]				not used				(algorithm independent)
 *			results[5]				u						(algorithm specific)					mean number of tracks contributing to representative point
 *			results[6]				u_start					(algorithm specific)					low starting value for mean number of tracks, where linearisation is applied
 *			results[7]				n_convolutions			(algorithm specific)					number of convolutions performed
 *	 		results[8]				not used				(algorithm specific)
 *	 		results[9]				not used				(algorithm specific)
 * @return							none
 */
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

void SGP_efficiency_Katz(	long*	n,
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
							float*	results);

#endif /* SGPARTICLE_H_ */
