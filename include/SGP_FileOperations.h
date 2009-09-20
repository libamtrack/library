#ifndef SGP_FILEOPERATIONS_H_
#define SGP_FILEOPERATIONS_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "SGP_Utils.h"

extern int indent_counter;
extern char isp[];
extern FILE * debf;

void SGP_browseInput(	char*	fileName,
						// return:
						long*	nLines,
						long*	n_gamma_parameter);

void SGP_browseInputS(	char**	fileName,
						// return:
						long*	nLines,
						long*	n_gamma_parameter);


void SGP_readInput(	char*	fileName,
					long*	nLines,
					long*	n_gamma_parameter,
					// return:
					float*	E_MeV_u,
					long*	particle_no,
					float*	fluence_cm2,
					long*	slab_no,
					float*	parameter,
					long*	N2,
					char*	material_name,
					long*	n_slabs,
					long*	gamma_model,
					float*	gamma_parameter);

void SGP_readInputS(char**	fileName,
					long*	nLines,
					long*	n_gamma_parameter,
					// return:
					float*	E_MeV_u,
					long*	particle_no,
					float*	fluence_cm2,
					long*	slab_no,
					float*	parameter,
					long*	N2,
					char**	material_name,
					long*	n_slabs,
					long*	gamma_model,
					float*	gamma_parameter);

void SGP_browseSpectrum(	char*	fileName,
							// return:
							long*	nLines);

void SGP_browseSpectrumS(	char**	fileName,
							// return:
							long*	nLines);


void SGP_readSpectrum(	char*	fileName,
						long*	nLines,
						// return:
						float*	E_MeV_u,
						long*	particle_no,
						float*	fluence_cm2,
						long*	slab_no);

void SGP_readSpectrumS(	char**	fileName,
						long*	nLines,
						// return:
						float*	E_MeV_u,
						long*	particle_no,
						float*	fluence_cm2,
						long*	slab_no);


/////////////////////////////////////////////////////////////////////////////////

#endif // SGP_FILEOPERATIONS_H_
