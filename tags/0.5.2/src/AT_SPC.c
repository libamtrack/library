/**
 * @brief SPC reader
 */

/*
 *    AT_SPC.c
 *    ========
 *
 *    Copyright 2006, 2010 The libamtrack team
 *
 *    This file is part of the AmTrack program (libamtrack.sourceforge.net).
 *
 *    AmTrack is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    AmTrack is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    aint with AmTrack (file: copying.txt).
 *    If not, see <http://www.gnu.org/licenses/>
 */


#include "AT_SPC.h"


int AT_SPC_get_number_of_bytes_in_file(const char* filename){
	int fd = open(filename, O_RDONLY);
	if (fd == -1)
		return -1;
	struct stat sb;
	if (fstat(fd, &sb) == -1){
		close(fd);
		return -2;
	}
	int result = sb.st_size;
	close(fd);
	return result;
}


int AT_SPC_fast_read_buffer(const char* filename,
		int content_size,
		int32_t* content){

	/* open the file */
	int fd = open(filename, O_RDONLY);
	if (fd == -1){
#ifndef NDEBUG
		printf("open");
#endif
		return EXIT_FAILURE;
	}

	/* obtain file size */
	struct stat sb;
	if (fstat(fd, &sb) == -1){
#ifndef NDEBUG
		printf("fstat");
#endif
		close(fd);
		return EXIT_FAILURE;
	}

	/* check output array size */
	size_t length = sb.st_size;
	if( length != content_size * sizeof(int32_t)){
#ifndef NDEBUG
		printf("content has wrong size\n");
#endif
		close(fd);
		return EXIT_FAILURE;
	}
	close(fd);

	/* read whole file into addr pointer */
	const char * mode = "r";
	FILE * fs = fopen( filename , mode);
	size_t no_items_read = fread(content , sizeof(int32_t), content_size, fs);

#ifndef NDEBUG
	if (no_items_read != content_size)
		printf("error, expected to be at the end of file %s, read %lu items (out of %d)\n", filename, (unsigned long)no_items_read, content_size);
#endif

	fclose(fs);

	return EXIT_SUCCESS;
}


void decomposeStructIntoString( const int32_t content[],
		char * string,
		int * length ){
	(*length) = content[1];
	string = (char*)calloc(sizeof(char),*length); // TODO move allocation outside
	memcpy( string, content+2, *length);
}


void decomposeStructIntoDouble( const int32_t content[],
		double * value,
		int * length ){
	(*length) = content[1];
	memcpy( value, content+2, *length);
}


void decomposeStructIntoInteger( const int32_t content[],
		uint64_t * value,
		int * length ){
	(*length) = content[1];
	memcpy( value, content+2, *length);
}


int skipStruct(int32_t ** content){
	int length = (*content)[1];
	(*content) += (length / sizeof(int32_t)) + 2;
	return length;
}


int decomposeTag(const int32_t content[]){
	return content[0];
}


int decomposeLength(const int32_t content[]){
	return content[1];
}


int AT_SPC_decompose_size(const int content_size,
		int32_t content_orig[]){
	int size   = 0;
	int length = 0;

	int32_t * content = content_orig;

	skipStruct(&content); // filetype
	skipStruct(&content); // fileversion
	skipStruct(&content); // filedate
	skipStruct(&content); // targname
	skipStruct(&content); // projname
	skipStruct(&content); // beamenergy
	skipStruct(&content); // peak position
	skipStruct(&content); // normalisation

	uint64_t numberOfDepthSteps = 0;
	decomposeStructIntoInteger(content, &numberOfDepthSteps, &length);
	skipStruct(&content);

	uint64_t stepNumber = 0;
	for( stepNumber = 0 ; stepNumber < numberOfDepthSteps ; stepNumber++ ){

//		printf("-----> Step nr %d\n", stepNumber);
		skipStruct(&content); // depth
		skipStruct(&content); // normalization

		uint64_t numberOfSpecies;
		decomposeStructIntoInteger(content, &numberOfSpecies, &length);
		skipStruct(&content);

		uint64_t specieNumber = 0;
		for( specieNumber = 0 ; specieNumber < numberOfSpecies ; specieNumber++ ){

//			printf("--------------> Specie nr %d\n", specieNumber);
			skipStruct(&content); // NZ
			skipStruct(&content); // Cum
			skipStruct(&content); // nC

			uint64_t nE;
			decomposeStructIntoInteger(content, &nE, &length);
			skipStruct(&content);

			size += nE;

			skipStruct(&content); // energy bins
			skipStruct(&content); // histo
			skipStruct(&content); // running sum
		}
	}

	return size;
}


int AT_SPC_decompose_header( const int content_size,
		int32_t   content_orig[],
		double*   E_MeV_u,
		double*   peak_position_g_cm2,
		long*     particle_no,
		int*      material_no,
		double*   normalisation,
		int*      depth_steps_no){

	int length = 0;
	int32_t * content = content_orig;

	skipStruct(&content); // filetype
	skipStruct(&content); // fileversion
	skipStruct(&content); // filedate

	char * targetName = NULL;
	decomposeStructIntoString(content, targetName, &length);
	skipStruct(&content); // targname
	*material_no = Water_Liquid; // TODO write converter from targetName to material_no

	char * projectileName = NULL;
	decomposeStructIntoString(content, projectileName, &length);
	skipStruct(&content); // projname
	*particle_no = 6012;  // TODO write converter from projectileName to particle_no

	decomposeStructIntoDouble(content, E_MeV_u, &length);
	skipStruct(&content); // beamenergy

	decomposeStructIntoDouble(content, peak_position_g_cm2, &length);
	skipStruct(&content); // peak position

	decomposeStructIntoDouble(content, normalisation, &length);
	skipStruct(&content); // normalisation

    uint64_t numberOfDepthSteps = 0;
    decomposeStructIntoInteger(content, &numberOfDepthSteps, &length);
    skipStruct(&content);
    *depth_steps_no = numberOfDepthSteps;

	return EXIT_SUCCESS;
}



int AT_SPC_decompose_data( const int content_size,
		int32_t   content_orig[],
		int*      depth_step[],
		double*   depth_g_cm2[],
		double*   E_MeV_u[],
		double*   DE_MeV_u[],
		long*     particle_no[],
		double*   fluence_cm2[]){

	int index  = 0;
	int length = 0;

	int32_t * content = content_orig;

	char * filetype = NULL;
	decomposeStructIntoString(content, filetype, &length);
	bool switchEndian=false;
	if (strcmp(filetype,"SPCM")==0)
		switchEndian=true;
	free(filetype);
	skipStruct(&content);

	char * fileversion = NULL;
	decomposeStructIntoString(content, fileversion, &length);
	// TODO for future use
	free(fileversion);
	skipStruct(&content);

	char * filedate = NULL;
	decomposeStructIntoString(content, filedate, &length);
	// TODO for future use
	free(filedate);
	skipStruct(&content);

	char * targname = NULL;
	decomposeStructIntoString(content, targname, &length);
	// TODO for future use
	free(targname);
	skipStruct(&content);

	char * projname = NULL;
	decomposeStructIntoString(content, projname, &length);
	// TODO for future use
	free(projname);
	skipStruct(&content);

	double beamEnergy = 0.0;
	decomposeStructIntoDouble(content, &beamEnergy, &length);
	// TODO for future use
	skipStruct(&content);

	double peakPosition = 0.0;
	decomposeStructIntoDouble(content, &peakPosition, &length);
	// TODO for future use
	skipStruct(&content);

	double normalization = 0.0;
	decomposeStructIntoDouble(content, &normalization, &length);
	// TODO for future use
	skipStruct(&content);

	uint64_t numberOfDepthSteps = 0;
	decomposeStructIntoInteger(content, &numberOfDepthSteps, &length);
	skipStruct(&content);

	uint64_t stepNumber = 0;

	for( stepNumber = 0 ; stepNumber < numberOfDepthSteps ; stepNumber++){

		//printf("-----> Step nr %llu\n", stepNumber);

		double depth = 0.0;
		decomposeStructIntoDouble(content, &depth, &length);
		skipStruct(&content);

		double normalizationStep = 0.0;
		decomposeStructIntoDouble(content, &normalizationStep, &length);
		// TODO for future use
		skipStruct(&content);

		uint64_t numberOfSpecies = 0;
		decomposeStructIntoInteger(content, &numberOfSpecies, &length);
		skipStruct(&content);

		uint64_t specieNumber = 0;

		double ** binPointersE = (double**)calloc(sizeof(double*), numberOfSpecies);
		double ** binPointersDE = (double**)calloc(sizeof(double*), numberOfSpecies);

		for( specieNumber = 0 ; specieNumber < numberOfSpecies ; specieNumber++){

			// TODO comment needed
			binPointersE[specieNumber] = (*E_MeV_u) + index;
			binPointersDE[specieNumber] = (*DE_MeV_u) + index;

			//printf("--------------> Specie nr %llu\n", specieNumber);

			// TODO comment needed
			double tmp[2];
			decomposeStructIntoDouble(content, &(tmp[0]), &length);
			skipStruct(&content);
			long particle_no_spc = AT_particle_no_from_Z_and_A_single((int)(tmp[0]),(int)(tmp[1]));

			double Cum = 0.0;
			decomposeStructIntoDouble(content, &Cum, &length);
			// TODO for future use
			skipStruct(&content);

			uint64_t nC = 0;
			decomposeStructIntoInteger(content, &nC, &length);
			// TODO for future use
			skipStruct(&content);

			uint64_t nE = 0;
			decomposeStructIntoInteger(content, &nE, &length);
			skipStruct(&content);

			int i = 0;
			for( i = index ; i < index + nE; i++){
				(*depth_step)[i] = stepNumber;
				(*depth_g_cm2)[i] = depth;
				(*particle_no)[i] = particle_no_spc;
			}

			if( decomposeTag(content) == TRPSPCDTAG_E ){
				length = decomposeLength(content);
				if( length != (nE + 1)*sizeof(double) ) {
#ifndef NDEBUG
						printf("problem nE  = %llu, length = %d\n", (long long unsigned int)nE, length);
#endif
					} else {
					double * tempBins = (double*)calloc(sizeof(double), nE+1);
					decomposeStructIntoDouble(content, tempBins, &length);
					for( i = index ; i < index + nE; i++){
						(*E_MeV_u)[i] = 0.5 * (tempBins[i+1-index] + tempBins[i-index]);
						(*DE_MeV_u)[i] = tempBins[i+1-index] - tempBins[i-index];
					}
					free(tempBins);
				}
			}
			if( decomposeTag(content) == TRPSPCDTAG_EREF ){
				uint64_t lSRef;
				decomposeStructIntoInteger(content, &lSRef, &length);
				memcpy( (*E_MeV_u)+index , binPointersE[lSRef], nE * sizeof(double));
				memcpy( (*DE_MeV_u)+index , binPointersDE[lSRef], nE * sizeof(double));
			}
			skipStruct(&content);

			decomposeStructIntoDouble(content, (*fluence_cm2) + index, &length);
			skipStruct(&content);

			skipStruct(&content); // running cumulated spectrum bin values
			// TODO for future use

			index += nE;
		}

		free(binPointersE);
		free(binPointersDE);

	}
	return EXIT_SUCCESS;
}


long AT_SPC_get_number_of_bins_from_filename_fast( const char filename[FILE_NAME_NCHAR] )
{

	int nb = AT_SPC_get_number_of_bytes_in_file( filename );
	if( nb <= 0){
		return nb-1;
	}

	int size = nb / sizeof(int32_t);

	int32_t * content = (int32_t*)calloc(sizeof(int32_t), size);
	int status = AT_SPC_fast_read_buffer(filename, size, content);

	int res;
	if( status == EXIT_FAILURE ){
		res = -1;
	} else {
		res = AT_SPC_decompose_size( size,
				content);
	}
	free(content);

	return res;
}


int AT_SPC_read_header_from_filename_fast( const char filename[FILE_NAME_NCHAR],
		double*   E_MeV_u,
		double*   peak_position_g_cm2,
		long*     particle_no,
		int*      material_no,
		double*   normalisation,
		int*      depth_steps_no)
{
	int nb = AT_SPC_get_number_of_bytes_in_file( filename );
	if( nb <= 0){
		return nb-1;
	}

	int size = nb / sizeof(int32_t);

	int32_t * content = (int32_t*)calloc(sizeof(int32_t), size);
	int status = AT_SPC_fast_read_buffer(filename, size, content);
	int res;
	if( status == EXIT_FAILURE ){
		res = -1;
	} else {
		res = AT_SPC_decompose_header(
				size,
				content,
				E_MeV_u,
				peak_position_g_cm2,
				particle_no,
				material_no,
				normalisation,
				depth_steps_no);
	}

	free(content);

	return res;
}


int AT_SPC_read_data_from_filename_fast( const char filename[FILE_NAME_NCHAR],
		int n,
		int    depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		long   particle_no[],
		double fluence_cm2[])
{
	int nb = AT_SPC_get_number_of_bytes_in_file( filename );
	if( nb <= 0){
		return nb-1;
	}

	int size = nb / sizeof(int32_t);

	int32_t * content = (int32_t*)calloc(sizeof(int32_t), size);
	int status = AT_SPC_fast_read_buffer(filename, size, content);
	int res;
	if( status == EXIT_FAILURE ){
		res = -1;
	} else {
		res = AT_SPC_decompose_data(
				size,
				content,
				&depth_step,
				&depth_g_cm2,
				&E_MeV_u,
				&DE_MeV_u,
				&particle_no,
				&fluence_cm2);
	}

	free(content);

	return res;
}


int AT_SPC_read_from_filename_fast( const char filename[FILE_NAME_NCHAR],
		int 	  n,
		double*   E_MeV_u_initial,
		double*   peak_position_g_cm2,
		long*     particle_no_initial,
		int*      material_no,
		double*   normalisation,
		int*      depth_steps_no,
		int       depth_step[],
		double    depth_g_cm2[],
		double    E_MeV_u[],
		double    DE_MeV_u[],
		long      particle_no[],
		double    fluence_cm2[]){

	AT_SPC_read_header_from_filename_fast( filename,
			E_MeV_u_initial,
			peak_position_g_cm2,
			particle_no_initial,
			material_no,
			normalisation,
			depth_steps_no);

	return (AT_SPC_read_data_from_filename_fast( filename,
			n,
			depth_step,
			depth_g_cm2,
			E_MeV_u,
			DE_MeV_u,
			particle_no,
			fluence_cm2));

}

