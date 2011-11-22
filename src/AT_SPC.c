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

void endian_swap2(unsigned short* x)
{

  *x = (*x>>8) | (*x<<8);
}


void endian_swap4(unsigned int* x)
{
  *x = (*x>>24) |
    ((*x<<8) & 0x00FF0000) |
    ((*x>>8) & 0x0000FF00) |
    (*x<<24);
}


// __int64 for MSVC, "int int" for gcc
void endian_swap8(uint64_t * x)
{
  *x = (*x>>56) |
    ((*x<<40) & 0x00FF000000000000) |
    ((*x<<24) & 0x0000FF0000000000) |
    ((*x<<8)  & 0x000000FF00000000) |
    ((*x>>8)  & 0x00000000FF000000) |
    ((*x>>24) & 0x0000000000FF0000) |
    ((*x>>40) & 0x000000000000FF00) |
    (*x<<56);
}


void readStruct(FILE *fp,
		struct STRPSPCBTAG * str,
		bool switchEndian)
{
  size_t n_read_elements = fread(str, sizeof(struct STRPSPCBTAG), 1, fp);
  if( n_read_elements != 1)
	  printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
  if (switchEndian){
      endian_swap4(&(*str).ulLen);
      endian_swap4(&(*str).ulTag);
  }
}


int AT_SPC_get_size_from_filename(const char filename[FILE_NAME_NCHAR])
{
	FILE* fp;
	fp = fopen(filename, "rb");
    int size = AT_SPC_get_size(fp);
    fclose(fp);

	return(size);
}


int AT_SPC_get_size(FILE *fp){
	char string[80];
	bool switchEndian=false;
	int size=0;
	struct STRPSPCBTAG  filetype;
	readStruct(fp,&filetype,false);
	//Skip structure. Just going to assume the tags are right here.
    size_t n_read_elements = fread(string, sizeof(char), 80, fp);
    if( n_read_elements != 80)
    	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 80);
	if (strcmp(string,"SPCM")==0)
		switchEndian=true;
	fseek (fp ,  88*2, SEEK_CUR ); //Skip Fileversion and filedate

	struct STRPSPCBTAG targetName;
	readStruct(fp,&targetName,switchEndian);
	fseek (fp ,  targetName.ulLen, SEEK_CUR );
	struct STRPSPCBTAG projectileName;
	readStruct(fp,&projectileName,switchEndian);
	fseek (fp , projectileName.ulLen, SEEK_CUR );
	fseek(fp,16*3, SEEK_CUR); //Skip beam energy, peak position and normalization

	struct STRPSPCBTAG  nStepsTag;
	uint64_t nSteps;
	readStruct(fp,&nStepsTag,switchEndian);
	n_read_elements = fread(&nSteps,sizeof(uint64_t),1,fp);
    if( n_read_elements != 1)
    	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
	if (switchEndian) endian_swap8(&nSteps);
	uint64_t binsize;
	int i;
	for (i=0; i < nSteps; i++){
		fseek(fp,16*2, SEEK_CUR); //Skip depth and normalization
		struct STRPSPCBTAG  nSpeciesTag;
		uint64_t nSpecies;
		readStruct(fp,&nSpeciesTag,switchEndian);
		if (nSpeciesTag.ulTag != TRPSPCDTAG_NS){
			printf("Species tag corrupt: %u \n",nSpeciesTag.ulTag);
			return 0;
		}
		n_read_elements = fread(&nSpecies,sizeof(int64_t),1,fp);
	    if( n_read_elements != 1)
	    	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
		if (switchEndian) endian_swap8(&nSpecies);
		//    printf("Species: %u \n",nSpecies);
		int j;
		for (j=0; j < nSpecies; j++){
			struct STRPSPCBTAG  zA;
			readStruct(fp,&zA,switchEndian);
			if (zA.ulTag != TRPSPCDTAG_S){
				printf("ZA tag corrupt: %u \n",zA.ulTag);
				return 0;
			}
			fseek(fp,8*3,SEEK_CUR);

			struct STRPSPCBTAG  cum;
			readStruct(fp,&cum,switchEndian);
			if (cum.ulTag != TRPSPCDTAG_CUM){
				printf("Cum tag corrupt: %u \n",cum.ulTag);
				return 0;
			}
			fseek(fp,8,SEEK_CUR);
			struct STRPSPCBTAG  nC;
			readStruct(fp,&nC,switchEndian);
			if (nC.ulTag != TRPSPCDTAG_NC){
				printf("Cum tag corrupt: %u \n",cum.ulTag);
				return 0;
			}
			fseek(fp,8,SEEK_CUR);
			struct STRPSPCBTAG  nBins;
			readStruct(fp,&nBins,switchEndian);
			if (nBins.ulTag != TRPSPCDTAG_NE){
				printf("Energy bin tag corrupt: %u \n",nBins.ulTag);
				return 0;
			}
			n_read_elements = fread(&binsize,sizeof(uint64_t),1,fp); //Read number of NE bins
		    if( n_read_elements != 1)
		    	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
			if (switchEndian) endian_swap8(&binsize);
			size += binsize;

			struct STRPSPCBTAG  energyBins;
			readStruct(fp,&energyBins,switchEndian);
			if (energyBins.ulTag==TRPSPCDTAG_EREF){ //If it's a reference, skip just 8 bytes
				fseek(fp,energyBins.ulLen,SEEK_CUR);
			}
			else if (energyBins.ulTag==TRPSPCDTAG_E){
				fseek(fp,energyBins.ulLen,SEEK_CUR);
			} //Else skip energybins
			else {
				printf("SPC file corrupt! \n");
				printf("%ld \n",ftell(fp));
				printf("%u \n",energyBins.ulTag);
				return 0;
			}
			struct STRPSPCBTAG  histoBins;
			readStruct(fp,&histoBins,switchEndian);
			if (histoBins.ulTag != TRPSPCDTAG_HISTO){
				printf("History bin tag corrupt: %u \n",histoBins.ulTag);
				return 0;
			}
			fseek(fp,histoBins.ulLen,SEEK_CUR);

			struct STRPSPCBTAG  runSum;
			readStruct(fp,&runSum,switchEndian);
			if (runSum.ulTag != TRPSPCDTAG_RUNNINGSUM){
				printf("History bin tag corrupt: %u \n",runSum.ulTag);
				return 0;
			}
			fseek(fp,runSum.ulLen,SEEK_CUR);
			//      printf("RunSum: %u %u \n",runSum.ulTag,runSum.ulLen);
		}

	}
	return size;
}


int AT_SPC_get_number_of_bytes_in_file(const char * filename){
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


int AT_SPC_fast_read_buffer(const char * filename, int content_size, int32_t * content){

	/* open the file */
	int fd = open(filename, O_RDONLY);
	if (fd == -1){
		printf("open");
		return EXIT_FAILURE;
	}

	/* obtain file size */
	struct stat sb;
	if (fstat(fd, &sb) == -1){
		printf("fstat");
		close(fd);
		return EXIT_FAILURE;
	}

	/* check output array size */
	size_t length = sb.st_size;
	if( length != content_size * sizeof(int32_t)){
		printf("content has wrong size\n");
		close(fd);
		return EXIT_FAILURE;
	}
	close(fd);

	/* read whole file into addr pointer */
	const char * mode = "r";
	FILE * fs = fopen( filename , mode);
	size_t no_items_read = fread(content , sizeof(int32_t), length, fs);

	if (!feof(fs))
		printf("error, expected to be at the end of file %s, read %d items\n", filename, no_items_read);

	fclose(fs);

	return EXIT_SUCCESS;
}


void decomposeStructIntoString( const int32_t content[], char * string, int * length ){
	(*length) = content[1];
	string = (char*)calloc(sizeof(char),*length); // TODO move allocation outside
	memcpy( string, content+2, *length);
}


void decomposeStructIntoDouble( const int32_t content[], double * value, int * length ){
	(*length) = content[1];
	memcpy( value, content+2, *length);
}


void decomposeStructIntoInteger( const int32_t content[], uint64_t * value, int * length ){
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


int AT_SPC_decompose_size(const int content_size, int32_t content_orig[]){
	int size=0;
	int length=0;

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


int AT_SPC_decompose_header(
		const int content_size,
		int32_t   content_orig[],
		double*   E_MeV_u,
		double*   peak_position_g_cm2,
		double*   normalisation){

	int length = 0;
	int32_t * content = content_orig;

	skipStruct(&content); // filetype
	skipStruct(&content); // fileversion
	skipStruct(&content); // filedate
	skipStruct(&content); // targname
	skipStruct(&content); // projname

	decomposeStructIntoDouble(content, E_MeV_u, &length);
	skipStruct(&content); // beamenergy

	decomposeStructIntoDouble(content, peak_position_g_cm2, &length);
	skipStruct(&content); // peak position

	decomposeStructIntoDouble(content, normalisation, &length);
	skipStruct(&content); // normalisation

	return EXIT_SUCCESS;
}



int AT_SPC_decompose_data(
		const int content_size,
		int32_t   content_orig[],
		int*      depth_step[],
		double*   depth_g_cm2[],
		double*   E_MeV_u[],
		double*   DE_MeV_u[],
		long*     particle_no[],
		double*   fluence_cm2[]){

	int index=0;
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
					printf("problem nE  = %llu, length = %d\n", nE, length);
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

int AT_SPC_read(const char filename[FILE_NAME_NCHAR],
		int*    depth_step[],
		double* depth_g_cm2[],
		double* E_MeV_u[],
		double* DE_MeV_u[],
		long*   particle_no[],
		double* fluence_cm2[]){

	FILE* fp;
	fp = fopen(filename, "rb");

	int total_n_bins = AT_SPC_get_size(fp);

	*depth_step    =    (int*)realloc(*depth_step,  total_n_bins * sizeof(int));
	*depth_g_cm2   = (double*)realloc(*depth_g_cm2, total_n_bins * sizeof(double));
	*E_MeV_u       = (double*)realloc(*E_MeV_u,     total_n_bins * sizeof(double));
	*DE_MeV_u      = (double*)realloc(*DE_MeV_u,    total_n_bins * sizeof(double));
	*particle_no   =   (long*)realloc(*particle_no, total_n_bins * sizeof(long));
	*fluence_cm2   = (double*)realloc(*fluence_cm2, total_n_bins * sizeof(double));

	int n_bins_read = AT_SPC_read_data(fp,
			total_n_bins,
			*depth_step,
			*depth_g_cm2,
			*E_MeV_u,
			*DE_MeV_u,
			*particle_no,
			*fluence_cm2);

	fclose(fp);

    return(n_bins_read);
}


int AT_SPC_read_data_from_filename( const char filename[FILE_NAME_NCHAR],
		int    n,
		int    depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		long   particle_no[],
		double fluence_cm2[])
{
	FILE* fp;
	fp = fopen(filename, "rb");

	int n_bins_read = AT_SPC_read_data(fp,
			n,
			depth_step,
			depth_g_cm2,
			E_MeV_u,
			DE_MeV_u,
			particle_no,
			fluence_cm2);

	fclose(fp);

	return(n_bins_read);
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
	int size = nb / sizeof(int32_t);

	int32_t * content = (int32_t*)calloc(sizeof(int32_t), size);
	AT_SPC_fast_read_buffer(filename, size, content);

	int res = AT_SPC_decompose_data(
			size,
			content,
			&depth_step,
			&depth_g_cm2,
			&E_MeV_u,
			&DE_MeV_u,
			&particle_no,
			&fluence_cm2);

	free(content);

	return res;
}


int AT_SPC_read_data(FILE* fp,
		const int n,
		int    depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		long   particle_no[],
		double fluence_cm2[]){

	char string[80];
	bool switchEndian=false;
	int size=0;
	unsigned int i=0;
	unsigned int k =0;
	int particles[1000];
	rewind(fp);
	struct STRPSPCBTAG  filetype;
	readStruct(fp,&filetype,false);
	//Skip structure. Just going to assume the tags are right here.
	size_t n_read_elements = fread(string, sizeof(char), 80, fp);
    if( n_read_elements != 80)
    	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 80);
	if (strcmp(string,"SPCM")==0)
		switchEndian=true;
	fseek (fp ,  88*2, SEEK_CUR ); //Skip Fileversion and filedate

	struct STRPSPCBTAG targetName;
	readStruct(fp,&targetName,switchEndian);
	fseek (fp ,  targetName.ulLen, SEEK_CUR );
	struct STRPSPCBTAG projectileName;
	readStruct(fp,&projectileName,switchEndian);
	fseek (fp , projectileName.ulLen, SEEK_CUR );
	fseek(fp,16*3, SEEK_CUR); //Skip beam energy, peak position and normalization

	struct STRPSPCBTAG  nStepsTag;
	uint64_t nSteps;
	readStruct(fp,&nStepsTag,switchEndian);
	n_read_elements = fread(&nSteps,sizeof(uint64_t),1,fp);
	if( n_read_elements != 1)
	   	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
	if (switchEndian) endian_swap8(&nSteps);

	uint64_t binsize;
	for (i=0; i < nSteps; i++){
		struct STRPSPCBTAG depthTag;
		readStruct(fp,&depthTag,switchEndian);
		double curDepth;
		n_read_elements = fread(&curDepth,sizeof(double),1,fp);
		if( n_read_elements != 1)
		   	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
		if (switchEndian) endian_swap8((uint64_t*)(&curDepth));

		fseek(fp,16, SEEK_CUR); //Skip normalization
		struct STRPSPCBTAG  nSpeciesTag;
		uint64_t nSpecies;
		readStruct(fp,&nSpeciesTag,switchEndian);
		if (nSpeciesTag.ulTag != TRPSPCDTAG_NS){
			printf("Species tag corrupt: %u \n",nSpeciesTag.ulTag);
			return 0;
		}
		n_read_elements = fread(&nSpecies,sizeof(uint64_t),1,fp);
		if( n_read_elements != 1)
		   	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
		if (switchEndian) endian_swap8(&nSpecies);

		int j;
		for (j=0; j < nSpecies; j++){
			struct STRPSPCBTAG  zA;
			readStruct(fp,&zA,switchEndian);
			if (zA.ulTag != TRPSPCDTAG_S){
				printf("ZA tag corrupt: %u \n",zA.ulTag);
				return 0;
			}
			fseek(fp,8*2,SEEK_CUR);//Skip double version of A and Z
			uint32_t z,a;
			n_read_elements = fread(&z,sizeof(uint32_t),1,fp);//Read A and Z values
			if( n_read_elements != 1)
			   	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
			n_read_elements = fread(&a,sizeof(uint32_t),1,fp);
			if( n_read_elements != 1)
			   	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
			fseek(fp,8*4,SEEK_CUR);//Skip Cum & nC
			struct STRPSPCBTAG  nBins;
			readStruct(fp,&nBins,switchEndian);
			if (nBins.ulTag != TRPSPCDTAG_NE){
				printf("Energy bin tag corrupt: %u \n",nBins.ulTag);
				return 0;
			}
			n_read_elements = fread(&binsize,sizeof(uint64_t),1,fp); //Read number of NE bins
			if( n_read_elements != 1)
			   	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
			if (switchEndian) endian_swap8(&binsize);
			size += binsize;

			struct STRPSPCBTAG  energyBins;
			readStruct(fp,&energyBins,switchEndian);
			if (energyBins.ulTag==TRPSPCDTAG_EREF){ //If it's a reference, copy the referred particle
				uint64_t pRef;
				n_read_elements = fread(&pRef,sizeof(uint64_t),1,fp);
				if( n_read_elements != 1)
				   	printf("Read %zu element(s), instead of %d !\n", n_read_elements, 1);
				if (pRef > k){
					printf("SPC file corrupt \n");
					return 0;
				}
				particles[j]=particles[pRef];
				int n;
				for (n=0; n < binsize; n++){ // Copy reference
					E_MeV_u[k+n]=E_MeV_u[particles[pRef]+n];
					DE_MeV_u[k+n]=DE_MeV_u[particles[pRef]+n];
				}
			} else if (energyBins.ulTag==TRPSPCDTAG_E){
				double tempBins[binsize];
				n_read_elements = fread(tempBins, sizeof(double), binsize+1, fp);
				if( n_read_elements != binsize+1)
				   	printf("Read %zu element(s), instead of %llu !\n", n_read_elements, binsize+1);
				int n;
				for (n=0; n < binsize; n++){
					E_MeV_u[k+n]=(tempBins[n+1]+tempBins[n])/2;
					DE_MeV_u[k+n]=(tempBins[n+1]-tempBins[n]);
				}
				particles[j]=k;

			} //Else skip energybins
			else {
				printf("SPC file corrupt! \n");
				printf("%ld \n",ftell(fp));
				printf("%u \n",energyBins.ulTag);
				return 0;
			}


			struct STRPSPCBTAG  histoBins;
			readStruct(fp,&histoBins,switchEndian);
			if (histoBins.ulTag != TRPSPCDTAG_HISTO){
				printf("History bin tag corrupt: %u \n",histoBins.ulTag);
				return 0;
			}
			n_read_elements = fread(&fluence_cm2[k],sizeof(double),binsize,fp);
			if( n_read_elements != binsize)
			   	printf("Read %zu element(s), instead of %llu !\n", n_read_elements, binsize);
			int n;
			for (n = 0; n < binsize; n++){
				depth_g_cm2[k+n]=curDepth;
				depth_step[k+n]=i;
				particle_no[k+n]=z*1000+a;
			}
			k = k+binsize;
			struct STRPSPCBTAG  runSum;
			readStruct(fp,&runSum,switchEndian);
			if (runSum.ulTag != TRPSPCDTAG_RUNNINGSUM){
				printf("History bin tag corrupt: %u \n",runSum.ulTag);
				return 0;
			}
			fseek(fp,runSum.ulLen,SEEK_CUR);
		}

	}
	return size;
}
