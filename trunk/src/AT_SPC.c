#include "AT_SPC.h"

void endian_swap2(unsigned short* x)
{

  *x = (*x>>8) |
    (*x<<8);
}

void endian_swap4(unsigned int* x)
{
  *x = (*x>>24) |
    ((*x<<8) & 0x00FF0000) |
    ((*x>>8) & 0x0000FF00) |
    (*x<<24);
}

// __int64 for MSVC, "long long" for gcc
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
  fread(str, sizeof(struct STRPSPCBTAG), 1, fp);
  if (switchEndian){
      endian_swap4(&(*str).ulLen);
      endian_swap4(&(*str).ulTag);
  }
}

int AT_SPC_get_size(FILE *fp){
  char string[80];
  bool switchEndian=false;
  unsigned int size=0;
  struct STRPSPCBTAG  filetype;
  readStruct(fp,&filetype,false);
  //Skip structure. Just going to assume the tags are right here.
  fread(string, sizeof(char), 80, fp);
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
  fread(&nSteps,sizeof(uint64_t),1,fp);
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
    fread(&nSpecies,sizeof(int64_t),1,fp);
    if (switchEndian) endian_swap8(&nSpecies);
    printf("Species: %u \n",nSpecies);
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
      fread(&binsize,sizeof(uint64_t),1,fp); //Read number of NE bins
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
	printf("%u \n",ftell(fp));
	printf("%u \n",energyBins.ulTag);
	return 0;
      }
      struct STRPSPCBTAG  histoBins;
      readStruct(fp,&histoBins,switchEndian);
      if (histoBins.ulTag != TRPSPCDTAG_HISTO){
	printf("Hisotry bin tag corrupt: %u \n",histoBins.ulTag);
	return 0;
      }
      fseek(fp,histoBins.ulLen,SEEK_CUR);

      struct STRPSPCBTAG  runSum;
      readStruct(fp,&runSum,switchEndian);
      if (runSum.ulTag != TRPSPCDTAG_RUNNINGSUM){
	printf("Hisotry bin tag corrupt: %u \n",runSum.ulTag);
	return 0;
      }
      fseek(fp,runSum.ulLen,SEEK_CUR);
      printf("RunSum: %u %u \n",runSum.ulTag,runSum.ulLen);
    }

  }
  return size;



}


int AT_SPC_read_data(FILE *fp,
		long depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		long particle_no[],
		double fluence_cm2[]){

  char string[80];
  bool switchEndian=false;
  unsigned int size=0;
  unsigned int i=0;
  unsigned int k =0;
  int particles[1000];
  rewind(fp);
  struct STRPSPCBTAG  filetype;
  readStruct(fp,&filetype,false);
  //Skip structure. Just going to assume the tags are right here.
  fread(string, sizeof(char), 80, fp);
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
  fread(&nSteps,sizeof(uint64_t),1,fp);
  if (switchEndian) endian_swap8(&nSteps);

  uint64_t binsize;
  for (i=0; i < nSteps; i++){
    struct STRPSPCBTAG depthTag;
    readStruct(fp,&depthTag,switchEndian);
    double curDepth;
    fread(&curDepth,sizeof(double),1,fp);
    if (switchEndian) endian_swap8((uint64_t*)(&curDepth));

    fseek(fp,16, SEEK_CUR); //Skip normalization
    struct STRPSPCBTAG  nSpeciesTag;
    uint64_t nSpecies;
    readStruct(fp,&nSpeciesTag,switchEndian);
    if (nSpeciesTag.ulTag != TRPSPCDTAG_NS){
      printf("Species tag corrupt: %u \n",nSpeciesTag.ulTag);
      return 0;
    }
    fread(&nSpecies,sizeof(uint64_t),1,fp);
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
      fread(&z,sizeof(uint32_t),1,fp);//Read A and Z values
      fread(&a,sizeof(uint32_t),1,fp);
      fseek(fp,8*4,SEEK_CUR);//Skip Cum & nC
      struct STRPSPCBTAG  nBins;
      readStruct(fp,&nBins,switchEndian);
      if (nBins.ulTag != TRPSPCDTAG_NE){
	printf("Energy bin tag corrupt: %u \n",nBins.ulTag);
	return 0;
      }
      fread(&binsize,sizeof(uint64_t),1,fp); //Read number of NE bins
      if (switchEndian) endian_swap8(&binsize);
      size += binsize;

      struct STRPSPCBTAG  energyBins;
      readStruct(fp,&energyBins,switchEndian);
      if (energyBins.ulTag==TRPSPCDTAG_EREF){ //If it's a reference, copy the referred particle
	uint64_t pRef;
	fread(&pRef,sizeof(uint64_t),1,fp);
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
	fread(tempBins, sizeof(double), binsize+1, fp);
	int n;
	for (n=0; n < binsize; n++){
	  E_MeV_u[k+n]=(tempBins[n+1]+tempBins[n])/2;
	  DE_MeV_u[k+n]=(tempBins[n+1]-tempBins[n]);
	}
	particles[j]=k;

      } //Else skip energybins
      else {
	printf("SPC file corrupt! \n");
	printf("%u \n",ftell(fp));
	printf("%u \n",energyBins.ulTag);
	return 0;
      }


      struct STRPSPCBTAG  histoBins;
      readStruct(fp,&histoBins,switchEndian);
      if (histoBins.ulTag != TRPSPCDTAG_HISTO){
	printf("Hisotry bin tag corrupt: %u \n",histoBins.ulTag);
	return 0;
      }
      fread(&fluence_cm2[k],sizeof(double),binsize,fp);
      int n;
      for (n = 0; n < binsize; n++){
	depth_g_cm2[k+n]=curDepth;
	depth_step[k+n]=i;
	particle_no[k+n]=z*1000000+a;
      }
      k = k+binsize;
      struct STRPSPCBTAG  runSum;
      readStruct(fp,&runSum,switchEndian);
      if (runSum.ulTag != TRPSPCDTAG_RUNNINGSUM){
	printf("Hisotry bin tag corrupt: %u \n",runSum.ulTag);
	return 0;
      }
      fseek(fp,runSum.ulLen,SEEK_CUR);
    }

  }
  return size;
}


//void main()
//{
//  FILE *fp;
//  //  fp=fopen("TRIP_GSI_12C.H2O.MeV27000.spc", "r");
//  fp=fopen("out.spc", "r");
//  int size = AT_SPC_get_size(fp);
//  printf("%u \n",size);
//  int particleNo[size];
//  int depthStep[size];
//  double energy[size];
//  double dEnergy[size];
//  double depth[size];
//  double fluence[size];
//
//  AT_SPC_read_data(fp,depthStep,particleNo,energy,dEnergy,depth,fluence);
//
//
//}


