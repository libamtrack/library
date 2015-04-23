#include "AT_StoppingPowerDataFromFile.h"


int AT_FromFile_wrapper( const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		double mass_stopping_power_MeV_cm2_g[]){

	FILE *CSV;
  	CSV = fopen("FLUKA_DEDX_WATER_76.8eV.txt", "r");

	 if(NULL == CSV) {
		 fprintf(stderr, "Error opening stopping power data file ...\n");
		 return EXIT_FAILURE;
	 }


	 char 	line[1000], test[1000];

	 int  	i = 0, ZZ;
	 float	EE, SS;

	 int	Z[10000];
	 float	E[10000];
	 float  S[10000];

	 while (fgets(line, sizeof line, CSV)) {
	     if (*line == EOF) break;
	     if (isalpha(*line)) continue; /* ignore comment line */
	     	 sscanf(line, "%s", test);
	     	 if ((sscanf(line, "%d %e %e", &ZZ, &EE, &SS) != 3)&&(strlen(test) != 0)&&(strlen(test) != 1)) {
	         /* handle error */
				 printf("Error reading stopping power data (format correct?)...");
				 return EXIT_FAILURE;
			 } else {
			 /* handle variables */
				 Z[i] = ZZ;
				 E[i] = EE;
				 S[i++] = SS;
			 }
	 	 }

	 int* Z_array = (int*)malloc(i*sizeof(int));
	 float* E_array = (float*)malloc(i*sizeof(float));
	 float* S_array = (float*)malloc(i*sizeof(float));

	 i--;
	 memcpy(Z_array, Z, i*sizeof(int));
	 memcpy(E_array, E, i*sizeof(float));
	 memcpy(S_array, S, i*sizeof(float));

	 mass_stopping_power_MeV_cm2_g[0] = 666.6;

	 return -1;
}

