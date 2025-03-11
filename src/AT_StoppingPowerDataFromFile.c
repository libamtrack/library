#include "AT_StoppingPowerDataFromFile.h"


int AT_FromFile_wrapper( const long n,
		const double E_MeV_u[],
		const long particle_no[],
		const long material_no,
		const char info[],
		double mass_stopping_power_MeV_cm2_g[]){

#define MAX_NUMBER_OF_LINES 10000
/*
 * stopping_power_source not defined
 * trying to assign value from environment variable
 */
	if(NULL == info) {
		info = getenv("AT_stopping_power_source");
	}

	if(NULL == info) {
#ifndef NDEBUG
		printf("stopping power data file is not defined ...\nThe file name can be defined through the environment variable AT_stopping_power_source\n");
#endif
		return EXIT_FAILURE;
	}

	/**
	 * Read data from file
	 */
	FILE *CSV;
  	CSV = fopen(info, "r");

	 if(NULL == CSV) {
#ifndef NDEBUG
		 printf("Error opening stopping power data file ...\n");
#endif
		 return EXIT_FAILURE;
	 }


	 char   line[255], test[255];
	 long  	i = 0, ZZ;
	 double	EE, SS;

	 long*	  Z = (long*)calloc(MAX_NUMBER_OF_LINES, sizeof(long));
	 double*  E = (double*)calloc(MAX_NUMBER_OF_LINES, sizeof(double));
	 double*  S = (double*)calloc(MAX_NUMBER_OF_LINES, sizeof(double));;

	 while (fgets(line, sizeof line, CSV)) {
	     if (*line == EOF) break;
	     if (!isdigit(*line)){
	    	 continue; /* ignore comment line */
	     }else{
			 sscanf(line, "%s", test);
			 if ((sscanf(line, "%ld %le %le", &ZZ, &EE, &SS) != 3)&&(strlen(test) != 0)&&(strlen(test) != 1)) {
			 /* handle error */
#ifndef NDEBUG
				 printf("Error reading stopping power data (format correct?)...");
#endif
				 return EXIT_FAILURE;
			 } else {
			 /* handle variables */
				 Z[i] = ZZ;
				 E[i] = EE;
				 S[i++] = SS;
			 }
	     }
	 }
	 i--;
         fclose(CSV);
         CSV = NULL;

	 /**
	  * Get stopping power values
	  */
	 long j, k, m, n_matches = 0;
	 long curZ = 0, lastZ = -1;

	 double* curE = (double*)calloc(1, sizeof(double));
	 double* curS = (double*)calloc(1, sizeof(double));

	 bool* matches = (bool*)calloc(i, sizeof(bool));

	 for( j = 0; j < n; j++){
		 if((lastZ == -1)||(lastZ == curZ)){

			 curZ  = AT_Z_from_particle_no_single(particle_no[j]);
			 lastZ = curZ;

			 n_matches = is_element_int( curZ,
		 		Z,
				i,
		 		matches);

			 free(curE);
		 	 free(curS);

		 	 curE = (double*)calloc(n_matches, sizeof(double));
		 	 curS = (double*)calloc(n_matches, sizeof(double));

		 	 m = 0;
		 	 for(k = 0; k < i; k++){
		 		 if(matches[k] == true){
		 			 curE[m] = E[k];
		 			 curS[m++] = S[k];
		 		 }
		 	 }
		 }
		 mass_stopping_power_MeV_cm2_g[j] = AT_get_interpolated_y_from_input_table(
				 curE,
				 curS,
		 		 n_matches,
		 		 E_MeV_u[j]);
	 }

	 free(curE);
	 free(curS);

	 free(Z);
	 free(E);
	 free(S);

	 free(matches);
	 return AT_Success;
}

