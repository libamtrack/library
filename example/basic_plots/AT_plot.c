/**
 * @brief Plot generator
 */

/*
 *    AT_plot.c
 *    ===================
 *
 *    Created on: 2010-07-02
 *    Creator: kongruencja
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
 *    along with AmTrack (file: copying.txt).
 *    If not, see <http://www.gnu.org/licenses/>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>

// Some headers are found in different places in Mac OS X
#ifdef __APPLE__
	#include <sys/malloc.h>
#else
	#include <malloc.h>
#endif

#include "AT_RDD.h"
#include "AT_StoppingPower.h"
#include "AT_PhysicsRoutines.h"
#include "AT_DataRange.h"
#include "AT_Version.h"

char * plottypes[] = {"ER","RDD","LET","CSDArange"};
int plottypes_nr = 4;


void plottype_usage(){
	int i;
	printf("Plottypes supported: ");
	for( i = 0 ; i < plottypes_nr ; i++){
		printf("%s ",plottypes[i]);
	}
	printf("\n");
}


void usage(int argc, char* argv[]){
	printf("Usage: %s --plottype PLOTTYPE [ --modeltype MODELTYPE --submodeltype SUBMODELTYPE --xmin XMIN --xmax XMAX --xlogscale --npoints NPOINTS --particle PARTICLE]\n", argv[0]);
	printf("Usage (short): %s -t PLOTTYPE [ -s MODELTYPE -y SUBMODELTYPE -n XMIN -x XMAX -l -m NPOINTS -p PARTICLE]\n", argv[0]);
	plottype_usage();
}


int main( int argc, char* argv[]){

	int c;

	char *plottype = NULL;
	double x_start = -1;
	double x_stop = -1;
	int number_of_points_on_x_axis = -1;

	bool xlogscale_flag = false;

	long particle_no = -1;
	char particle_name[PARTICLE_NAME_NCHAR];

	long modeltype = -1;
	long submodeltype = -1;

	/* decode user parameters passed to program */
	while (1){
		static struct option long_options[] = {
				{"plottype",  required_argument,    0, 't'},
				{"modeltype",  required_argument,    0, 's'},
				{"submodeltype",  required_argument,    0, 'y'},
				{"xmin",   required_argument,    0, 'n'},
				{"xmax",   required_argument,    0, 'x'},
				{"xlogscale", no_argument  , 0, 'l'},
				{"npoints", required_argument  , 0, 'm'},
				{"particle", required_argument  , 0, 'p'},
                {"version", no_argument, 0, 'V'},
				{0, 0, 0, 0}
		};

		int option_index = 0; /* getopt_long stores the option index here. */

		c = getopt_long (argc, argv, "t:n:x:lm:p:s:y:V:", long_options, &option_index);

		if (c == -1) /* Detect the end of the options. */
			break;

		char * test = 0;
		switch (c) {
		case 't':
			plottype = optarg;
			break;

		case 'n':
			x_start = strtod(optarg,&test);
			if( ((x_start == 0) && (*test != 0)) || (*test != '\0')){
				fprintf(stderr, "Error in decoding xmin (--xmin option)\n");
			}
			break;

		case 'x':
			x_stop = strtod(optarg,&test);
			if( ((x_stop == 0) && (*test != 0)) || (*test != '\0')){
				fprintf(stderr, "Error in decoding xmax (--xmax option)\n");
			}
			break;

		case 'l':
			xlogscale_flag = true;
			break;

		case 'm':
			number_of_points_on_x_axis = strtol(optarg,&test,10);
			if( ((number_of_points_on_x_axis == 0) && (*test != 0)) || (*test != '\0')){
				fprintf(stderr, "Error in decoding npoints (--npoints option)\n");
			}
			break;

		case 'p':
			strncpy(particle_name, optarg, PARTICLE_NAME_NCHAR-1);
			particle_no = AT_particle_no_from_particle_name_single( particle_name );
			if( particle_no < 0)
				fprintf(stderr, "Error in decoding particle type (--particle option)\n");
			break;

		case 's':
			modeltype = strtol(optarg,&test,10);
			if( ((modeltype == 0) && (*test != 0)) || (*test != '\0')){
				fprintf(stderr, "Error in decoding modeltype (--modeltype option)\n");
			}
			break;

		case 'y':
			submodeltype = strtol(optarg,&test,10);
			if( ((submodeltype == 0) && (*test != 0)) || (*test != '\0')){
				fprintf(stderr, "Error in decoding submodeltype (--submodeltype option)\n");
			}
			break;

        case 'V':
            printf("Version: %s\n", AT_VERSION_GIT);
            return EXIT_SUCCESS;

		case '?': /* getopt_long already printed an error message. */
			break;

		default:
			abort ();
			break;
		}
	}

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		fprintf(stderr, "non-option ARGV-elements: ");
		while (optind < argc)
			fprintf(stderr, "%s ", argv[optind++]);
		fprintf(stderr, "\n");
	}

	/* prepare input */
	const long material_no = Water_Liquid;  // water

	/* prepare default values */
	if( plottype == NULL){
		fprintf(stderr, "Please specify plottype using -t (--plottype) option\n");
		usage(argc,argv);
		exit(EXIT_FAILURE);
	}

	if( strcmp( plottype , "ER") == 0 ){
		char test[200];
		int status_code = -1;
		if( modeltype >= 0 )
			status_code = getERName( modeltype, test);
		if( status_code != AT_Success){
			fprintf(stderr,"Specify ER model using --modeltype (-s) option\n");
			int i;
			for( i = 0 ; i < AT_ER_Data.n ; i++ ){
				fprintf(stderr,"%d - %s\n", AT_ER_Data.ER_no[i], AT_ER_Data.ER_name[i]);
			}
			exit(EXIT_FAILURE);
		}
	} else if( strcmp( plottype , "RDD") == 0 ){
		char test[200];
		int status_code = -1;
		if( modeltype >= 0 )
			status_code = AT_RDD_name_from_number( modeltype, test);
		if( status_code != AT_Success){
			fprintf(stderr,"Specify RDD model using --modeltype (-s) option\n");
			int i;
			for( i = 0 ; i < AT_RDD_Data.n ; i++ ){
				fprintf(stderr,"%ld - %s\n", AT_RDD_Data.RDD_no[i], AT_RDD_Data.RDD_name[i]);
			}
			exit(EXIT_FAILURE);
		}
		status_code = -1;
		if( submodeltype >= 0 )
			status_code = getERName( submodeltype, test);
		if( status_code != AT_Success){
			fprintf(stderr,"Specify ER model using --submodeltype (-y) option\n");
			int i;
			for( i = 0 ; i < AT_ER_Data.n ; i++ ){
				fprintf(stderr,"%d - %s\n", AT_ER_Data.ER_no[i], AT_ER_Data.ER_name[i]);
			}
			exit(EXIT_FAILURE);
		}
	} else if( strcmp( plottype , "LET") == 0 ){
		int status_code = -1;
		if( modeltype >= 0 )
			//status_code = AT_stopping_power_source_model_name_from_number( modeltype, test);
			status_code = 0;
		if( status_code != AT_Success){
			fprintf(stderr,"Specify LET source data type using --modeltype (-s) option\n");
			int i;
			for( i = 0 ; i < STOPPING_POWER_SOURCE_N ; i++ ){
//				fprintf(stderr,"%ld - %s\n", AT_stopping_power_sources.stopping_power_source_no[i], AT_stopping_power_sources.stopping_power_source_name[i]);
			}
			exit(EXIT_FAILURE);
		}
	}

	if( particle_no < 0){
		particle_no = 6012;
		AT_particle_name_from_particle_no_single(particle_no, particle_name);
		fprintf(stderr, "particle not set, setting default value %s\n", particle_name);
	}


	if( strcmp( plottype , "ER") == 0 || strcmp( plottype , "LET") == 0 || strcmp( plottype , "CSDArange") == 0){
		if(x_start < 0){
			x_start = 0.1;	    fprintf(stderr, "xmin not set, setting default value %g [MeV]\n", x_start);
		}
		if(x_stop < 0){
			x_stop = 500.0;	    fprintf(stderr, "xmax not set, setting default value %g [MeV]\n", x_stop);
		}
	} else if( strcmp( plottype , "RDD") == 0){
		if(x_start < 0){
			x_start = 1e-12;	fprintf(stderr, "xmin not set, setting default value %g [m]\n", x_start);
		}
		if(x_stop < 0){
			x_stop = 1e-6;		fprintf(stderr, "xmax not set, setting default value %g [m]\n", x_stop);
		}
	} else {
		if(x_start < 0){
			x_start = 1e-12;	fprintf(stderr, "xmin not set, setting default value %g\n", x_start);
		}
		if(x_stop < 0){
			x_stop = 1e-6;		fprintf(stderr, "xmax not set, setting default value %g\n", x_stop);
		}
	}

	if( number_of_points_on_x_axis < 0 ){
		number_of_points_on_x_axis = 10;	fprintf(stderr,"npoints not set, setting default value %d\n", number_of_points_on_x_axis);
	}

	double * x = (double*)calloc(number_of_points_on_x_axis,sizeof(double));
	double * y = (double*)calloc(number_of_points_on_x_axis,sizeof(double));
	int i;
	x[0] = x_start;

	if( xlogscale_flag ){
		const double x_multiplier = exp(((log( x_stop ) - log( x_start)) / ((double)number_of_points_on_x_axis-1.)));
		for( i = 1 ; i < number_of_points_on_x_axis ; i++){
			x[i] = x[i-1] * x_multiplier;
		}
	} else {
		const double x_increment = (x_stop - x_start) / ((double)number_of_points_on_x_axis-1.);
		for( i = 1 ; i < number_of_points_on_x_axis ; i++){
			x[i] = x[i-1] + x_increment;
		}
	}

	/* calculate results */
	if( strcmp( plottype , "ER") == 0){
		long er_model = modeltype;
		char er_name[200];
		getERName( er_model, er_name);
		printf("#ER: %s\n", er_name);
		printf("#particle: %s (code: %ld)\n", particle_name, particle_no);
		printf("#ionenergy[MeV] range[m]\n");
		AT_max_electron_ranges_m( number_of_points_on_x_axis , x, material_no, er_model, y);
	} else if( strcmp( plottype , "RDD") == 0){
		long rdd_model = modeltype;
		char rdd_name[200];
		AT_RDD_name_from_number(rdd_model, rdd_name);
		long  index = AT_RDD_index_from_RDD_number( (const long)rdd_model );
		double rdd_parameter[RDD_MAX_NUMBER_OF_PARAMETERS];
		int i;
		for( i = 0 ; i < RDD_MAX_NUMBER_OF_PARAMETERS ; i++)
			rdd_parameter[i] = AT_RDD_Data.parameter_default[index][i];
		double E_MeV_u = 100.0;
		long er_model = submodeltype;
		char er_name[200];
		getERName( er_model, er_name);
		printf("#ER: %s\n", er_name);
		printf("#RDD: %s\n", rdd_name);
		for( i = 0 ; i < AT_RDD_number_of_parameters(rdd_model) ; i++)
			printf("# RDD parameter %d : %s = %g\n", i , AT_RDD_Data.parameter_name[index][i], rdd_parameter[i]);
		printf("#particle: %s (code: %ld)\n", particle_name, particle_no);
		long source_no = PSTAR;
		char source_name[STOPPING_POWER_SOURCE_NAME_LENGTH] = "JAJA!";
        AT_stopping_power_source_model_name_from_number( source_no, source_name);
		printf("#data source: %s (code: %ld)\n", source_name, source_no);
		printf("#E: %g [MeV/u]\n", E_MeV_u);
		printf("#r[m] D[Gy]\n");
		int status = AT_D_RDD_Gy( number_of_points_on_x_axis, x, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, source_no, y);
		if( status != AT_Success ){
			fprintf(stderr, "Incompatible ER model (%s) used with RDD model (%s)\n", er_name, rdd_name);
			free(x);
			free(y);
			exit(EXIT_FAILURE);
		}
	} else if( strcmp( plottype , "LET") == 0){
		printf("#LET vs primary ion energy\n");
		printf("#particle: %s (code: %ld)\n", particle_name, particle_no);
		long source_no = modeltype;
		char source_name[STOPPING_POWER_SOURCE_NAME_LENGTH] = "JAJA!";
		AT_stopping_power_source_model_name_from_number( source_no, source_name);
	    printf("#data source: %s (code: %ld)\n", source_name, source_no);
		printf("#E[MeV/u] LET[MeV/cm2g]\n");
		for( i = 0 ; i < number_of_points_on_x_axis ; i++){
			AT_Mass_Stopping_Power_with_no( source_no,
					1,
					&x[i],
					&particle_no,
					material_no,
					&y[i]);

		}
	}  else if( strcmp( plottype , "CSDArange") == 0){
		printf("#CSDArange vs primary ion energy\n");
		printf("#particle: %s (code: %ld)\n", particle_name, particle_no);
		printf("#E[MeV/u] CSDArange[m]\n");
		for( i = 0 ; i < number_of_points_on_x_axis ; i++){
			y[i] = AT_CSDA_range_g_cm2_single( x[i], 0.0, particle_no,material_no);
			y[i] *= AT_density_g_cm3_from_material_no(material_no);
		}
	}  else {
		printf("Plottype %s not supported\n", plottype);
		plottype_usage();
		free(x);
		free(y);
		exit(EXIT_FAILURE);
	}

	/* save to file (stdout) */
	for( i = 0 ; i < number_of_points_on_x_axis ; i++){
		printf("%g %g\n",x[i],y[i]);
	}

	free(x);
	free(y);

	return EXIT_SUCCESS;
}
