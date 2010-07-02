/**
 * @file
 * @brief Plot generator
 */

/*
 *    AT_plot.c
 *    ===================
 *
 *    Created on: 2010-07-02
 *    Author: grzanka
 *
 *    Copyright 2006, 2009 Steffen Greilich / the libamtrack team
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
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>

#include "AmTrack.h"

char * plottypes[] = {"ER","RDD","LET"};
int plottypes_nr = 3;


void plottype_usage(){
	int i;
	printf("Plottypes supported: ");
	for( i = 0 ; i < plottypes_nr ; i++){
		printf("%s ",plottypes[i]);
	}
	printf("\n");
}


void usage(int argc, char* argv[]){
	printf("Usage: %s --plottype TYPE [--xmin XMIN --xmax XMAX --xlogscale]\n", argv[0]);
	plottype_usage();
}



int main( int argc, char* argv[]){

	int c;

	char *plottype = NULL;
	double x_start = -1;
	double x_stop = -1;
	int number_of_points_on_x_axis = -1;

	bool xlogscale_flag = false;

	while (1)
	{
		static struct option long_options[] = {
				{"plottype",  required_argument,    0, 'p'},
				{"xmin",   required_argument,    0, 'n'},
				{"xmax",   required_argument,    0, 'x'},
				{"xlogscale", no_argument  , 0, 'l'},
				{"npoints", required_argument  , 0, 'm'},
				{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "p:n:x:lm:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c) {
		case 0:
			/* If this option set a flag, do nothing else now. */
			if (long_options[option_index].flag != 0)
				break;
			printf ("option %s", long_options[option_index].name);
			if (optarg)
				printf (" with arg %s", optarg);
			printf ("\n");
			break;

		case 'p':
			plottype = optarg;
			break;

		case 'n':
		{
			char * test1 = 0;
			x_start = strtod(optarg,&test1);
			if( ((x_start == 0) && (*test1 != 0)) || (*test1 != '\0')){
				printf("Error in decoding xmin (--xmin option)\n");
			}
			break;
		};

		case 'x':
		{
			char * test2 = 0;
			x_stop = strtod(optarg,&test2);
			if( ((x_stop == 0) && (*test2 != 0)) || (*test2 != '\0')){
				printf("Error in decoding xmax (--xmax option)\n");
			}
			break;
		};

		case 'l':
			xlogscale_flag = true;
			break;

		case 'm':
		{
			char * test3 = 0;
			number_of_points_on_x_axis = strtol(optarg,&test3,10);
			if( ((number_of_points_on_x_axis == 0) && (*test3 != 0)) || (*test3 != '\0')){
				printf("Error in decoding npoints (--npoints option)\n");
			}
			break;
		};

		case '?': /* getopt_long already printed an error message. */
			break;

		default:
			abort ();
		}
	}

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}

	if( plottype == NULL){
		printf("Please specify plottype using -p option\n");
		usage(argc,argv);
		exit(EXIT_FAILURE);
	}

	/* prepare input */
	const long particle_no = 6012;  // carbon ion
	char particle_name[PARTICLE_NAME_NCHAR];
	AT_particle_name_from_particle_no_single(particle_no, particle_name);

	const long material_no = Water_Liquid;  // water

	if( strcmp( plottype , "ER") == 0 || strcmp( plottype , "LET") == 0){
		if(x_start < 0){
			x_start = 0.1;			printf("xmin not set, setting default value %g [MeV]\n", x_start);
		}
		if(x_stop < 0){
			x_stop = 500.0;			printf("xmax not set, setting default value %g [MeV]\n", x_stop);
		}
	} else if( strcmp( plottype , "RDD") == 0){
		if(x_start < 0){
			x_start = 1e-12;			printf("xmin not set, setting default value %g [m]\n", x_start);
		}
		if(x_stop < 0){
			x_stop = 1e-6;			printf("xmax not set, setting default value %g [m]\n", x_stop);
		}
	} else {
		if(x_start < 0){
			x_start = 1e-12;			printf("xmin not set, setting default value %g\n", x_start);
		}
		if(x_stop < 0){
			x_stop = 1e-6;			printf("xmax not set, setting default value %g\n", x_stop);
		}
	}

	if( number_of_points_on_x_axis < 0 ){
		number_of_points_on_x_axis = 10;			printf("npoints not set, setting default value %d\n", number_of_points_on_x_axis);
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
		long er_model = ER_Tabata;
		char er_name[200];
		getERName( er_model, er_name);
		printf("#ionenergy[MeV] range[m]\n");
		AT_max_electron_ranges_m( number_of_points_on_x_axis , x, material_no, er_model, y);
	} else if( strcmp( plottype , "RDD") == 0){
		long rdd_model = RDD_KatzExtTarget;
		char rdd_name[200];
		AT_RDD_name_from_number(rdd_model, rdd_name);
		printf("#RDD: %s; particle: %s (code: %ld)\n", rdd_name, particle_name, particle_no);
		printf("#r[m] D[Gy]\n");
		double rdd_parameter[RDD_MAX_NUMBER_OF_PARAMETERS] = {5e-11,1e-8,1e-10};
		double E_MeV_u = 100.0;
		long er_model = ER_Tabata;
		AT_D_RDD_Gy( number_of_points_on_x_axis, x, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, y);
	} else if( strcmp( plottype , "LET") == 0){
		printf("#LET\n");
		printf("#E[MeV] LET[MeV/cm2g]\n");
		for( i = 0 ; i < number_of_points_on_x_axis ; i++){
			y[i] = AT_LET_MeV_cm2_g_single(x[i],particle_no, material_no);
		}
	}  else {
		printf("Plottype %s not supported\n", plottype);
		plottype_usage();
	}

	/* save to file (stdout) */
	for( i = 0 ; i < number_of_points_on_x_axis ; i++){
		printf("%g %g\n",x[i],y[i]);
	}

	free(x);
	free(y);

	return EXIT_SUCCESS;
}
