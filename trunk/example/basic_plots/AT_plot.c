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
	printf("Usage: %s -p plottype\n", argv[0]);
	plottype_usage();
}



int main( int argc, char* argv[]){

	int c;

	char *plottype = NULL;

	while (1)
	{
		static struct option long_options[] = {
				{"plottype",   required_argument,    0, 'p'},
				{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "p:", long_options, &option_index);

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
	const long material_no = Water_Liquid;  // water

	int number_of_points_on_x_axis = 10;
	const double x_start = 1e-1;
	const double x_stop = 500.0;
	const double x_multiplier = exp(((log( x_stop ) - log( x_start)) / ((double)number_of_points_on_x_axis-1.)));
	double * x = (double*)calloc(number_of_points_on_x_axis,sizeof(double));
	double * y = (double*)calloc(number_of_points_on_x_axis,sizeof(double));
	int i;
	x[0] = x_start;
	for( i = 1 ; i < number_of_points_on_x_axis ; i++){
		x[i] = x[i-1] * x_multiplier;
	}

	/* calculate results */
	if( strcmp( plottype , "ER") == 0){
		printf("#ER: Tabata\n");
		printf("#ionenergy[MeV] range[m]\n");
		long er_model = ER_Tabata;
		AT_max_electron_ranges_m( number_of_points_on_x_axis , x, material_no, er_model, y);
	} else if( strcmp( plottype , "RDD") == 0){
		printf("#RDD: Katz\n");
		printf("#r[m] D[Gy]\n");
		long rdd_model = RDD_KatzExtTarget;
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
