/*
 * @file
 * @brief Dummy file to save dose pattern to output file.
 */

/*
 *    AT_test_GSM.c
 *    ===================
 *
 *    Created on: 2010-06-11
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
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "AmTrack.h"
#include "AT_ElectronRange.h"

void AT_IGK_save_survival_curve_to_file( const long  number_of_datapoints,
    const double   max_dose_Gy,
    const double   E_MeV_u,
    const long     particle_no,
    const long     material_no,
    const long     rdd_model,
    const double   rdd_parameter[],
    const long     er_model,
    const double   D0_characteristic_dose_Gy,
    const double   m_number_of_targets,
    const double   kappa,
    const char *   filename){

	FILE * outfile = fopen( filename, "w" );

	double dose_Gy;
	for( dose_Gy = 0 ; dose_Gy < max_dose_Gy ; dose_Gy += max_dose_Gy/number_of_datapoints ){
		double fluence_cm2;
		AT_fluence_cm2(  1,
		    &E_MeV_u,
		    &particle_no,
		    &dose_Gy,
		    material_no,
		    &fluence_cm2);
		double survival = AT_KatzModel_single_field_survival( fluence_cm2, E_MeV_u, particle_no, material_no, rdd_model, rdd_parameter, er_model, D0_characteristic_dose_Gy, m_number_of_targets, kappa );
		fprintf( outfile, "%g %g\n", dose_Gy, survival);
	}

	fclose(outfile);
}


int main( int argc, char * argv[] ){

	char filename1[200] = "sf.dat";

	if( argc != 8){
		printf("Usage: %s nPoints maxDose_Gy E_MeV D0_Gy m a0_m kappa\n", argv[0]);
		printf("\t\texample: %s 50 15.0 10.0 5.0 2.0 1e-8 2000\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	printf("Saving SF to file %s\n", filename1);

	long     number_of_points = atoi(argv[1]);
	double   max_dose_Gy = atof(argv[2]);
	double   E_MeV_u   = atof(argv[3]);

    double   D0_characteristic_dose_Gy = atof(argv[4]);
    double   m_number_of_targets = atof(argv[5]);
    double   a0_m = atof(argv[6]);
    double   kappa = atof(argv[7]);

	long     particle_no = 6012;
	long     material_no = Water_Liquid;
	long     rdd_model   = RDD_KatzExtTarget;
	double   r_min_m     = 1e-10;
	double   D_cutoff_Gy = 1e-40;
	double   rdd_parameter[] = {r_min_m,a0_m,D_cutoff_Gy};
	long     er_model    = ER_Waligorski;

	AT_IGK_save_survival_curve_to_file( number_of_points,
		max_dose_Gy,
	    E_MeV_u,
	    particle_no,
	    material_no,
	    rdd_model,
	    rdd_parameter,
	    er_model,
	    D0_characteristic_dose_Gy,
	    m_number_of_targets,
	    kappa,
	    filename1);

	return EXIT_SUCCESS;
};
