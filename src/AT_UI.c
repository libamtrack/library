/**
 * @file
 * @brief ...
 */

/*
*    AT_UI.c
*    ==============
*
*    A simple user interface to libamtrack
*    Created on: 20.09.2009
*    Author: greilich
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

#include "AmTrack.h"
#include "AT_ElectronRange.h"
#include "AT_RDD.h"
#include "AT_GammaResponse.h"

int main(){

  long    i;
  long*   particle_no     = NULL;
  long    cur_particle_no = 1;
  float*  E_MeV_u         = NULL;
  float   cur_E_MeV_u;
  float*  fluence_cm2     = NULL;
  float   cur_fluence_cm2;
  long    n_particles     = 0;
  int     retval;

  printf("############################################################\n");
  printf("This is AmTrack & AmTrack user interface, version 2009/09/24\n");
  printf("Copyright 2009, libamtrack project, licensed under GNU v3\n");
  printf("Please see documentation (libamtrack.sourceforge.net) for\nparticle and material indices etc.\n");
  printf("############################################################\n");

  do{
    printf("\n> %ld. particle index (0 for end, 999 for list): ", n_particles + 1);
    retval = scanf("%ld", &cur_particle_no);
    if( retval != 1 ){
      printf("Error ! Provide correct number\n");
      exit(EXIT_FAILURE);
    }

    if(cur_particle_no > 0){
      if(cur_particle_no != 999){
        printf("> %ld. particle energy (MeV/u): ", n_particles + 1);

        retval = scanf("%f", &cur_E_MeV_u);
        if( retval != 1 ){
          printf("Error ! Provide correct number\n");
          exit(EXIT_FAILURE);
        }

        printf("> %ld. particle fluence (1/cm2) or dose (Gy; enter negative number): ", n_particles + 1);
        retval = scanf("%f", &cur_fluence_cm2);
        if( retval != 1 ){
          printf("Error ! Provide correct number\n");
          exit(EXIT_FAILURE);
        }
        n_particles += 1;

        particle_no          =   (long*)realloc(particle_no,  n_particles * sizeof(long));
        E_MeV_u              =  (float*)realloc(E_MeV_u,      n_particles * sizeof(float));
        fluence_cm2          =  (float*)realloc(fluence_cm2,  n_particles * sizeof(float));
        particle_no[n_particles-1]  =  cur_particle_no;
        E_MeV_u[n_particles-1]      =  cur_E_MeV_u;
        fluence_cm2[n_particles-1]  =  cur_fluence_cm2;
      }else{
        for (i = 0; i < PARTICLE_DATA_N; i++){
          printf("\nparticle index: %2ld --> %s", AT_Particle_Data.particle_no[i],AT_Particle_Data.particle_name[i]);
        }
      }
    }

  }while(cur_particle_no > 0);

  if(n_particles == 0){
    printf("\n####\nbye.\n####\n");
    return EXIT_FAILURE;
  }

  long   material_no;
  long   RDD_model;
  float  RDD_parameters[10];
  long   ER_model;
  float  ER_parameters[10];
  long   gamma_model;
  float  gamma_parameters[10];

  char   output_dummy[100];
  float  float_dummy;

  do{
    printf("\n> ** Select material **");
    long  i;
    for(i = 0; i < MATERIAL_DATA_N;i++){
      printf("\n index: %ld --> %s", AT_Material_Data.material_no[i], AT_Material_Data.material_name[i]);
    }
    printf("\n> material index: ");
    retval = scanf("%ld", &material_no);
    if( retval != 1 ){
      printf("Error ! Provide correct number\n");
      exit(EXIT_FAILURE);
    }
    getMaterialName(material_no, output_dummy);
    printf("%s selected.\n", output_dummy);
  }while(strcmp(output_dummy, "*** invalid choice ***") == 0);

  do{
    printf("\n> ** Select radial dose distribution **");
    long  i;
    for(i = 0; i < RDD_DATA_N;i++){
      printf("\n index: %ld --> %s", AT_RDD_Data.RDD_no[i], AT_RDD_Data.RDD_name[i]);
    }
    printf("\n> radial-dose model index: ");
    retval = scanf("%ld", &RDD_model);
    if( retval != 1 ){
      printf("Error ! Provide correct number\n");
      exit(EXIT_FAILURE);
    }
    getRDDName(&RDD_model, output_dummy);
    printf("%s selected.\n", output_dummy);
  }while(strcmp(output_dummy, "*** invalid choice ***") == 0);

  //TODO not fully correct, works only if RDD_index = RDD_no !!

  long index = RDD_model-1;
  if(AT_RDD_Data.n_parameters[index]>0){
    printf("\n> ** Select radial dose distribution parameters **");
    for(i = 0; i < AT_RDD_Data.n_parameters[index];i++){
      printf("\n %s [0 for default: %g]: ", AT_RDD_Data.parameter_name[index][i],
          AT_RDD_Data.parameter_default[index][i]);
      retval = scanf("%g", &float_dummy);
      if( retval != 1 ){
        printf("Error ! Provide correct number\n");
        exit(EXIT_FAILURE);
      }
      if (float_dummy == 0.0f){float_dummy = AT_RDD_Data.parameter_default[index][i];}
      printf("%g understood.\n", float_dummy);
      RDD_parameters[i] = float_dummy;
    }
  }

  do{
    printf("\n> *** select electron-range model ***");
    for(i = 0; i < ER_DATA_N;i++){
      printf("\n index: %d --> %s", AT_ER_Data.ER_no[i], AT_ER_Data.ER_name[i]);
    }
    printf("\n> electron-range model index: ");
    retval = scanf("%ld", &ER_model);
    if( retval != 1 ){
      printf("Error ! Provide correct number\n");
      exit(EXIT_FAILURE);
    }
    getERName(ER_model, output_dummy);
    printf("%s selected.\n", output_dummy);
  }while(strcmp(output_dummy, "*** invalid choice ***") == 0);

  do{
    printf("\n> *** select gamma-response model ***");
    for(i = 0; i < GR_DATA_N;i++){
      printf("\n index: %ld --> %s", AT_GR_Data.GR_no[i], AT_GR_Data.GR_name[i]);
    }
    printf("\n> gamma response model index: ");
    retval = scanf("%ld", &gamma_model);
    if( retval != 1 ){
      printf("Error ! Provide correct number\n");
      exit(EXIT_FAILURE);
    }
    getGammaName(&gamma_model, output_dummy);
    printf("%s selected.\n", output_dummy);
  }while(strcmp(output_dummy, "*** invalid choice ***") == 0);

  if(AT_GR_Data.n_parameters[index]>0){
    printf("\n> ** Select gamma response parameters **");
    index = gamma_model-1; //TODO
    for(i = 0; i < AT_GR_Data.n_parameters[index];i++){
      printf("\n %s [0 for default: %g]: ", AT_GR_Data.parameter_name[index][i], AT_GR_Data.parameter_default[index][i]);
      retval = scanf("%g", &float_dummy);
      if( retval != 1 ){
        printf("Error ! Provide correct number\n");
        exit(EXIT_FAILURE);
      }
      if (float_dummy == 0.0f){float_dummy = AT_GR_Data.parameter_default[index][i];}
      printf("%g understood.\n", float_dummy);
      gamma_parameters[i] = float_dummy;
    }
  }

  printf("\n\n >>> Computing efficiency now, using SPIFF algorithm with default settings.");
  float results[10];
  long  N2 = 10;
  float fluence_factor     = 1.0f;
  bool   write_output      = true;
  bool   shrink_tails      = true;
  bool   adjust_N2         = true;
  float shrink_tails_under = 1e-30f;
  bool  lethal_events_mode = false;

  AT_SPIFF(  &n_particles,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      &material_no,
      &RDD_model,
      RDD_parameters,
      &ER_model,
      ER_parameters,
      &gamma_model,
      gamma_parameters,
      &N2,
      &fluence_factor,
      &write_output,
      &shrink_tails,
      &shrink_tails_under,
      &adjust_N2,
      &lethal_events_mode,
      results);

  printf("\n\n >>> Done.");

  printf("\nEfficiency:          %g", results[0]);
  printf("\nDose check / Gy:     %g", results[1]);
  printf("\nParticle response:   %g", results[2]);
  printf("\nGamma response:       %g", results[3]);
  printf("\n");
  printf("\nMean impact number:  %g", results[5]);
  printf("\nStart impact number: %g", results[6]);
  printf("\nNo. of convolutions: %d", (int)(results[7]));

  printf("\n\n####\nbye.\n####\n");
  return EXIT_SUCCESS;
};
