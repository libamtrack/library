/**
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

int main(){

  long*  particle_no = NULL;
  long  cur_particle_no;
  float*  E_MeV_u = NULL;
  float  cur_E_MeV_u;
  float*  fluence_cm2 = NULL;
  float  cur_fluence_cm2;
  long  n_particles = 0;

  printf("############################################################\n");
  printf("This is AmTrack & AmTrack user interface, version 2009/09/20\n");
  printf("Copyright 2009, libamtrack project, licensed under GNU v3\n");
  printf("Please see documentation (libamtrack.sourceforge.net) for\nparticle and material indices etc.\n");
  printf("############################################################\n");

  do{
    printf("\n> %d. particle index (0 for end): ", n_particles + 1);
    scanf("%d", &cur_particle_no);

    if(cur_particle_no > 0){
      printf("> %d. particle energy (MeV/u): ", n_particles + 1);
      scanf("%f", &cur_E_MeV_u);

      printf("> %d. particle fluence (1/cm2) or dose (Gy; enter negative number): ", n_particles + 1);
      scanf("%f", &cur_fluence_cm2);

      n_particles += 1;

      particle_no          =  (long*)realloc(particle_no,   n_particles * sizeof(long));
      E_MeV_u            =  (float*)realloc(E_MeV_u,     n_particles * sizeof(float));
      fluence_cm2          =  (float*)realloc(fluence_cm2,  n_particles * sizeof(float));
      particle_no[n_particles-1]  =  cur_particle_no;
      E_MeV_u[n_particles-1]    =  cur_E_MeV_u;
      fluence_cm2[n_particles-1]  =  cur_fluence_cm2;
    }

  }while(cur_particle_no > 0);

  if(n_particles == 0){
    printf("\n####\nbye.\n####\n");
    return EXIT_FAILURE;
  }

  long   material_no;
  long  RDD_model;
  float*  RDD_parameters;
  long*  ER_model;
  float*  ER_parameters;
  long  gamma_model;
  float*  gamma_parameters;

  char  output_dummy[100];

  printf("\n> material index: ");
  scanf("%d", &material_no);
  getMaterialName(material_no, output_dummy);
  printf("%s selected.\n", output_dummy);

  printf("\n> radial-dose model: ");
  scanf("%d", &RDD_model);


  printf("\n> electron-range model: ");
  scanf("%d", &ER_model);

  printf("\n> gamma-reponse model: ");
  scanf("%d", &gamma_model);



  long i;
  for (i = 0; i < n_particles; i++){
    printf("\nO %d particle index: %d, Energy / (MeV/u): %g",   i,
                                  particle_no[i],
                                  E_MeV_u[i]);
  }


  return EXIT_SUCCESS;
};
