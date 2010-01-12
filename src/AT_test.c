/**
 * @file
 * @brief ...
 */

/*
*    AT_test.c
*    ===================
*
*    Created on: 2009-06-08
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

#include "AmTrack.h"

void test_AT_SPISS(){

  // INPUT :
  long          n                       = 1;
  float 	E_MeV_u[] 		= {10.,50};
  long 		particle_no[] 		= {1, 1};
  float 	fluence_cm2[] 		= {-3.f, -0.005f};
  long		material_no 		= 1;			// Water

  long		RDD_model		= 3;			// Gei�
  float 	RDD_parameters[] 	= {5e-8};
  long		ER_model		= 3;			// Walig�rski
  float		ER_parameters[]		= {0.0f};
  long		GR_model		= 4;			// Exp-sat
  float		GR_parameters[]		= {1, 10};

  long		n_runs			= 1e5;
  long 		N2 			= 40;
  float		fluence_factor		= 1.0f;
  int		write_output		= 1;
  long		importance_sampling	= 0;

  float		results[10];

  AT_SPISS(	&n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      &material_no,
      &RDD_model,
      RDD_parameters,
      &ER_model,
      ER_parameters,
      &GR_model,
      GR_parameters,
      &n_runs,
      &N2,
      &fluence_factor,
      &write_output,
      &importance_sampling,
      results);

}

void test_AT_SPIFF(){

  // INPUT :
  long 		n 			= 1;
  float 	E_MeV_u[] 		= {10.,50};
  long 		particle_no[] 		= {1, 1};
  float 	fluence_cm2[] 		= {-3.f, -0.005f};
  long		material_no 		= 1;			// Water

  long		RDD_model		= 3;			// Gei�
  float 	RDD_parameters[] 	= {5e-8};
  long		ER_model		= 3;			// Walig�rski
  float		ER_parameters[]		= {0.0f};
  long		GR_model		= 4;			// Exp-sat
  float		GR_parameters[]		= {1, 10};

  long 		N2 			= 40;
  float		fluence_factor		= 1.0f;
  int		write_output		= 1;
  int		shrink_tails		= 1;
  float		shrink_tails_under	= 1e-30;
  int		adjust_N2		= 1;

  float		results[10];

  AT_SPIFF(	&n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      &material_no,
      &RDD_model,
      RDD_parameters,
      &ER_model,
      ER_parameters,
      &GR_model,
      GR_parameters,
      &N2,
      &fluence_factor,
      &write_output,
      &shrink_tails,
      &shrink_tails_under,
      &adjust_N2,
      results);

}

void test_AT_GSM(){

  // INPUT :
  long 		n 			= 1;
  float 	E_MeV_u[] 		= {10};
  long 		particle_no[] 		= {1, 1, 1};
  float 	fluence_cm2[] 		= {-3.f, 1e4, 1e4};
  long		material_no 		= 1;			// Water

  long		RDD_model		= 3;			// Gei�
  float 	RDD_parameters[] 	= {5e-8};
  long		ER_model		= 3;			// Walig�rski
  float		ER_parameters[]		= {0.0f};
  long		GR_model		= 4;			// Exp-sat
  float		GR_parameters[]		= {1, 10};

  long		n_runs			= 1;
  long 		N2 			= 40;
  float		fluence_factor		= 1.0f;
  bool		write_output		= true;

  long		nX			= 200;
  float		grid_size_m		= 5e-8;
  bool		lethal_events_mode	= false;

  float		results[10];

  AT_GSM(	&n,
      E_MeV_u,
      particle_no,
      fluence_cm2,
      &material_no,
      &RDD_model,
      RDD_parameters,
      &ER_model,
      ER_parameters,
      &GR_model,
      GR_parameters,
      &n_runs,
      &N2,
      &fluence_factor,
      &write_output,
      &nX,
      &grid_size_m,
      &lethal_events_mode,
      results);
}

int main(){
  //  test_AT_SPISS();
  test_AT_GSM();
  //	test_AT_SPIFF();

  return 0;
};
