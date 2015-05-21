/**
 * @brief Radial Dose Distribution models
 */

/*
 *    AT_RDD_Simple.c
 *    ========
 *
 *    Created on: 05.04.2010
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

#include "AT_RDD_Tabulated.h"

/* --------------------------------------------------- Radical Diffusion RDD ---------------------------------------------------*/

 long    AT_RDD_RadicalDiffusion_get_energy_idx(const double E_MeV_u){
	 return(locate( AT_RadDiff_RDD.E_MeV_u,
	 			 N_ENERGY,
	 			 E_MeV_u));
 }

 double  AT_RDD_RadicalDiffusion_Gy( const double r_m,
		 const double E_MeV_u,
		 const long particle_no,
		 const long material_no){

	 long    E_idx = AT_RDD_RadicalDiffusion_get_energy_idx(E_MeV_u);

	 assert(E_idx >= 1);
	 assert(E_idx <= N_ENERGY);
	 assert(r_m >= 0.0);

	 if(r_m <= AT_r_min_RadicalDiffusion_m(E_idx)){
		 return (AT_d_max_RadicalDiffusion_Gy( E_MeV_u,
				 particle_no,
				 material_no));
	 }

	 if(r_m > AT_r_max_RadicalDiffusion_m(E_MeV_u)){
		 return (0.0);
	 }

	 double LET_keV_um;
	 AT_Stopping_Power_with_no( PSTAR,
				1,
				&E_MeV_u,
				&particle_no,
				material_no,
				&LET_keV_um);

	 return( 0.160219 *
			 LET_keV_um *
			 AT_get_interpolated_y_from_input_table( AT_RadDiff_RDD.r_m[E_idx-1],
			 AT_RadDiff_RDD.d_Gy[E_idx-1],
			 AT_n_bins_RadicalDiffusion(E_MeV_u),
	 	 	 r_m));

}


 double  AT_inverse_RadicalDiffusion_m( const double d_Gy,
		 const double E_MeV_u,
		 const long particle_no,
		 const long material_no){

	 long    E_idx = AT_RDD_RadicalDiffusion_get_energy_idx(E_MeV_u);

	 assert(E_idx >= 1);
	 assert(E_idx <= 40);
	 assert(d_Gy >= 0.0);

	 if((d_Gy < AT_d_min_RadicalDiffusion_Gy(E_MeV_u, particle_no, material_no)) ||
		(d_Gy > AT_d_max_RadicalDiffusion_Gy(E_MeV_u, particle_no, material_no))){
		 return (0.0);
	 }

	 double LET_keV_um;
	 AT_Stopping_Power_with_no( PSTAR,
	 				1,
	 				&E_MeV_u,
	 				&particle_no,
	 				material_no,
	 				&LET_keV_um);

	 if(LET_keV_um == 0){
		 return 0.0;
	 }


	 return( AT_get_interpolated_y_from_input_table( AT_RadDiff_RDD.d_Gy[E_idx-1],
			 AT_RadDiff_RDD.r_m[E_idx-1],
			 AT_n_bins_RadicalDiffusion(E_MeV_u),
	 	 	 d_Gy / (0.160219 * LET_keV_um)));

 }

 double	 AT_E_RadicalDiffusion_MeV_u( const double E_MeV_u){
	 long    E_idx = AT_RDD_RadicalDiffusion_get_energy_idx(E_MeV_u);

	 assert(E_idx >= 1);
	 assert(E_idx <= 40);

     return AT_RadDiff_RDD.E_MeV_u[E_idx - 1];
 }

 double	 AT_r_min_RadicalDiffusion_m( const double E_MeV_u){
	 long    E_idx = AT_RDD_RadicalDiffusion_get_energy_idx(E_MeV_u);

	 assert(E_idx >= 1);
	 assert(E_idx <= 40);

     return AT_RadDiff_RDD.r_min_m[E_idx - 1];
 }

 double	 AT_r_max_RadicalDiffusion_m( const double E_MeV_u){
	 long    E_idx = AT_RDD_RadicalDiffusion_get_energy_idx(E_MeV_u);

	 assert(E_idx >= 1);
	 assert(E_idx <= 40);

	 return AT_RadDiff_RDD.r_max_m[E_idx - 1];
 }

 double	 AT_d_min_RadicalDiffusion_Gy( const double E_MeV_u,
		 const long particle_no,
		 const long material_no){
	 long    E_idx = AT_RDD_RadicalDiffusion_get_energy_idx(E_MeV_u);

	 assert(E_idx >= 1);
	 assert(E_idx <= 40);

	 double LET_keV_um;
	 AT_Stopping_Power_with_no( PSTAR,
	 				1,
	 				&E_MeV_u,
	 				&particle_no,
	 				material_no,
	 				&LET_keV_um);

	 return 0.160219 * LET_keV_um * AT_RadDiff_RDD.d_min_Gy[E_idx - 1];
 }

 double	 AT_d_max_RadicalDiffusion_Gy( const double E_MeV_u,
		 const long particle_no,
		 const long material_no){
	 long    E_idx = AT_RDD_RadicalDiffusion_get_energy_idx(E_MeV_u);

	 assert(E_idx >= 1);
	 assert(E_idx <= 40);

	 double LET_keV_um;
	 AT_Stopping_Power_with_no( PSTAR,
	 	 				1,
	 	 				&E_MeV_u,
	 	 				&particle_no,
	 	 				material_no,
	 	 				&LET_keV_um);

	 return 0.160219 * LET_keV_um * AT_RadDiff_RDD.d_max_Gy[E_idx - 1];
 }

 long	 AT_n_bins_RadicalDiffusion( const double E_MeV_u){
	 long    E_idx = AT_RDD_RadicalDiffusion_get_energy_idx(E_MeV_u);

	 assert(E_idx >= 1);
	 assert(E_idx <= 40);

	 return (AT_RadDiff_RDD.n_bins[E_idx - 1]);
 }


