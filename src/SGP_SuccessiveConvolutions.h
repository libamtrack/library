#ifndef SGP_SUCCESSIVECONVOLUTIONS_H_
#define SGP_SUCCESSIVECONVOLUTIONS_H_

#include <stdbool.h>

void	SGP_SC_get_f1_array_size(	/* radiation field parameters */
									long*	n,
									float*	E_MeV_u,
									long*	particle_no,
									/* detector parameters */
									long*	material_no,
									long*	rdd_model,
									float*	rdd_parameter,
									/* electron range model */
									long*	er_model,
									float*	er_parameter,
									/* algorith parameters*/
									long*	N2,
									// from here: return values
									long*	n_bins_f1,
									float*	f1_parameters);


void	SGP_SC_get_f1(	/* radiation field parameters */
						long*	n,
						float*	E_MeV_u,
						long*	particle_no,
						float*	fluence_cm2,
						/* detector parameters */
						long*	material_no,
						long*	rdd_model,
						float*	rdd_parameter,
						/* electron range model */
						long*	er_model,
						float*	er_parameter,
						/* algorith parameters*/
						long*	N2,
						long*	n_bins_f1,
						/* f1 parameters*/
						float*	f1_parameters,
						// from here: return values
						float*	norm_fluence,
						float*	dose_contribution_Gy,
						float*	f_parameters,
						/*  1 - total fluence_cm2
						 *  2 - total_dose_Gy
						 *  3 - ave_E_MeV
						 *  4 - dw_E_MeV
						 *  5 - ave_LET_MeV_cm2_g
						 *  6 - dw_LET_MeV_cm2_g
						 *  0 - u
						 */
						float*	f1_d_Gy,
						float*	f1_dd_Gy,
						float*	f1);

void SGP_SC_get_f_array_size(	float*	u,
		float*	fluence_factor,
		long*	N2,
		long*	n_bins_f1,
		float*	f1_d_Gy,
		float*	f1_dd_Gy,
		float*	f1,
		// from here: return values
		long*	n_bins_f,
		float*	u_start,
		long*	n_convolutions);


void	SGP_SC_get_f_start(	float*	u_start,
		long*	n_bins_f1,
		long*	N2,
		float*	f1_d_Gy,
		float*	f1_dd_Gy,
		float*	f1,
		long*	n_bins_f,
		// from here: return values
		float*	f_d_Gy,
		float*	f_dd_Gy,
		float*	f_start);


void SGP_SuccessiveConvolutions(	float*	u,
		long*	n_bins_f,
		long*	N2,
		// input + return values
		long*	n_bins_f_used,
		float*	f_d_Gy,
		float*	f_dd_Gy,
		float*	f,
		// return values
		float*	f0,
		float*	fdd,
		float*	dfdd,
		float*	d,
		bool*	write_output,
		bool*	shrink_tails,
		float*	shrink_tails_under,
		bool*	adjust_N2);



#endif // SGP_SUCCESSIVECONVOLUTIONS_H_
