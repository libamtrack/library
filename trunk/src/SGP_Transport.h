#ifndef SGP_TRANSPORT_H_
#define SGP_TRANSPORT_H_


void SGP_BortfeldTransportProton(	float*	E_initial_MeV,
									float*	sE_initial_MeV,
									float*	fluence_initial_cm2,
									char*	plateau_dose_material_name,
									long*	n_shielding_slabs,
									float*	shielding_thickness_m,
									char**	shielding_material_name,
									long*	n_detector_slabs,
									float*	detector_thickness_m,
									char*	detector_material_name,
									/* return values */
									float*  plateau_dose_Gy,
									float*  plateau_dose_noNuc_Gy,
									float*	detector_z_cm,
									float*	detector_E_MeV,
									float*	detector_sE_MeV,
									float*	detector_fluence_cm2,
									float*	detector_dfluencedz_cm,
									float*	detector_LET_MeV_g_cm2,
									float*	detector_dose_Gy,
									float*	detector_dose_noNuc_Gy,
									float*	geom_factor);


void SGP_BortfeldTransportProtonL(	float*	E_initial_MeV,
									float*	sE_initial_MeV,
									float*	fluence_initial_cm2,
									char*	plateau_dose_material_name,
									long*	n_shielding_slabs,
									float*	shielding_thickness_m,
									char*	shielding_material_nameL,
									long*	n_detector_slabs,
									float*	detector_thickness_m,
									char*	detector_material_name,
									/* return values */
									float*  plateau_dose_Gy,
									float*  plateau_dose_noNuc_Gy,
									float*	detector_z_cm,
									float*	detector_E_MeV,
									float*	detector_sE_MeV,
									float*	detector_fluence_cm2,
									float*	detector_dfluencedz_cm,
									float*	detector_LET_MeV_g_cm2,
									float*	detector_dose_Gy,
									float*	detector_dose_noNuc_Gy,
									float*	geom_factor);


void SGP_BortfeldTransportProtonS(	float*	E_initial_MeV,
									float*	sE_initial_MeV,
									float*	fluence_initial_cm2,
									char**	plateau_dose_material_name,
									long*	n_shielding_slabs,
									float*	shielding_thickness_m,
									char**	shielding_material_name,
									long*	n_detector_slabs,
									float*	detector_thickness_m,
									char**	detector_material_name,
									/* return values */
									float*  plateau_dose_Gy,
									float*  plateau_dose_noNuc_Gy,
									float*	detector_z_cm,
									float*	detector_E_MeV,
									float*	detector_sE_MeV,
									float*	detector_fluence_cm2,
									float*	detector_dfluencedz_cm,
									float*	detector_LET_MeV_g_cm2,
									float*	detector_dose_Gy,
									float*	detector_dose_noNuc_Gy,
									float*	geom_factor);

void	SGP_Dose_Gy(	float*	density_g_cm3,
						float*	fluence_cm2,
						float*	LET_MeV_g_cm2,
						float*	dfluencedz_cm,
						float*	E_MeV,
						float*	dose_Gy,
						float*	dose_noNuc_Gy);
#endif // SGP_TRANSPORT_H_
