#!/bin/bash

if [ "$1" != "" ]; then
    echo "WASM parameter set to $1"
    WASM=$1
else
    echo "WASM parameter EMPTY !!! Default set to 1"
    WASM=1
fi

#emcmake="/home/osboxes/emsdk/emscripten/1.37.36/./emcmake"
#emmake="/home/osboxes/emsdk/emscripten/1.37.36/./emmake"
#emcc="/home/osboxes/emsdk/emscripten/1.37.36/./emcc"

cd .
mkdir _build
cd _build
cp ../libgsl.a .
ls -al .

emcmake cmake .. -DGSL_INCLUDE_DIRS=$GSL_INCLUDE_DIRS -DGSL_LIBRARY=$GSL_LIBRARY -DGSL_CBLAS_LIBRARY=$GSL_CBLAS_LIBRARY
emmake make -j4

funs='['	
 
  #----AT_Error.h	
 funs+='"_AT_get_error_msg",'	
 funs+='"_AT_check_energy_range_single_particle",'	
 funs+='"_AT_check_energy_range_single_field",'	
 funs+='"_AT_check_particle_no_single_particle",'	
 funs+='"_AT_check_particle_no_single_field",'	
 
  #----AT_StoppingPowerDataBethe.h	
 funs+='"_AT_Bethe_wrapper",'	
 funs+='"_AT_el_energy_loss_leading_term_MeV_cm2_g",'	
 funs+='"_AT_Bethe_Stopping_Number",'	
 funs+='"_AT_Bethe_energy_loss_MeV_cm2_g_single",'	
 funs+='"_AT_Bethe_energy_loss_MeV_cm2_g",'	

  #----AT_Histograms.h
 funs+='"_AT_histo_linear_left_limit",'	
 funs+='"_AT_histo_logarithmic_left_limit",'	
 funs+='"_AT_histo_left_limit",'	
 funs+='"_AT_histo_left_limits",'	
 funs+='"_AT_histo_linear_bin_width",'	
 funs+='"_AT_histo_logarithmic_bin_width",'	
 funs+='"_AT_histo_bin_width",'	
 funs+='"_AT_histo_bin_widths",'	
 funs+='"_AT_histo_linear_midpoint",'	
 funs+='"_AT_histo_logarithmic_midpoint",'	
 funs+='"_AT_histo_midpoint",'	
 funs+='"_AT_histo_midpoints",'	
 funs+='"_AT_histo_linear_step",'	
 funs+='"_AT_histo_logarithmic_step",'	
 funs+='"_AT_histo_step",'	
 funs+='"_AT_histo_linear_n_bins",'	
 funs+='"_AT_histo_logarithmic_n_bins",'	
 funs+='"_AT_histo_n_bins",'	
 funs+='"_AT_histo_linear_bin_no",'	
 funs+='"_AT_histo_logarithmic_bin_no",'	
 funs+='"_AT_histo_bin_no",'	
 funs+='"_AT_histo_add_single",'	
 funs+='"_AT_histo_add_multi",'	
 funs+='"_AT_histo_sum",'	
 funs+='"_AT_histo_normalize",'	
 funs+='"_AT_N2_to_step",'	
 funs+='"_AT_step_to_N2",'	
 funs+='"_AT_histoOld_log_bin_width",'	
 funs+='"_AT_histoOld_lower_bin_limit",'	
 funs+='"_AT_histoOld_upper_bin_limit",'	
 funs+='"_AT_histoOld_get_bin_width",'	
 funs+='"_AT_histoOld_get_bin_widths",'	
 funs+='"_AT_histoOld_bin_no",'	
 
  #----AT_MultipleCoulombScattering.h	
 funs+='"_AT_characteristic_single_scattering_angle_single",'	
 funs+='"_AT_characteristic_single_scattering_angle",'	
 funs+='"_AT_screening_angle_single",'	
 funs+='"_AT_screening_angle",'	
 funs+='"_AT_effective_collision_number_single",'	
 funs+='"_AT_effective_collision_number",'	
 funs+='"_AT_reduced_target_thickness_single",'	
 funs+='"_AT_reduced_target_thickness",'	
 funs+='"_AT_characteristic_multiple_scattering_angle_single",'	
 funs+='"_AT_characteristic_multiple_scattering_angle",'	
 funs+='"_AT_Moliere_function_f0",'	
 funs+='"_AT_Moliere_function_f1",'	
 funs+='"_AT_Moliere_function_f2",'	
 funs+='"_AT_scattering_angle_distribution_single",'	
 funs+='"_AT_scattering_angle_distribution",'	
 funs+='"_AT_Highland_angle_single",'	
 funs+='"_AT_Highland_angle",'	
 
  #----AT_RDD.h	
 funs+='"_AT_RDD_index_from_RDD_number",'	
 funs+='"_AT_RDD_name_from_number",'	
 funs+='"_AT_RDD_number_from_name",'	
 funs+='"_AT_RDD_number_of_parameters",'	
 funs+='"_AT_D_RDD_Gy",'	
 funs+='"_AT_r_RDD_m",'	
 funs+='"_AT_RDD_r_min_m",'	
 funs+='"_AT_RDD_a0_m",'	
 funs+='"_AT_RDD_precalculated_constant_Gy",'	
 funs+='"_AT_RDD_d_min_Gy",'	
 funs+='"_AT_RDD_d_max_Gy",'	
 funs+='"_AT_RDD_f1_parameters_single_field",'	
 funs+='"_AT_RDD_f1_parameters_mixed_field",'	
 
  #----AT_StoppingPowerDataPSTAR.h	
 funs+='"_AT_PSTAR_wrapper",'	
 
  #----AT_EnergyLoss.h	
 funs+='"_AT_mean_energy_loss_keV",'	
 funs+='"_AT_xi_keV",'	
 funs+='"_AT_kappa_multi",'	
 funs+='"_AT_kappa_single",'	
 funs+='"_AT_Landau_PDF",'	
 funs+='"_AT_Landau_IDF",'	
 funs+='"_AT_lambda_landau_from_energy_loss_multi",'	
 funs+='"_AT_lambda_landau_from_energy_loss_single",'	
 funs+='"_AT_lambda_mean_multi",'	
 funs+='"_AT_lambda_mean_single",'	
 funs+='"_AT_lambda_max_multi",'	
 funs+='"_AT_lambda_max_single",'	
 funs+='"_AT_lambda_Landau_Mode",'	
 funs+='"_AT_lambda_Landau_Mean",'	
 funs+='"_AT_lambda_Landau_FWHM_left",'	
 funs+='"_AT_lambda_Landau_FWHM_right",'	
 funs+='"_AT_lambda_Landau_FWHM",'	
 funs+='"_AT_energy_loss_keV_Landau_FWHM",'	
 funs+='"_AT_energy_loss_keV_Landau_Mode",'	
 funs+='"_AT_energy_loss_from_lambda_landau_multi",'	
 funs+='"_AT_energy_loss_from_lambda_landau_single",'	
 funs+='"_AT_Landau_energy_loss_distribution",'	
 funs+='"_AT_Vavilov_PDF",'	
 funs+='"_AT_Vavilov_IDF",'	
 funs+='"_AT_lambda_vavilov_from_energy_loss_multi",'	
 funs+='"_AT_lambda_vavilov_from_energy_loss_single",'	
 funs+='"_AT_lambda_Vavilov_Mode",'	
 funs+='"_AT_lambda_Vavilov_Mean",'	
 funs+='"_AT_lambda_Vavilov_Variance",'	
 funs+='"_AT_lambda_Vavilov_Skewness",'	
 funs+='"_AT_lambda_Vavilov_FWHM_left",'	
 funs+='"_AT_lambda_Vavilov_FWHM_right",'	
 funs+='"_AT_lambda_Vavilov_FWHM",'	
 funs+='"_AT_energy_loss_keV_Vavilov_FWHM",'	
 funs+='"_AT_energy_loss_from_lambda_vavilov_multi",'	
 funs+='"_AT_Vavilov_energy_loss_distribution",'
 funs+='"_AT_Gauss_energy_loss_distribution",'
 funs+='"_AT_Gauss_PDF",'	
 funs+='"_AT_Gauss_IDF",'	
 funs+='"_AT_energy_loss_from_lambda_gauss_multi",'	
 funs+='"_AT_Gauss_Mode",'	
 funs+='"_AT_Gauss_Mean",'	
 funs+='"_AT_Gauss_FWHM",'	
 funs+='"_AT_energy_loss_distribution",'	
 funs+='"_AT_energy_loss_mode",'	
 funs+='"_AT_energy_loss_FWHM",'	
 
  #----AT_NumericalRoutines.h	
 funs+='"_AT_range_straggling_convolution",'	
 funs+='"_AT_Dyx",'	
 funs+='"_AT_gamma_",'	
 funs+='"_AT_sum",'	
 funs+='"_AT_normalize",'	
 funs+='"_AT_get_interpolated_y_from_input_table",'	
 funs+='"_AT_get_interpolated_y_from_input_2d_table",'	
 funs+='"_AT_get_interpolated_x_from_input_2d_table",'	
 funs+='"_AT_get_interpolated_y_from_interval",'	
 
  #----AT_GammaResponse.h	
 funs+='"_AT_Gamma_index_from_material_number",'	
 funs+='"_AT_Gamma_name_from_number",'	
 funs+='"_AT_Gamma_number_of_parameters",'	
 funs+='"_AT_gamma_response",'	
 funs+='"_AT_get_gamma_response_for_average_dose",'	
 funs+='"_AT_get_response_distribution_from_dose_distribution",'	
 funs+='"_AT_get_ion_response_from_response_distribution",'	
 funs+='"_AT_get_ion_response_from_dose_distribution",'	
 funs+='"_AT_get_ion_efficiency_from_dose_distribution",'	
 funs+='"_AT_get_ion_efficiency_from_response_distribution",'	
 funs+='"_AT_get_gamma_response",'	
 
  #----AT_Constants.h	
 
  #----AT.h	
 funs+='"_AT_test_fun",'	
 
  #----AT_StoppingPowerDataFromFile.h	
 funs+='"_AT_FromFile_wrapper",'	
 
  #----AT_StoppingPower.h	
 funs+='"_AT_Mass_Stopping_Power",'	
 funs+='"_AT_Stopping_Power",'	
 funs+='"_AT_Mass_Stopping_Power_with_no",'	
 funs+='"_AT_Stopping_Power_with_no",'	
 funs+='"_AT_Energy_MeV_u_from_Stopping_Power_single",'	

 #----AT_ProtonAnalyticalModels.h
 funs+='"_AT_dose_Bortfeld_Gy_single",'
 funs+='"_AT_dose_Bortfeld_Gy_multi",'
 funs+='"_AT_LET_t_Wilkens_keV_um_single",'
 funs+='"_AT_LET_t_Wilkens_keV_um_multi",'
 funs+='"_AT_LET_d_Wilkens_keV_um_single",'
 funs+='"_AT_LET_d_Wilkens_keV_um_multi",'
 funs+='"_AT_proton_RBE_single",'
 funs+='"_AT_proton_RBE_multi",'

 #----AT_ProtonAnalyticalBeamParameters.h
 funs+='"_AT_range_Bortfeld_cm",'
 funs+='"_AT_fwhm_Bortfeld_cm",'
 funs+='"_AT_max_plateau_Bortfeld",'
 funs+='"_AT_energy_Bortfeld_MeV_u",'


  #----AT_SuccessiveConvolutions.h	
 funs+='"_AT_n_bins_for_single_impact_local_dose_distrib",'	
 funs+='"_AT_single_impact_local_dose_distrib",'	
 funs+='"_AT_n_bins_for_low_fluence_local_dose_distribution",'	
 funs+='"_AT_low_fluence_local_dose_distribution",'	
 funs+='"_AT_SuccessiveConvolutions",'	
 funs+='"_AT_Kellerer_normalize",'	
 funs+='"_AT_Kellerer_interpolation",'	
 funs+='"_AT_Kellerer_reset",'	
 funs+='"_AT_Kellerer_zero",'	
 funs+='"_AT_Kellerer_shrink",'	
 funs+='"_AT_Kellerer_folding",'	
 funs+='"_AT_n_bins_for_DSB_distribution",'	
 funs+='"_AT_get_DSB_distribution",'	
 funs+='"_AT_translate_dose_into_DSB_distribution",'	
 
  #----AT_Algorithms_GSM.h	
 funs+='"_AT_GSM_sample_particle_positions",'	
 funs+='"_AT_GSM_dose_grid_from_particles_positions",'	
 funs+='"_AT_GSM_local_dose_distrib_from_dose_grid",'	
 funs+='"_AT_GSM_response_grid_from_dose_grid",'	
 funs+='"_AT_run_GSM_method",'	
 funs+='"_AT_GSM_local_dose_distrib",'	
 funs+='"_AT_GSM_multiple_local_dose_distrib",'	
 
  #----AT_KatzModel.h	
 funs+='"_AT_KatzModel_KatzExtTarget_inactivation_probability",'	
 funs+='"_AT_KatzModel_CucinottaExtTarget_inactivation_probability",'	
 funs+='"_AT_KatzModel_inactivation_probability",'	
 funs+='"_AT_KatzModel_KatzExtTarget_inactivation_cross_section_integrand_m",'	
 funs+='"_AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2",'	
 funs+='"_AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_integrand_m",'	
 funs+='"_AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2",'	
 funs+='"_AT_KatzModel_inactivation_cross_section_m2",'	
 funs+='"_AT_KatzModel_KatzExtTarget_ButtsKatz_TrackWidth",'	
 funs+='"_AT_KatzModel_KatzExtTarget_Zhang_TrackWidth",'	
 funs+='"_AT_KatzModel_single_field_survival_from_inactivation_cross_section",'	
 funs+='"_AT_KatzModel_inactivation_cross_section_approximation_m2",'	
 funs+='"_AT_KatzModel_single_field_survival",'	
 funs+='"_AT_KatzModel_mixed_field_survival",'	
 funs+='"_AT_KatzModel_single_field_survival_optimized_for_fluence_vector",'	
 funs+='"_AT_P_RDD",'	
 funs+='"_AT_sI_int",'	
 funs+='"_AT_D_RDD_Gy_int",'	
 
  #----AT_PhysicsRoutines.h	
 funs+='"_AT_beta_from_E_single",'	
 funs+='"_AT_beta_from_E",'	
 funs+='"_AT_E_from_beta_single",'	
 funs+='"_AT_E_from_beta",'	
 funs+='"_AT_E_from_gamma_single",'	
 funs+='"_AT_E_from_gamma",'	
 funs+='"_AT_E_MeV_u_from_momentum_single",'	
 funs+='"_AT_E_MeV_u_from_momentum_MeV_c_u",'	
 funs+='"_AT_gamma_from_E_single",'	
 funs+='"_AT_gamma_from_E",'	
 funs+='"_AT_effective_charge_from_beta_single",'	
 funs+='"_AT_effective_charge_from_beta",'	
 funs+='"_AT_energy_straggling_MeV2_cm2_g",'	
 funs+='"_AT_energy_straggling_after_slab_E_MeV_u",'	
 funs+='"_AT_effective_charge_from_E_MeV_u_single",'	
 funs+='"_AT_effective_charge_from_E_MeV_u",'	
 funs+='"_AT_mean_excitation_energy_eV_from_Z_single",'	
 funs+='"_AT_mean_excitation_energy_eV_from_Z",'	
 funs+='"_AT_mass_correction_terms_new",'	
 funs+='"_AT_max_relativistic_E_transfer_MeV_new_single",'	
 funs+='"_AT_max_classic_E_transfer_MeV_new_single",'	
 funs+='"_AT_max_E_transfer_MeV_new_single",'	
 funs+='"_AT_max_E_transfer_MeV_new",'	
 funs+='"_AT_mass_correction_terms",'	
 funs+='"_AT_max_relativistic_E_transfer_MeV_single",'	
 funs+='"_AT_max_classic_E_transfer_MeV_single",'	
 funs+='"_AT_max_E_transfer_MeV_single",'	
 funs+='"_AT_max_E_transfer_MeV",'	
 funs+='"_AT_momentum_from_E_MeV_c_u_single",'	
 funs+='"_AT_momentum_MeV_c_u_from_E_MeV_u",'	
 funs+='"_AT_dose_Gy_from_fluence_cm2_single",'	
 funs+='"_AT_dose_Gy_from_fluence_cm2",'	
 funs+='"_AT_fluence_cm2_from_dose_Gy_single",'	
 funs+='"_AT_fluence_cm2_from_dose_Gy",'	
 funs+='"_AT_beam_par_physical_to_technical",'	
 funs+='"_AT_beam_par_technical_to_physical",'	
 funs+='"_AT_interparticleDistance_m",'	
 funs+='"_AT_inv_interparticleDistance_Gy",'	
 funs+='"_AT_single_impact_fluence_cm2_single",'	
 funs+='"_AT_single_impact_fluence_cm2",'	
 funs+='"_AT_single_impact_dose_Gy_single",'	
 funs+='"_AT_single_impact_dose_Gy",'	
 funs+='"_AT_total_D_Gy",'	
 funs+='"_AT_total_fluence_cm2",'	
 funs+='"_AT_fluence_weighted_E_MeV_u",'	
 funs+='"_AT_dose_weighted_E_MeV_u",'	
 funs+='"_AT_fluence_weighted_LET_MeV_cm2_g",'	
 funs+='"_AT_dose_weighted_LET_MeV_cm2_g",'	
 funs+='"_AT_stopping_power_ratio",'	
 funs+='"_AT_mean_number_of_tracks_contrib",'	
 funs+='"_AT_kinetic_variable_single",'	
 funs+='"_AT_Rutherford_SDCS",'	
 funs+='"_AT_Rutherford_scatter_cross_section",'	
 funs+='"_AT_gyroradius_m",'	
 
  #----AT_DataMaterial.h	
 funs+='"_AT_material_index_from_material_number",'	
 funs+='"_AT_material_name_from_number",'	
 funs+='"_AT_material_number_from_name",'	
 funs+='"_AT_density_g_cm3_from_material_no",'	
 funs+='"_AT_I_eV_from_material_no",'	
 funs+='"_AT_alpha_g_cm2_MeV_from_material_no",'	
 funs+='"_AT_p_MeV_from_material_no",'	
 funs+='"_AT_m_g_cm2_from_material_no",'	
 funs+='"_AT_average_A_from_material_no",'	
 funs+='"_AT_average_Z_from_material_no",'	
 funs+='"_AT_phase_from_material_no",'	
 funs+='"_AT_get_material_data",'	
 funs+='"_AT_get_materials_data",'	
 funs+='"_AT_electron_density_m3_from_material_no_single",'	
 funs+='"_AT_electron_density_m3_from_material_no_multi",'	
 funs+='"_AT_plasma_energy_J_from_material_no",'	
 funs+='"_AT_electron_density_m3_single",'	
 funs+='"_AT_electron_density_m3_multi",'	
 funs+='"_AT_plasma_energy_J_single",'	
 funs+='"_AT_electron_density_m3_from_composition",'	
 funs+='"_AT_average_A_from_composition",'	
 funs+='"_AT_average_Z_from_composition",'	
 funs+='"_AT_effective_Z_from_composition",'	
 funs+='"_AT_I_eV_from_composition",'	
 funs+='"_AT_set_user_material",'	
 funs+='"_AT_set_user_material_from_composition",'	
 
  #----AT_SPC.h	
 funs+='"_AT_SPC_get_number_of_bytes_in_file",'	
 funs+='"_AT_SPC_fast_read_buffer",'	
 funs+='"_AT_SPC_decompose_size",'	
 funs+='"_AT_SPC_decompose_header",'	
 funs+='"_AT_SPC_decompose_data",'	
 funs+='"_AT_SPC_get_number_of_bins_from_filename_fast",'	
 funs+='"_AT_SPC_read_header_from_filename_fast",'	
 funs+='"_AT_SPC_read_data_from_filename_fast",'	
 funs+='"_AT_SPC_read_from_filename_fast",'	
 funs+='"_AT_SPC_number_of_bins_at_range",'	
 funs+='"_AT_SPC_spectrum_at_range",'	
 
  #----AT_StoppingPowerData.h	
 funs+='"_AT_stopping_power_source_model_name_from_number",'	
 funs+='"_AT_stopping_power_source_model_number_from_name",'	
 
  #----AT_Algorithms_CPP.h	
 funs+='"_AT_run_CPPSC_method",'	
 funs+='"_AT_run_CPPSS_method",'	
 
  #----AT_RDD_ExtendedTarget.h	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_integrand_Gy",'	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration",'	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_Gy",'	
 funs+='"_AT_inverse_RDD_ExtendedTarget_KatzPoint_solver_function_Gy",'	
 funs+='"_AT_inverse_RDD_ExtendedTarget_KatzPoint_m",'	
 funs+='"_AT_RDD_ExtendedTarget_CucinottaPoint_integrand_Gy",'	
 funs+='"_AT_RDD_ExtendedTarget_CucinottaPoint_Gy_by_integration",'	
 funs+='"_AT_RDD_ExtendedTarget_CucinottaPoint_Gy",'	
 funs+='"_AT_inverse_RDD_ExtendedTarget_CucinottaPoint_solver_function_Gy",'	
 funs+='"_AT_inverse_RDD_ExtendedTarget_CucinottaPoint_m",'	
 
  #----AT_ElectronRange.h	
 funs+='"_AT_ER_ButtsKatz_range_g_cm2",'	
 funs+='"_AT_ER_Waligorski_range_g_cm2",'	
 funs+='"_AT_ER_Edmund_range_g_cm2",'	
 funs+='"_AT_ER_Geiss_range_g_cm2",'	
 funs+='"_AT_ER_Scholz_range_g_cm2",'	
 funs+='"_AT_ER_Tabata_range_g_cm2",'	
 funs+='"_AT_ER_PowerLaw_alpha",'	
 funs+='"_AT_ER_Tabata_constants",'	
 funs+='"_AT_ER_Scholz_new_range_g_cm2",'	
 funs+='"_AT_ER_AM_RadDiff_range_g_cm2",'	
 funs+='"_AT_max_electron_ranges_m",'	
 funs+='"_AT_max_electron_range_m",'	
 
  #----AT_StoppingPowerDataICRU.h	
 funs+='"_AT_ICRU_wrapper",'	
 
  #----AT_Algorithms_IGK.h	
 funs+='"_AT_run_IGK_method",'	
 
  #----AT_DataParticle.h	
 funs+='"_AT_particle_no_from_Z_and_A_single",'	
 funs+='"_AT_particle_no_from_Z_and_A",'	
 funs+='"_AT_A_from_particle_no_single",'	
 funs+='"_AT_A_from_particle_no",'	
 funs+='"_AT_atomic_weight_from_Z",'	
 funs+='"_AT_Z_from_particle_no_single",'	
 funs+='"_AT_Z_from_particle_no",'	
 funs+='"_AT_atomic_weight_from_particle_no",'	
 funs+='"_AT_I_eV_from_particle_no",'	
 funs+='"_AT_nuclear_spin_from_particle_no_multi",'	
 funs+='"_AT_nuclear_spin_from_particle_no_single",'	
 funs+='"_AT_nuclear_spin_from_Z_and_A",'	
 funs+='"_AT_particle_name_from_particle_no_single",'	
 funs+='"_AT_particle_no_from_particle_name_single",'	
 funs+='"_AT_particle_name_from_particle_no",'	
 funs+='"_AT_particle_no_from_particle_name",'	
 funs+='"_AT_Z_from_element_acronym_single",'	
 funs+='"_AT_Z_from_element_acronym",'	
 funs+='"_AT_element_acronym_from_Z_single",'	
 funs+='"_AT_element_acronym_from_Z",'	
 funs+='"_AT_atomic_weight_from_element_acronym_single",'	
 funs+='"_AT_atomic_weight_from_element_acronym",'	
 funs+='"_AT_density_g_cm3_from_element_acronym_single",'	
 funs+='"_AT_density_g_cm3_from_element_acronym",'	
 funs+='"_AT_I_eV_from_element_acronym_single",'	
 funs+='"_AT_I_eV_from_element_acronym",'	
 funs+='"_AT_electron_density_cm3_from_element_acronym_single",'	
 funs+='"_AT_electron_density_cm3_from_element_acronym",'	
 
  #----AT_RDD_Tabulated.h	
 funs+='"_AT_RDD_RadicalDiffusion_get_energy_idx",'	
 funs+='"_AT_RDD_RadicalDiffusion_Gy",'	
 funs+='"_AT_inverse_RadicalDiffusion_m",'	
 funs+='"_AT_r_min_RadicalDiffusion_m",'	
 funs+='"_AT_r_max_RadicalDiffusion_m",'	
 funs+='"_AT_d_min_RadicalDiffusion_Gy",'	
 funs+='"_AT_d_max_RadicalDiffusion_Gy",'	
 funs+='"_AT_E_RadicalDiffusion_MeV_u",'	
 funs+='"_AT_n_bins_RadicalDiffusion",'	
 
  #----AT_DataRange.h	
 funs+='"_AT_Stopping_Power_Mass_MeV_cm2_g_int",'	
 funs+='"_AT_CSDA_range_g_cm2_multi",'	
 funs+='"_AT_CSDA_range_g_cm2_single",'	
 funs+='"_AT_CSDA_range_difference_solver",'	
 funs+='"_AT_CSDA_energy_after_slab_E_MeV_u_single",'	
 funs+='"_AT_CSDA_energy_after_slab_E_MeV_u_multi",'	
 funs+='"_AT_WEPL_multi",'	
 funs+='"_AT_WEPL_single",'	
 
  #----AT_RDD_Simple.h	
 funs+='"_AT_RDD_Test_Gy",'	
 funs+='"_AT_inverse_RDD_Test_m",'	
 funs+='"_AT_RDD_Geiss_Gy",'	
 funs+='"_AT_inverse_RDD_Geiss_m",'	
 funs+='"_AT_RDD_Katz_coeff_Gy",'	
 funs+='"_AT_RDD_Katz_coeff_Gy_general",'	
 funs+='"_AT_RDD_Katz_LinearER_Dpoint_Gy",'	
 funs+='"_AT_RDD_Katz_PowerLawER_Dpoint_Gy",'	
 funs+='"_AT_RDD_KatzPoint_Gy",'	
 funs+='"_AT_inverse_RDD_KatzPoint_LinearER_m",'	
 funs+='"_AT_inverse_RDD_KatzPoint_PowerLawER_solver_function_Gy",'	
 funs+='"_AT_inverse_RDD_KatzPoint_m",'	
 funs+='"_AT_RDD_Cucinotta_f_shortRange",'	
 funs+='"_AT_RDD_Cucinotta_f_longRange",'	
 funs+='"_AT_RDD_Cucinotta_Ddelta_Gy",'	
 funs+='"_AT_RDD_Cucinotta_Dexc_Gy",'	
 funs+='"_AT_RDD_CucinottaPoint_Gy",'	
 funs+='"_AT_inverse_RDD_Cucinotta_solver_function_Gy",'	
 funs+='"_AT_inverse_RDD_Cucinotta_m",'	
 
  #----AT_RDD_ShellAveraged.h	
 funs+='"_AT_RDD_Katz_LinearER_Daverage_Gy",'	
 funs+='"_AT_RDD_Katz_PowerLawER_DaverageKernel",'	
 funs+='"_AT_RDD_Katz_PowerLawER_DaverageKernel_approx",'	
 funs+='"_AT_RDD_Katz_PowerLawER_Daverage_Gy",'	
 funs+='"_AT_RDD_Cucinotta_Ddelta_average_integrand_m",'	
 funs+='"_AT_RDD_Cucinotta_Ddelta_average_Gy",'	
 funs+='"_AT_RDD_Cucinotta_Dexc_average_Gy",'	
 funs+='"_AT_RDD_Cucinotta_Cnorm",'	
 funs+='"_AT_RDD_Geiss_average_Gy",'	
 funs+='"_AT_RDD_Katz_LinearER_dEdx_J_m",'	
 funs+='"_AT_RDD_Katz_PowerLawER_dEdx_J_m",'	
 funs+='"_AT_RDD_Katz_LinearER_DSite_Gy",'	
 funs+='"_AT_RDD_Katz_PowerLawER_DSite_Gy",'	
 funs+='"_AT_RDD_KatzSite_Gy",'	
 funs+='"_AT_inverse_RDD_KatzSite_m"'	
 
  funs+=']'


emcc libat.a libgsl.a -o libat.html -s WASM=$WASM -s EXPORTED_FUNCTIONS="$funs" -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]'

mkdir -p ../output/
rm ../output/*
cp libat.a ../output/
cp libat.html ../output/
cp libat.wasm ../output/ 2>/dev/null || : #ignore error when build with -s WASM=0
cp libat.js ../output/

cd ..
rm -r _build
