#!/bin/bash

if [ "$1" != "" ]; then
    echo "WASM parameter set to $1"
    WASM=$1
else
    echo "WASM parameter EMPTY !!! Default set to 1"
    WASM=1
fi

cd .
mkdir _build
cd _build
cp ../libgsl.a .
cp ../libgslcblas.a .
ls -al .

emcmake cmake .. -DGSL_INCLUDE_DIRS=$GSL_INCLUDE_DIRS -DGSL_LIBRARY=$GSL_LIBRARY -DGSL_CBLAS_LIBRARY=$GSL_CBLAS_LIBRARY
emmake make -j4

funs='['	
 
  #----AT_Error.h	

  #----AT_StoppingPowerDataBethe.h	

  #----AT_MultipleCoulombScattering.h	

  #----AT_RDD.h	

  #----AT_StoppingPowerDataPSTAR.h	

  #----AT_EnergyLoss.h
 funs+='"_AT_mean_energy_loss_keV",'	
 funs+='"_AT_xi_keV",'	
 funs+='"_AT_kappa_single",'	
 funs+='"_AT_Landau_PDF",'	
 funs+='"_AT_Landau_IDF",'	
 funs+='"_AT_lambda_landau_from_energy_loss_multi",'	
 funs+='"_AT_lambda_landau_from_energy_loss_single",'	
 funs+='"_AT_lambda_mean_single",'	
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
 
  #----AT_GammaResponse.h	

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
 funs+='"_AT_energy_Bortfeld_MeV",'
 funs+='"_AT_fit_Bortfeld",'

  #----AT_SuccessiveConvolutions.h

 #----AT_KatzModel.h
 funs+='"_AT_KatzModel_sigma_um2_single",'
 funs+='"_AT_KatzModel_sigma_um2",'

  #----AT_KatzModel_Implementation.h
 funs+='"_AT_KatzModel_KatzExtTarget_inactivation_probability",'	
 funs+='"_AT_KatzModel_CucinottaExtTarget_inactivation_probability",'	
 funs+='"_AT_KatzModel_inactivation_probability",'	
 funs+='"_AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2",'	
 funs+='"_AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2",'	
 funs+='"_AT_KatzModel_inactivation_cross_section_m2",'	
 funs+='"_AT_KatzModel_KatzExtTarget_ButtsKatz_TrackWidth",'	
 funs+='"_AT_KatzModel_KatzExtTarget_Zhang_TrackWidth",'	
 funs+='"_AT_KatzModel_single_field_survival_from_inactivation_cross_section",'	
 funs+='"_AT_KatzModel_inactivation_cross_section_approximation_m2",'	
 funs+='"_AT_KatzModel_single_field_survival",'	

  #----AT_PhysicsRoutines.h	
 funs+='"_AT_E_MeV_u_from_E_MeV",'	
 funs+='"_AT_E_MeV_from_E_MeV_u",'	
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

  #----AT_DataMaterial.h	

  #----AT_StoppingPowerData.h	

  #----AT_Algorithms_CPP.h	

  #----AT_RDD_ExtendedTarget.h	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_integrand_Gy",'	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration",'	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_Gy",'	
 funs+='"_AT_inverse_RDD_ExtendedTarget_KatzPoint_m",'	
 funs+='"_AT_RDD_ExtendedTarget_CucinottaPoint_Gy_by_integration",'	
 funs+='"_AT_RDD_ExtendedTarget_CucinottaPoint_Gy",'	
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
 
  #----AT_DataParticle.h	
 funs+='"_AT_atomic_weight_from_particle_no_single",'	
 funs+='"_AT_I_eV_from_particle_no",'	
 funs+='"_AT_nuclear_spin_from_particle_no_single",'	
 
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
 funs+='"_AT_RDD_Geiss_Gy",'	
 funs+='"_AT_inverse_RDD_Geiss_m",'	
 funs+='"_AT_RDD_Katz_coeff_Gy",'	
 funs+='"_AT_RDD_Katz_coeff_Gy_general",'	
 funs+='"_AT_RDD_Katz_LinearER_Dpoint_Gy",'	
 funs+='"_AT_RDD_Katz_PowerLawER_Dpoint_Gy",'	
 funs+='"_AT_RDD_KatzPoint_Gy",'	
 funs+='"_AT_inverse_RDD_KatzPoint_LinearER_m",'	
 funs+='"_AT_inverse_RDD_KatzPoint_m",'	
 funs+='"_AT_RDD_Cucinotta_f_shortRange",'	
 funs+='"_AT_RDD_Cucinotta_f_longRange",'	
 funs+='"_AT_RDD_Cucinotta_Ddelta_Gy",'	
 funs+='"_AT_RDD_Cucinotta_Dexc_Gy",'	
 funs+='"_AT_RDD_CucinottaPoint_Gy",'	

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

  funs+=']'

emcc libat.a libgsl.a libgslcblas.a -o libat.html -s WASM=$WASM -s EXPORTED_FUNCTIONS="$funs" -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]'

mkdir -p ../output/
rm ../output/*
cp libat.a ../output/
cp libat.html ../output/
cp libat.wasm ../output/ 2>/dev/null || : #ignore error when build with -s WASM=0
cp libat.js ../output/

cd ..
rm -r _build
