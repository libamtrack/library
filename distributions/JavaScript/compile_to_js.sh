#!/bin/bash

if [ "$1" != "" ]; then
    echo "WASM parameter set to $1"
    WASM=$1
else
    echo "WASM parameter EMPTY !!! Default set to 1"
    WASM=1
fi

rm -rf build 

emcmake cmake -DGSL_INCLUDE_DIRS=$GSL_INCLUDE_DIRS -DGSL_LIBRARY=$GSL_LIBRARY -DGSL_CBLAS_LIBRARY=$GSL_CBLAS_LIBRARY -S . -B build || exit 1
emmake cmake --build build || exit 1

funs='['	
 

 #----AT_DataParticle.h
 funs+='"_AT_atomic_weight_from_particle_no_single",'

  #----AT_EnergyLoss.h
 funs+='"_AT_Vavilov_energy_loss_distribution",'
 funs+='"_AT_Gauss_energy_loss_distribution",'
 funs+='"_AT_Landau_energy_loss_distribution",'	
 funs+='"_AT_energy_loss_distribution",'	
 
  #----AT_StoppingPower.h	
 funs+='"_AT_Mass_Stopping_Power",'	
 funs+='"_AT_Stopping_Power",'	
 funs+='"_AT_Mass_Stopping_Power_with_no",'	
 funs+='"_AT_Stopping_Power_with_no",'	
 funs+='"_AT_Energy_MeV_u_from_Stopping_Power_single",'	

 #----AT_ProtonAnalyticalModels.h
 funs+='"_AT_dose_Bortfeld_Gy_multi",'
 funs+='"_AT_dose_Bortfeld_Gy_single",'
 funs+='"_AT_LET_t_Wilkens_keV_um_multi",'
 funs+='"_AT_LET_d_Wilkens_keV_um_multi",'
 funs+='"_AT_proton_RBE_multi",'

 #----AT_ProtonAnalyticalBeamParameters.h
 funs+='"_AT_range_Bortfeld_cm",'
 funs+='"_AT_fwhm_Bortfeld_cm",'
 funs+='"_AT_max_plateau_Bortfeld",'
 funs+='"_AT_energy_Bortfeld_MeV",'
 funs+='"_AT_fit_Bortfeld",'
 
  #----AT_KatzModel_Implementation.h
 funs+='"_AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2",'
 funs+='"_AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2",'
 funs+='"_AT_KatzModel_inactivation_cross_section_m2",'

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
 funs+='"_AT_momentum_from_E_MeV_c_u_single",'	
 funs+='"_AT_momentum_MeV_c_u_from_E_MeV_u",'	
 funs+='"_AT_dose_Gy_from_fluence_cm2_single",'	
 funs+='"_AT_dose_Gy_from_fluence_cm2",'	
 funs+='"_AT_fluence_cm2_from_dose_Gy_single",'	
 funs+='"_AT_fluence_cm2_from_dose_Gy",'	
 
  #----AT_RDD_ExtendedTarget.h	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_integrand_Gy",'	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_Gy_by_integration",'	
 funs+='"_AT_RDD_ExtendedTarget_KatzPoint_Gy",'	
 funs+='"_AT_RDD_ExtendedTarget_CucinottaPoint_Gy_by_integration",'	
 funs+='"_AT_RDD_ExtendedTarget_CucinottaPoint_Gy",'	
 
  #----AT_ElectronRange.h	
 funs+='"_AT_max_electron_ranges_m",'	
 funs+='"_AT_max_electron_range_m",'	
 
  #----AT_DataRange.h	
 funs+='"_AT_CSDA_range_g_cm2_multi",'	
 funs+='"_AT_CSDA_range_g_cm2_single",'	
 funs+='"_AT_CSDA_energy_after_slab_E_MeV_u_single",'	
 funs+='"_AT_CSDA_energy_after_slab_E_MeV_u_multi",'	
 
  #----AT_D_RDD_Gy.h	
 funs+='"_AT_D_RDD_Gy"'

  funs+=']'

emcc build/libat.a $GSL_LIBRARY $GSL_CBLAS_LIBRARY -o output/libat.html -sWASM=$WASM -sEXPORTED_FUNCTIONS="$funs" -sEXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]'  || exit 1
