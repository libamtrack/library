import numpy as np
import pyamtrack.libAT as libam

def GSM_test(E_MeV_u, fluence_cm2):
    
    # this part is probably okay
    particle_no = libam.AT_particle_no_from_Z_and_A_single(1, 1) # Z, A
    material_no = libam.AT_material_number_from_name("Water, Liquid")
    rdd_model_no = libam.RDDModels["RDD_Geiss"].value
    a0_m = 10e-9
    rdd_parameter = [a0_m, 0., 0.]
    er_model = libam.AT_ERModels["ER_Edmund"].value
    stopping_power_source_no = libam.stoppingPowerSource_no["PSTAR"].value
    
    # this part is probably wrong ...
    nX = 100 # what is nX?
    pixel_size_m = 1e-6 # too coarse, just a test ..
    number_of_bins = 100 # hmm?
    dose_bin_centers = 10 ** np.linspace(np.log10(0.001), np.log10(1), number_of_bins)
    random_number_generator_seed = [2705490069]
    zero_dose_fraction = [0]
    dose_frequency_Gy = [0] * number_of_bins
    

    libam.AT_GSM_local_dose_distrib(p_E_MeV_u=[E_MeV_u], 
                                    p_fluence_cm2=[fluence_cm2], 
                                    p_particle_no=[particle_no],
                                    p_material_no=material_no, 
                                    p_rdd_model=rdd_model_no,
                                    p_rdd_parameter=rdd_parameter,
                                    p_er_model=er_model,
                                    p_stopping_power_source_no=stopping_power_source_no,
                                    p_nX=nX,
                                    p_pixel_size_m=pixel_size_m,
#                                    p_number_of_bins=number_of_bins,
                                    p_dose_bin_centers_Gy=dose_bin_centers.tolist(),
                                    p_random_number_generator_seed=random_number_generator_seed,
                                    p_zero_dose_fraction=zero_dose_fraction,
                                    p_dose_frequency_Gy=dose_frequency_Gy,
                                   )

if __name__ == '__main__':
    GSM_test(60., 1e6)