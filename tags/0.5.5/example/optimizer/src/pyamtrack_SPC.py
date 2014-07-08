import src.pyamtrack
import os
import bisect

def extract_spectrum_at_range( spcdirpath, range_cmg):
    
    number_of_bins = src.pyamtrack.AT_SPC_number_of_bins_at_range(spcdirpath, range_cmg)
    print "Number of bins " + str(number_of_bins)
                
    if number_of_bins <= 0:
        print "Could not get number of bins "
        return None
    
    depth_step = [0] * number_of_bins
    depth_g_cm2 = [0] * number_of_bins
    E_MeV_u = [0] * number_of_bins
    DE_MeV_u =  [0] * number_of_bins
    particle_no =  [0] * number_of_bins
    fluence_cm2 = [0] * number_of_bins

    result = src.pyamtrack.AT_SPC_spectrum_at_range(spcdirpath, range_cmg, number_of_bins, depth_step, depth_g_cm2, E_MeV_u, DE_MeV_u, particle_no, fluence_cm2)
    
    if result[0] != 0:
        print "Could not get read data"
        return None
    
    print result
    
    depth_step, depth_g_cm2, E_MeV_u, DE_MeV_u, particle_no, fluence_cm2 = result[1:]    
    
    return depth_step, depth_g_cm2, E_MeV_u, DE_MeV_u, particle_no, fluence_cm2

def extract_range(filename):
    
    E_MeV_u = 0.
    peak_position_g_cm2 = 0.
    particle_no = 0
    material_no = 0
    normalisation = 1.
    steps_no = 1
    
    result = src.pyamtrack.AT_SPC_read_header_from_filename_fast(filename, E_MeV_u, peak_position_g_cm2, particle_no, material_no, normalisation, steps_no)
    
    return result[2]

    
def scale_fluence_to_input_dose( input_dose_Gy, spectrum_dict):
    
    material_no = 1
    stopping_power_source_no = 0
    
    min_depth = min(spectrum_dict.keys())
    spectrum_input_channel = spectrum_at_depth(min_depth, spectrum_dict)
        
    E_MeV_u_input, particle_no_input, fluence_cm2_input = spectrum_input_channel
    
    old_input_dose_Gy = src.pyamtrack.AT_total_D_Gy(len(E_MeV_u_input), E_MeV_u_input, particle_no_input, fluence_cm2_input, material_no, stopping_power_source_no)
    
    factor = input_dose_Gy / old_input_dose_Gy
    
    new_spectrum_dict = spectrum_dict
    
    for depth in spectrum_dict.keys():
        fluence_cm2 = spectrum_dict[depth][2]
        scaled_fluence_cm2 = [ f * factor for f in fluence_cm2]
        new_spectrum_dict[depth] = spectrum_dict[depth][0], spectrum_dict[depth][1], scaled_fluence_cm2
    
    return new_spectrum_dict


def scale_fluence_to_maximum_dose( maximum_dose_Gy, maximum_position, spectrum_dict):
    
    material_no = 1
    stopping_power_source_no = 0
    
    spectrum_at_maximum = spectrum_at_depth(maximum_position, spectrum_dict)
    
    E_MeV_u_input, particle_no_input, fluence_cm2_input = spectrum_at_maximum
    
    old_input_dose_Gy = src.pyamtrack.AT_total_D_Gy(len(E_MeV_u_input), E_MeV_u_input, particle_no_input, fluence_cm2_input, material_no, stopping_power_source_no)
    
    factor = maximum_dose_Gy / old_input_dose_Gy
    
    new_spectrum_dict = spectrum_dict
    
    for depth in spectrum_dict.keys():
        fluence_cm2 = spectrum_dict[depth][2]
        scaled_fluence_cm2 = [ f * factor for f in fluence_cm2]
        new_spectrum_dict[depth] = spectrum_dict[depth][0], spectrum_dict[depth][1], scaled_fluence_cm2
    
    return new_spectrum_dict

    
    
def get_spectrum_dictionary_step(spectrum):
    depth_step = spectrum[0]
    depth_step_uniq = sorted(list(set(depth_step)))

    spect_dict = {}
    for step in depth_step_uniq:
        depth_step_new, depth_g_cm2_new, E_MeV_u_new, DE_MeV_u_new, particle_no_new, fluence_cm2_new = spectrum_at_depth_step(step, spectrum) 
        spect_dict[step] =  depth_g_cm2_new, E_MeV_u_new, particle_no_new, fluence_cm2_new
    
    return spect_dict

def get_spectrum_dictionary_depth(spectrum):
    depth_step_v, depth_g_cm2, E_MeV_u, DE_MeV_u, particle_no, fluence_cm2 = spectrum
    depth_g_cm2_uniq = sorted(list(set(depth_g_cm2)))

    spect_dict = {}
    for depth in depth_g_cm2_uniq:
        indexes_for_given_depth = [i for i in range(len(depth_g_cm2)) if depth_g_cm2[i] == depth]
        E_MeV_u_new = [E_MeV_u[i] for i in indexes_for_given_depth]
        particle_no_new = [particle_no[i] for i in indexes_for_given_depth]
        fluence_cm2_new = [fluence_cm2[i] for i in indexes_for_given_depth]
        spect_dict[depth] =  E_MeV_u_new, particle_no_new, fluence_cm2_new
    
    return spect_dict

    
def spectrum_at_depth_step(depth_step, spectrum):
    depth_step_v, depth_g_cm2, E_MeV_u, DE_MeV_u, particle_no, fluence_cm2 = spectrum
    
    indexes_for_given_depth = [i for i in range(len(depth_step_v)) if depth_step_v[i] == depth_step]    
         
    depth_step_new = [depth_step] * len(indexes_for_given_depth)
    depth_g_cm2_new = [depth_g_cm2[indexes_for_given_depth[0]]] * len(indexes_for_given_depth)
    E_MeV_u_new = [E_MeV_u[i] for i in indexes_for_given_depth]
    DE_MeV_u_new = [DE_MeV_u[i] for i in indexes_for_given_depth]
    particle_no_new = [particle_no[i] for i in indexes_for_given_depth]
    fluence_cm2_new = [fluence_cm2[i] for i in indexes_for_given_depth]

    return depth_step_new, depth_g_cm2_new, E_MeV_u_new, DE_MeV_u_new, particle_no_new, fluence_cm2_new
    
def spectrum_at_depth(depth_g_cm2, spectrum):   
    result = None
    
    if depth_g_cm2 in spectrum:
        result = spectrum[depth_g_cm2]
    elif depth_g_cm2 > max(spectrum.keys()):
        result = [],[],[]
    else:
        depth_ind = bisect.bisect_left(sorted(spectrum.keys()),depth_g_cm2)
        depth_next = sorted(spectrum.keys())[depth_ind]
	ratio = 1
	if depth_ind > 0:
	        depth_prev = sorted(spectrum.keys())[depth_ind-1]
	        ratio = (depth_g_cm2-depth_prev)/(depth_next-depth_prev)
	else:
	        depth_prev = sorted(spectrum.keys())[0]
        E_MeV_u, particle_no, fluence_cm2_prev = spectrum[depth_prev]
        fluence_cm2_next = spectrum[depth_next][2]
        fluence_cm2 = [0] * len(fluence_cm2_next)
        for i in range(len(fluence_cm2_prev)):
            fluence_cm2[i] = (1.0-ratio)*fluence_cm2_prev[i] + ratio*fluence_cm2_next[i]
        result = E_MeV_u, particle_no, fluence_cm2
        
    return result
    
    
def dose_at_depth(depth_g_cm2, local_max, maximum, spectrum):

    material_no = 1
    stopping_power_source_no = 0
    
    shift = local_max - maximum

    E_MeV_u, particle_no, fluence_cm2 = spectrum_at_depth(depth_g_cm2 + shift, spectrum)
            
    non_zero_fluence_indexes = [i for i in range(len(fluence_cm2)) if fluence_cm2[i] > 0]    
    E_MeV_u_OK = [E_MeV_u[i] for i in non_zero_fluence_indexes]
    particle_no_OK = [particle_no[i] for i in non_zero_fluence_indexes]
    fluence_cm2_OK = [fluence_cm2[i] for i in non_zero_fluence_indexes]
    
    dose_Gy = src.pyamtrack.AT_total_D_Gy(len(E_MeV_u_OK), E_MeV_u_OK, particle_no_OK, fluence_cm2_OK, material_no, stopping_power_source_no)

    return dose_Gy

def fluence_at_depth(depth_g_cm2, local_max, maximum, spectrum):

    material_no = 1
    stopping_power_source_no = 0
    
    shift = local_max - maximum

    E_MeV_u, particle_no, fluence_cm2 = spectrum_at_depth(depth_g_cm2 + shift, spectrum)

    non_zero_fluence_indexes = [i for i in range(len(fluence_cm2)) if fluence_cm2[i] > 0]    
    E_MeV_u_OK = [E_MeV_u[i] for i in non_zero_fluence_indexes]
    particle_no_OK = [particle_no[i] for i in non_zero_fluence_indexes]
    fluence_cm2_OK = [fluence_cm2[i] for i in non_zero_fluence_indexes]
    
    res_fluence_cm2 = {}

    particle_no_uniq = set(particle_no_OK)

    for p in particle_no_uniq:
	    res_fluence_cm2[p] = 0

    for i in range(len(fluence_cm2_OK)):
	res_fluence_cm2[particle_no_OK[i]] += fluence_cm2_OK[i]

    return res_fluence_cm2
    
        
def survival(E_MeV_u, particle_no, fluence_cm2,  m_number_of_targets, D0_characteristic_dose_Gy, sigma0_m2, kappa, er_model = 2 ):
    
    non_zero_fluence_indexes = [i for i in range(len(fluence_cm2)) if fluence_cm2[i] > 0]    
    E_MeV_u_OK = [E_MeV_u[i] for i in non_zero_fluence_indexes]
    particle_no_OK = [particle_no[i] for i in non_zero_fluence_indexes]
    fluence_cm2_OK = [fluence_cm2[i] for i in non_zero_fluence_indexes]

        
    number_of_items = len(E_MeV_u_OK)
    
    material_no = 1
    rdd_model = 6
    rdd_parameters = [1e-10,1e-8,1e-10,0]
    use_approximation = True
    stopping_power_source_no = 0
    survival_v = 0.
    code, survival = src.pyamtrack.AT_KatzModel_mixed_field_survival ( number_of_items, 
                                                               fluence_cm2_OK, 
                                                               E_MeV_u_OK, 
                                                               particle_no_OK, 
                                                               material_no, 
                                                               rdd_model, 
                                                               rdd_parameters, 
                                                               er_model, 
                                                               D0_characteristic_dose_Gy, 
                                                               m_number_of_targets, 
                                                               sigma0_m2, 
                                                               use_approximation, 
                                                               kappa, 
                                                               stopping_power_source_no, 
                                                               survival_v )
    
    return survival
    
