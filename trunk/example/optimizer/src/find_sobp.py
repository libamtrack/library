import logging
import argparse

import matplotlib
matplotlib.use('PDF')

import src.common as common
import src.generation_library as generation_library 
import src.in_out as in_out
import src.plotting as plotting
import src.pyamtrack
import src.pyamtrack_SPC

from src.in_out import *
from src.generation_library import *
from src.plotting import *


def sampleDoseFromBPdictionary( spectrum ):
    result = [ (depth, src.pyamtrack_SPC.dose_at_depth(depth, 0, 0, spectrum)) for depth in sorted(spectrum.keys())]
    return result

def printBP( BP ):
    print "#depth dose"
    for (depth, dose_Gy) in BP:
        print depth, dose_Gy
    

def removeZeroFluenceItemsFromSpectrum( spectrum ):
    depth_step, depth_g_cm2, E_MeV_u, DE_MeV_u, particle_no, fluence_cm2 = spectrum
        
    non_zero_fluence_indexes = [i for i in range(len(fluence_cm2)) if fluence_cm2[i] > 0]    
    E_MeV_u_OK = [E_MeV_u[i] for i in non_zero_fluence_indexes]
    DE_MeV_u_OK = [DE_MeV_u[i] for i in non_zero_fluence_indexes]
    particle_no_OK = [particle_no[i] for i in non_zero_fluence_indexes]
    fluence_cm2_OK = [fluence_cm2[i] for i in non_zero_fluence_indexes]
    depth_step_OK =  [depth_step[i] for i in non_zero_fluence_indexes]
    depth_g_cm2_OK =  [depth_g_cm2[i] for i in non_zero_fluence_indexes]
    
    return [depth_step_OK, depth_g_cm2_OK, E_MeV_u_OK, DE_MeV_u_OK, particle_no_OK, fluence_cm2_OK]


def main():

    # argument parsing:
    # http://docs.python.org/library/argparse.html#module-argparse
    description_string = """ Script for designing modulators. Options are
    read from configuration file 'modulatory.cfg'. Options passed directly
    as program parameters have higher priority.
    """
    parser = argparse.ArgumentParser(description = description_string)

    parser.add_argument('--range', '-r', action = 'store',
                        help = 'range of SOBP in mm')
    parser.add_argument('--spread', '-s', action = 'store',
                        help = 'width of SOBP in mm')
    parser.add_argument('--mesh', '-m', action = 'store',
                        help = 'size of mesh used when minimizing')
    parser.add_argument('--precision', '-p', action = 'store',
                        help = 'precision for fitting algorithm')
    parser.add_argument('--interrupt_bound', action = 'store',
                        help = 'interrupt fitting when SOBP is flat within this range')

    parser.add_argument('--verbose', action = 'store_true',
                        help = 'be verbose')
    parser.add_argument('--debug', action = 'store_true',
                        help = 'display debug messages')
    parser.add_argument('--version', '-V', action = 'store_true',
                        help = 'show version and exit')
    parser.add_argument('--noplot', action = 'store_true',
                        help = 'do not generate plots')

    parser.add_argument('--logfile', action = 'store',
                        help = 'where write logs to', default = 'generation.log')
    parser.add_argument('--cfgfile', action = 'store',
                        help = 'config file', default = 'modulatory.cfg')
    parser.add_argument('--name', '-n', action = 'store',
                        help = 'name of modulator being generated')
    

    # doctest
    import doctest
    failure_count = doctest.testmod(common)[0]
    failure_count += doctest.testmod(generation_library)[0]
    failure_count += doctest.testmod(in_out)[0]
    failure_count += doctest.testmod(plotting)[0]
    failure_count += doctest.testmod()[0]
    if failure_count > 0 :
        exit(-1)

    args_dict = common.parse_command_line_arguments(parser)

    config_dict = read_config(filename = args_dict['cfgfile'])
    common.options = merge_options_dictionaries_fav_first(args_dict, config_dict)

    # create folder fo output
    # output folder = name of the folder with system folder separator at end
    output_folder = create_output_folder(common.options['name']) + os.sep

    # copy configuration file
    copy(args_dict['cfgfile'], output_folder + args_dict['cfgfile'])

    common.setup_logging(filename = output_folder + common.options['logfile'],
                  write_on_screen = args_dict['verbose'],
                  display_debug = args_dict['debug'],
                  logger_name = '')

    logging.debug('Options: %s', common.options)

    SOBP_range = float(common.options['range'])
    SOBP_spread = float(common.options['spread'])

    # GENERATION PROCESS
        
    spectrum_full = src.pyamtrack_SPC.extract_spectrum_at_range(common.options['spcdir'], SOBP_range)
        
    #spectrum = removeZeroFluenceItemsFromSpectrum( spectrum_full )
    spectrum = spectrum_full

    spectrum_dict = src.pyamtrack_SPC.get_spectrum_dictionary_depth(spectrum)

    BPsampled = sampleDoseFromBPdictionary( spectrum_dict )
    maximum = src.common.calculate_BP_maximum_position( BPsampled )
    logging.debug('BP maximum: %s', maximum)
    #printBP(BPsampled)

    base_BP_positions = calculate_base_BP_positions_given_number(
            first_BP_position = SOBP_range - SOBP_spread,
            last_BP_position = SOBP_range,
            n_of_base_BPs = int(common.options['number_of_base_bps']))
    logging.debug('Selected base BP positions: %s', base_BP_positions)

    scaled_spectrum_dict = {}
    initial_coefs = []
    bounds = []
    value_at_plateau = None
    
    # scale fluence to input dose
    if common.options['algorithm'] == 'survival':
        input_dose_Gy = 2.
        scaled_spectrum_dict = src.pyamtrack_SPC.scale_fluence_to_input_dose( input_dose_Gy, spectrum_dict)
        initial_coefs = [0.5/len(base_BP_positions)] * len(base_BP_positions)
        initial_coefs[-1] *= 8
        initial_coefs[-2] *= 5
        bounds = [[0.1/len(base_BP_positions),10.0/len(base_BP_positions)]] * len(base_BP_positions)
        value_at_plateau = 0.2
    elif common.options['algorithm'] == 'dose':
        maximum_dose_Gy = 1.
        scaled_spectrum_dict = src.pyamtrack_SPC.scale_fluence_to_maximum_dose( maximum_dose_Gy, maximum, spectrum_dict)
        initial_coefs = [1./len(base_BP_positions)] * len(base_BP_positions)
        initial_coefs[-1] = 1.
        bounds = [[0.0,1.3]] * len(base_BP_positions)
        value_at_plateau = 1.
    else:
        logging.debug('Wrong algorithm type: %s', common.options['algorithm'])
        return 1
        # scale fluence to maximum dose for flat dose problem ??
    

    # optimize both base BPs coefficients and their positions
    # big_mesh_size = (base_BP_positions[1] - base_BP_positions[0])/4.
    big_mesh_size = float(common.options['mesh'])
    
    print "Starting from ", initial_coefs
    Katz_coefs = [float(common.options['m']),float(common.options['d0']), float(common.options['sigma0']), float(common.options['kappa'])] 
    
    plateau_shape_coefs  = [float(common.options['plateau_shape_a0']),float(common.options['plateau_shape_a1']), float(common.options['plateau_shape_a2']), float(common.options['plateau_shape_a3'])]
    
    first_coefficients, chi2 = fit(
            algorithm = common.options['algorithm'],
            local_max = maximum,
            base_BP_positions = base_BP_positions,
            spectrum = scaled_spectrum_dict,
            value_at_plateau = value_at_plateau, 
            initial_coefficients = initial_coefs,
            Katz_coeff = Katz_coefs,
            two_beams = common.check_option('two_beams'),
            bounds = bounds,
            plateau_prox = SOBP_range - SOBP_spread,
            plateau_dist = SOBP_range,
            mesh_size = big_mesh_size,
            precision = float(common.options['precision']),
            fit_interrupt_bound = float(common.options['interrupt_bound']),
            plateau_shape_coefs = plateau_shape_coefs)

    coefficients = first_coefficients

    logging.info('Coefficients found: %s', coefficients)
    
    # PLOTTING
    if common.check_option('draw_plots') and common.options['noplot'] == False:
        src.plotting.plot_Dose(output_folder + "/dose", scaled_spectrum_dict, maximum, base_BP_positions, coefficients, common.check_option('two_beams'), 'pdf', common.check_option('show_plots'))
        factor = 1.
        er_model = int(common.options['er_model'])
        src.plotting.plot_Survival(output_folder + "/survival", factor, scaled_spectrum_dict, maximum, base_BP_positions, coefficients, Katz_coefs, er_model, common.check_option('two_beams'), 'pdf', common.check_option('show_plots'))
        
    logging.info('END OF SCRIPT')


# MAIN FUNCTION
if __name__ == "__main__":
    main()

