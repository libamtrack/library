import common
import generation_library
import logging
import ConfigParser
import pylab
from os import listdir, getcwd, mkdir
import datetime
import math


def read_config(filename):
    """Read config file and return dictionary of settings.
    If file is not found - exit with error status.
    TODO description
    """
    # check if configuration file exists
    try:
        f = open(filename)
        f.close()
    except IOError:
        #logging.error('Config file \'' + filename + '\' not found. Exiting.')
        exit(1)

    #logging.info('Parsing config file \'' + filename + '\'...')
    config = ConfigParser.RawConfigParser()
    config.read(filename)

    list_of_name_value_pairs = []

    for section in config.sections():
        list_for_current_section = config.items(section)
        list_of_name_value_pairs.extend( list_for_current_section )

    options_dict = dict(list_of_name_value_pairs)
    
    return options_dict


def merge_options_dictionaries_fav_first(dictionary1, dictionary2):
    """Merge two dictionaries. If same value is set in both
    choose the one from first. Return result dictionary.
    """
    result_dictionary = dictionary2

    for key in dictionary1.keys():
        if dictionary1[key] == None:
            continue
        else:
            result_dictionary[key] = dictionary1[key]

    return result_dictionary


def check_options_consistency(options):
    """Check if there are conflicts in options
    (passed as a dictionary).
    """
    logging.info('Checking options consistency...')

    error = False

    # provide Gottshalck rule OR number of base BPs
    if 'use_gottshalck_80_rule' in options.keys() and 'number_of_base_bps' in options.keys():
        if options['use_gottshalck_80_rule'] == 'True' and int(options['number_of_base_bps']) > 0:
            logging.error('Conflict in options. Either specify number_of_base_bps OR tell to use_gottshalck_80_rule.')
            error = True

    if not 'use_gottshalck_80_rule' in options.keys() or options['use_gottshalck_80_rule'] == 'False':
        if not 'number_of_base_bps' in options.keys() or int(options['number_of_base_bps']) < 2:
            logging.error('Conflict in options. Number of base BPs is not specified.')
            error = True

    # check if datafile is provided
    if not 'bp_database' in options.keys():
        logging.error('Conflict in options. File with BPs database not provided.')
        error = True

    # check if file with fluence drop coefficients is provided
    if not 'fluence_drop_coefficients_file' in options.keys():
        logging.error('Lack in options. File with fluence drop curve coefficients not provided.')
        error = True

    # check if plexi to water coefficient is provided
    if not 'mod_plexi_to_water' in options.keys():
        logging.error('Lack in options. Coefficient for changing water range to modulator plexi range not provided.')
        error = True

    # check if plexi to water coefficient is provided
    if not 'rs_plexi_to_water' in options.keys():
        logging.error('Lack in options. Coefficient for changing water range to rs plexi range not provided.')
        error = True

    # check if both correction and adjusting plateau were specified
    if common.check_option('fit_adjusting_plateau') and 'correction_to_spread' in common.options.keys() and float(common.options['correction_to_spread']) != 0:
        logging.warning('Manual correction to spread was made when fitting with adjusting plateau.')

    if error:
        logging.error('There were conflicts in options. Exiting.')
        exit(1)
    else:
        logging.info('Options consistency OK.')



def write_info(modulator_name, version, filename, RD_IN, SOBP_begin, SOBP_end, SOBP_proximal_100, SOBP_distal_100, SOBP_distal_10, chi2_at_minimum, sobp_begin_percentage, sobp_end_percentage, newline_character = '\n', separator = ' '):
    """TODO
    """
    logging.info('Writing info output to file \'' + filename + '\'...')

    # prepare output as list of lines
    list_of_lines = []

    list_of_lines.append('Modulator name:' + separator + modulator_name + newline_character)
    list_of_lines.append('Script version:' + separator + version + newline_character)
    now = datetime.datetime.now()
    date = str(now.day) + '.' + str(now.month) + '.' + str(now.year)
    date += ' ' + str(now.hour) + ':' + str(now.minute)
    list_of_lines.append('Date created:' + separator + date + newline_character)
    list_of_lines.append('RD_IN:' + separator + str(RD_IN) + newline_character)
    list_of_lines.append('SOBP begin (proximal ' + '%.1f' % sobp_begin_percentage + '%):' + separator + '%.3f' % SOBP_begin + newline_character)
    list_of_lines.append('SOBP proximal 100%:' + separator + '%.3f' % SOBP_proximal_100 + newline_character)
    list_of_lines.append('SOBP distal 100%:' + separator + '%.3f' % SOBP_distal_100 + newline_character)
    list_of_lines.append('SOBP end (distal ' + '%.1f' % sobp_end_percentage + '%):' + separator + '%.3f' % SOBP_end + newline_character)
    list_of_lines.append('SOBP distal 10%:' + separator + '%.3f' % SOBP_distal_10 + newline_character)
    list_of_lines.append('SOBP distal falloff:' + separator + '%.3f' % (SOBP_distal_10 - SOBP_end) + newline_character)

    # check if able to write to filename and write
    try:
        f = open(filename, 'w')
        f.writelines(list_of_lines)
        f.close()
        logging.info('Writing output done.')
    except IOError:
        logging.error('Cannot write output to \'' + filename + '\'. Exiting.')
        exit(1)
    

def write_BP(modulator_name, version, filename, BP_data, newline_character = '\n', separator = ' '):
    """TODO
    """
    logging.info('Writing BP output to file \'' + filename + '\'...')

    # prepare output as list of lines
    list_of_lines = []

    for (x,y) in BP_data:
        if (not math.isnan(x)) and (not math.isnan(y)):
            list_of_lines.append(str(x) + separator + str(y) + newline_character)

    # check if able to write to filename and write
    try:
        f = open(filename, 'w')
        f.writelines(list_of_lines)
        f.close()
        logging.info('Writing output done.')
    except IOError:
        logging.error('Cannot write output to \'' + filename + '\'. Exiting.')
        exit(1)



def create_output_folder(name):
    """Create output folder for modulator named name.
    Name is in format 'name_date_version'. If folder
    already exists number is added at the end.

    Parameters:
        name: string representing name of modulator

    Returns:
        folder_name: name of created folder
    """
    # get current date
    now = datetime.datetime.now()
    date = str(now.day) + '-' + str(now.month) + '-' + str(now.year)

    folder_name = name + '_' + date + '_v' + common.version

    if folder_name in listdir(getcwd()):
        additional_number = 1
        while folder_name + '_' + str(additional_number) in listdir(getcwd()):
            additional_number += 1
        folder_name += '_' + str(additional_number)
        
    mkdir(folder_name)
    return folder_name
