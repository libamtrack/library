import os
from shutil import copy
import logging
import argparse
import csv

import src.common as common
import src.pyamtrack
import src.pyamtrack_SPC

from src.in_out import *
from src.generation_library import *

def main():
        
    # argument parsing:
    # http://docs.python.org/library/argparse.html#module-argparse
    description_string = """ Script for designing modulators. Options are
    read from configuration file 'plot.cfg'. Options passed directly
    as program parameters have higher priority.
    """
    parser = argparse.ArgumentParser(description = description_string)

    parser.add_argument('--verbose', action = 'store_true',
                        help = 'be verbose')
    parser.add_argument('--debug', action = 'store_true',
                        help = 'display debug messages')
    parser.add_argument('--version', '-V', action = 'store_true',
                        help = 'show version and exit')

    parser.add_argument('--logfile', action = 'store',
                        help = 'where write logs to', default = 'generation.log')
    parser.add_argument('--cfgfile', action = 'store',
                        help = 'config file', default = 'plot.cfg')
    parser.add_argument('--name', '-n', action = 'store',
                        help = 'name of modulator being generated')
    


    # doctest
    import doctest
    failure_count = doctest.testmod(common)[0]
    failure_count += doctest.testmod(generation_library)[0]
    failure_count += doctest.testmod()[0]
    if failure_count > 0 :
        exit(-1)

    args_dict = common.parse_command_line_arguments(parser)

    config_dict = read_config(filename = args_dict['cfgfile'])
    common.options = merge_options_dictionaries_fav_first(args_dict, config_dict)

    # create folder for output
    # output folder = name of the folder with system folder separator at end
    output_folder = create_output_folder(common.options['name']) + os.sep

    # copy configuration file
    copy(args_dict['cfgfile'], output_folder + args_dict['cfgfile'])

    common.setup_logging(filename = output_folder + common.options['logfile'],
                  write_on_screen = args_dict['verbose'],
                  display_debug = args_dict['debug'],
                  logger_name = '')

    logging.debug('Options: %s', common.options)
   
    
    return
    

def save_datafile( filename, columns):
    writer = csv.writer(open( filename, "w"), delimiter=' ')
    zipped = zip(*columns)
    writer.writerows(zipped)


def executeParallelMap( functionName, argsT, inputV ):
    outputV = []

    # Create jobserver
    job_server = pp.Server( secret = 'abcd' )
    ncpus = int(common.options['ncpu'])     
    job_server.set_ncpus(ncpus)
    parts = ncpus
    step = (len(inputV)) / parts + 1
    jobs = []
    start_time = time.time()
             
    for index in xrange(parts):
	starti = index*step
	endi = min((index+1)*step, len(inputV))
        jobs.append(job_server.submit(functionName, args=(tuple([inputV[starti:endi]]) + argsT), depfuncs=(), modules=("src.pyamtrack_SPC",), globals = None))
    
    for job in jobs:
	outputV.extend(job())

	print "Time elapsed: ", time.time() - start_time, "s"                              
	print job_server.print_stats()        

    job_server.destroy()

    return outputV


def plot_Dose(filename, spectrum, maximum, base_BP_positions, coefficients, two_beams, plot_format = 'png', show = True):
    logging.info('Plotting Dose... ')

    # clear previous plots
    pylab.clf()

    # first take wide range to determine where should plot end
    dx = 0.01
    X_wide = [ n * dx for n in range( int((base_BP_positions[-1] + 4) / dx)) ]

    serial = common.check_option('serial')
    parallel = common.check_option('parallel')

    if parallel:
       coefficients_p = []
       for c in coefficients:
       	coefficients_p.append(float(c))

    Y_SOBP = []
           
    if two_beams:
        X_wide = [ n * dx for n in range( int((base_BP_positions[0] + base_BP_positions[-1]) / dx)) ]
        plateau_dist = base_BP_positions[-1]
        plateau_prox = base_BP_positions[0]

	def calculateYa(X_vect, maximum, base_BP_positions, spectrum, plateau_dist, plateau_prox, coefficients):
	        from src import pyamtrack_SPC
	        Y_vect = [sum( [coefficients[i] * (pyamtrack_SPC.dose_at_depth(x, maximum, base_BP_positions[i], spectrum) + 
                                           pyamtrack_SPC.dose_at_depth(plateau_dist +  plateau_prox - x, maximum, base_BP_positions[i], spectrum))
                                          for i in range(len(coefficients)) ] ) for x in X_vect]
		return Y_vect

    	if serial:
		Y_SOBP = calculateYa(X_wide, maximum, base_BP_positions, spectrum, plateau_dist, plateau_prox, coefficients)

	if parallel:
		args = (maximum, base_BP_positions, spectrum, plateau_dist, plateau_prox, coefficients_p)
		Y_SOBP = executeParallelMap( calculateYa, args, X_wide)

    else:

	def calculateYb(X_vect, maximum, base_BP_positions, spectrum, coefficients):
	        from src import pyamtrack_SPC
	        Y_vect = [sum( [coefficients[i] * pyamtrack_SPC.dose_at_depth(x, maximum, base_BP_positions[i], spectrum)
                                          for i in range(len(coefficients)) ] ) for x in X_vect]
		return Y_vect

    	if serial:
		Y_SOBP = calculateYb(X_wide, maximum, base_BP_positions, spectrum, coefficients)

	if parallel:
		args = (maximum, base_BP_positions, spectrum, coefficients_p)
		Y_SOBP = executeParallelMap( calculateYb, args, X_wide)

    # plot line at 1
    pylab.axhline(1, color = 'g')

    pylab.plot(X_wide, Y_SOBP, color = 'r', label = 'calculated SOBP')

    listPartPeaksY = []
    for i in range(len(coefficients)):
	def calculateYc(X_vect, coefficient, maximum, base_BP_position, spectrum):
	        from src import pyamtrack_SPC
		return [coefficient * pyamtrack_SPC.dose_at_depth(x, maximum, base_BP_position, spectrum) for x in X_vect]
	if serial:
		Yi = calculateYc(X_wide, coefficients[i], maximum, base_BP_positions[i], spectrum)
	if parallel:
		args = (coefficients_p[i], maximum, base_BP_positions[i], spectrum)
		Yi = executeParallelMap( calculateYc, args, X_wide)

        listPartPeaksY.append(Yi)
        pylab.plot(X_wide, Yi, color = 'g', label = 'BP nr ' + str(i))
        if two_beams:
            plateau_dist = base_BP_positions[-1]
            plateau_prox = base_BP_positions[0]
            Yis = [coefficients[i] * pyamtrack_SPC.dose_at_depth(plateau_dist +  plateau_prox - x, maximum, base_BP_positions[i], spectrum) for x in X_wide]
            pylab.plot(X_wide, Yis, color = 'g', label = 'BP nr ' + str(i) + ' bis')
        
                
    pylab.xlim(0.,X_wide[-1])

    pylab.grid(True)

    # save plot to file
    if filename != None:
        logging.info('Saving SOBP plot to file\'' + filename + '.' + plot_format + '\'...')
        pylab.savefig(filename + '.' + plot_format, format = plot_format)
        save_datafile(filename + '.dat', [X_wide,Y_SOBP]+listPartPeaksY)

    # show what has been plotted
    if show:
        pylab.show()


def plot_Survival(filename, scaling_factor, spectrum, maximum, base_BP_positions, coefficients, Katz_coeffs, er_model, two_beams, plot_format = 'png', show = True):
    logging.info('Plotting Survival... ')

    # clear previous plots
    pylab.clf()
    
    m, D0, sigma, kappa = Katz_coeffs

    # first take wide range to determine where should plot end
    dx = 0.05
    X_wide = [ n * dx for n in range( int((base_BP_positions[-1] + 4) / dx)) ]
    if two_beams:
        X_wide = [ n * dx for n in range( int((base_BP_positions[0] + base_BP_positions[-1]) / dx)) ]
    Y_wide = []
        
    def calculateY( X_vect, coefficients, maximum, base_BP_positions, spectrum, scaling_factor, two_beams, m, D0, sigma, kappa, er_model):
      Y_vect = []
      from src import pyamtrack_SPC
      for x in X_vect:
        E_MeV_u_total = []
        particle_no_total = []
        fluence_cm2_total = []
        for i in range(len(coefficients)):
            shift = maximum - base_BP_positions[i]
            E_MeV_u, particle_no, fluence_cm2 = pyamtrack_SPC.spectrum_at_depth(x + shift, spectrum)
            fluence_cm2_coef = [f * coefficients[i] * scaling_factor for f in fluence_cm2]
            E_MeV_u_total += E_MeV_u
            particle_no_total += particle_no
            fluence_cm2_total += fluence_cm2_coef
            
            if two_beams:                
                plateau_dist = base_BP_positions[-1]
                plateau_prox = base_BP_positions[0]
                y = plateau_dist +  plateau_prox - x
                
                E_MeV_u_second, particle_no_second, fluence_cm2_second = pyamtrack_SPC.spectrum_at_depth(y + shift, spectrum)
                fluence_cm2_coef_second = [f * coefficients[i] * scaling_factor for f in fluence_cm2_second] # TODO check                                
                E_MeV_u_total += E_MeV_u_second
                particle_no_total += particle_no_second
                fluence_cm2_total += fluence_cm2_coef_second

            
        y = pyamtrack_SPC.survival(E_MeV_u_total, particle_no_total, fluence_cm2_total, m, D0, sigma, kappa, er_model)
        Y_vect.append(y)
      return Y_vect

    serial = common.check_option('serial')
    parallel = common.check_option('parallel')

    # Create jobserver
    job_server = pp.Server( secret = 'abcd' )
    ncpus = int(common.options['ncpu'])
    job_server.set_ncpus(ncpus)

    if serial:
       Y_wide = calculateY( X_wide, coefficients, maximum, base_BP_positions, spectrum, scaling_factor, two_beams, m, D0, sigma, kappa, er_model)

    if parallel:
	coefficients_p = []
	for c in coefficients:
      		coefficients_p.append(float(c))
	args = (coefficients_p, maximum, base_BP_positions, spectrum, scaling_factor, two_beams, m, D0, sigma, kappa, er_model)
	Y_wide = executeParallelMap( calculateY, args, X_wide)

    # plot line at 1
    pylab.axhline(1, color = 'g')

    pylab.plot(X_wide, Y_wide, color = 'r', label = 'survival')
               
    pylab.xlim(0.,X_wide[-1])

    pylab.grid(True)

    # save plot to file
    if filename != None:
        logging.info('Saving survival plot to file\'' + filename + '.' + plot_format + '\'...')
        pylab.savefig(filename + '.' + plot_format, format = plot_format)
        save_datafile(filename + '.dat', [X_wide,Y_wide])

    # show what has been plotted
    if show:
        pylab.show()


# MAIN FUNCTION
if __name__ == "__main__":
    main()

