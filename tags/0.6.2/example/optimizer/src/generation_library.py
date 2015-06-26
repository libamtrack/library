import bisect
import pylab
import scipy.optimize
import logging
import numpy
import math
import sys

import common
from src import pyamtrack_SPC

import pp
import time

def calculate_base_BP_positions_given_number(first_BP_position, last_BP_position, n_of_base_BPs):
	"""Calculate base BP positions so that n_of_base_BPs are
	evenly spread between fist_BP_position and last_BP_position.

	Parameters:
		first_BP_position: position of first base BP in mm in water
		last_BP_position: position of last base BP in mm in water
		n_of_base_BPs: desired number of base BPs
	Returns:
		base_BP_positions: a list of BPs positions in mm (in water).
	"""
	position_step = (last_BP_position - first_BP_position) / n_of_base_BPs
	positions = []
	for n in range(n_of_base_BPs):
		positions.append(last_BP_position - n * position_step)
	positions.append(first_BP_position)

	# return sorted values
	return sorted(positions)
	

# TODO add a function to calculate base BP position uniformly for given distance between two adjacent BPs


def shift_BP(BP_to_be_shifted, desired_range):
	"""Returns BP as a list of (x,y) pairs with maximum shifted by given amount.
	Distance should be given assumed there's no range shifter.

	Parameters:
		BP_to_be_shifted: BP to be shifted (list of (x,y) pairs), will be modified
		desired_range: desired shift in range in mm in water

	TODO doctest missing
	"""

	shifted_BP = [ (XY[0] + desired_range, XY[1]) for XY in BP_to_be_shifted]

	return shifted_BP



# we count how many times function which is minimized is called
# this gives us insight into minimizing algorithm performance
call_count = 0
call_count_glob = 0


def fit(algorithm, local_max, base_BP_positions, spectrum, value_at_plateau, plateau_prox, plateau_dist, mesh_size, precision, initial_coefficients=[], Katz_coeff=[], two_beams=False, bounds=[], fit_interrupt_bound=0, plateau_shape_coefs=[1., 0., 0., 0.]):
	"""Base fitting procedure. It minimizes sum squared distances between
	SOBP and 1 (chi^2). Sum is over mesh on given plateau.
	Fitting is aborted when SOBP plateau is in area between
	1 - fit_interrupt_bound and 1 + fit_interrupt_bound.
	
	Parameters:
		algorithm: Algorithm for fitting. Only fmin_l_bfgs_b is available.
		base_BPs: List of base BPs (each BP is (x,y) pair). SOBP if formed
			from these BPs.
		plateau_prox, plateau_dist: SOBP is optimized to be small between
			these values.
		mesh_size: Distance between points of mesh. chi^2 is calculated
			only on mesh points.
		precision: how small should be DERIVATIVE of chi^2 (times machine
			precision) for fit to end. 1 is a extremely good fit, 1e10 -
			mediocre, 1e20 - poor fit.
		initial coefficients (default: [0,0,...,1]): starting point for
			fitting algorithm.
		fit_interrupt_bound (default: 0): Fitting is interrupted when whole
			SOBP's plateau is in area around 1 of this size.

	Returns:
		coefficients: list of optimal coefficients found
		chi2: chi2 value with coefficients found
	"""
	logging.info('Fitting...')

	if plateau_prox > plateau_dist:
		logging.error('In arguments to fit function plateau_prox is greater than plateau_dist. Exiting.')
		exit(1)

	# starting point for fitting algorithm
	coeffs = initial_coefficients
	
	# points where chi^2 will be evaluated
	function_evaluation_mesh = [ plateau_prox + n * mesh_size
								 for n in
								 range(int((plateau_dist - plateau_prox) / mesh_size) + 1)]
	
	print "mesh: ", function_evaluation_mesh
	print "base pos: ", base_BP_positions

	logging.debug('Fit mesh start: %f, fit mesh end: %f',
				  function_evaluation_mesh[0],
				  function_evaluation_mesh[-1])

	def prepareEmptyDict(x, spectrum, local_max, positions):
		pn = sorted([1001, 5011, 2004, 4009, 6012, 3007])
		d = dict.fromkeys(pn, {})
		
		xmin = x + local_max - max(positions)
		xmax = x + local_max - min(positions)
		
		depth_min_ind = bisect.bisect_left(sorted(spectrum.keys()), xmin)
		depth_max_ind = bisect.bisect_left(sorted(spectrum.keys()), xmax)
		
		for depth in sorted(spectrum.keys())[depth_min_ind - 1:depth_max_ind + 1]:
			E_MeV_u, particle_no, fluence_cm2 = pyamtrack_SPC.spectrum_at_depth(depth, spectrum)
			for j in range(len(fluence_cm2)):
				d[particle_no[j]][E_MeV_u[j]] = 0
		return d
	
	
	class FittingDone(Exception):
		"""Exception raised when fitting should be aborted.
		"""
		def __init__(self, coefficients, minimum, maximum, chi2):
			self.coefficients = coefficients
			self.minimum = minimum
			self.maximum = maximum
			self.chi2 = chi2
			
	def prepareLists(x, spectrum, local_max, positions, coefficients_of_base_BPs):
		from src import pyamtrack_SPC
		particle_no_total = []
		E_MeV_u_total = []
		fluence_cm2_total = []

  
		for i in range(len(coefficients_of_base_BPs)):
			shift = local_max - positions[i]
			E_MeV_u, particle_no, fluence_cm2 = pyamtrack_SPC.spectrum_at_depth(x + shift, spectrum)
			fluence_cm2_coef = [f * coefficients_of_base_BPs[i] for f in fluence_cm2]					
			particle_no_total.extend(particle_no)
			E_MeV_u_total.extend(E_MeV_u)
			fluence_cm2_total.extend(fluence_cm2_coef)
		
		return E_MeV_u_total, particle_no_total, fluence_cm2_total


	def chi2simple(coefficients_of_base_BPs, args):
		"""sum of square differences between plateau and sum_of_base_BPs
		in fixed values of x

		chi2simple(C), where coefficients_of_base_BPs = [C1, C2, C3, ...]
		
		TODO doctest missing
		"""
		total = 0
		# maximum and minimum are used for interrupt bound
		maximum = -1e20
		minimum = 1e20
				
		value_at_plateau, local_max, positions, Katz_coeffs, two_beams, spectrum, er_model = args
		
		m, D0, sigma, kappa = Katz_coeffs

		coefficients = []
		for c in coefficients_of_base_BPs:
			coefficients.append(float(c))

		def partial_sum(spectrum, function_evaluation_mesh, value_at_plateau, local_max, positions, coefficients_of_base_BPs, plateau_dist, plateau_prox, two_beams, m, D0, sigma, kappa, er_model):
			"""Calculates partial sum"""
		
			from src import pyamtrack_SPC
			maximum = -100
			minimum = 100

			total = 0
			for x in function_evaluation_mesh:
				E_MeV_u_total, particle_no_total, fluence_cm2_total = prepareLists(x, spectrum, local_max, positions, coefficients_of_base_BPs)

			##	two opposite beams
				if two_beams:
					y = plateau_dist + plateau_prox - x
					E_MeV_u_total_second, particle_no_total_second, fluence_cm2_total_second = prepareLists(y, spectrum, local_max, positions, coefficients_of_base_BPs)
					E_MeV_u_total.extend(E_MeV_u_total_second)
					particle_no_total.extend(particle_no_total_second)
					fluence_cm2_total.extend(fluence_cm2_total_second)
						
				survival_at_depth_x = pyamtrack_SPC.survival(E_MeV_u_total, particle_no_total, fluence_cm2_total, m, D0, sigma, kappa, er_model)			
			
				if value_at_plateau > 0 and survival_at_depth_x > 0:
					total += (math.log10(value_at_plateau) - math.log10(survival_at_depth_x)) ** 2
				else:
					total = 1e100
			
				if survival_at_depth_x > maximum:
					 maximum = survival_at_depth_x
				if survival_at_depth_x < minimum:
					 minimum = survival_at_depth_x 

			return total, minimum, maximum

		maximum = -100
		minimum = 100
		if parallel:
			parts = ncpus
			step = ((len(function_evaluation_mesh)) / parts) + 1

			jobs = []
			start_time = time.time()

			for index in xrange(parts):
				starti = index * step
				endi = min((index + 1) * step, len(function_evaluation_mesh))			
				jobs.append(job_server.submit(partial_sum, args=(spectrum, function_evaluation_mesh[starti:endi], value_at_plateau, local_max, positions, coefficients, plateau_dist, plateau_prox, two_beams, m, D0, sigma, kappa, er_model), depfuncs=(prepareLists,), modules=("src.pyamtrack_SPC", "math",)))

			for job in jobs:
				total_p, min_p, max_p = job()
				total += total_p
				if max_p > maximum:
					maximum = max_p
				if min_p < minimum:
					minimum = min_p

#			print "Time elapsed: ", time.time() - start_time, "s"							  
#			print job_server.print_stats()

		if serial:
			total_ser, minimum, maximum = partial_sum(spectrum, function_evaluation_mesh, value_at_plateau, local_max, positions, coefficients, plateau_dist, plateau_prox, two_beams, m, D0, sigma, kappa, er_model)
		
		if serial and not parallel:
			total = total_ser

		global call_count
		call_count += 1
		
		logging.info('%i calls, chi^2 = %.12f, max = %f, min = %f',
						 call_count, total, maximum, minimum)
		print coefficients


		if (maximum <= value_at_plateau + fit_interrupt_bound) and (minimum >= value_at_plateau - fit_interrupt_bound):
			raise FittingDone(coefficients=coefficients_of_base_BPs,
							minimum=minimum,
							maximum=maximum,
							chi2=total)

		return total


	def chi2withDerivative(coefficients_of_base_BPs, args):
		"""Sum of square differences between plateau and sum_of_base_BPs
		in fixed values of x (specified in function_evaluation_mesh).
		If at any point sum of BPs with coefficients passed is constrained
		around 1 in area less then specified in fit_interrupt_bound exception
		is thrown.

		chi2withDerivative(C), where coefficients_of_base_BPs = [C1, C2, C3, ...]
		
		TODO doctest missing
		"""
		total = 0
		# maximum and minimum are used for interrupt bound
		maximum = 0
		minimum = 100
		grad = [0] * len(coefficients_of_base_BPs)
				
		value_at_plateau, local_max, positions, spectrum, plateau_shape_coefs = args
		
		print "positions : ", positions
		coef_str = "[ "
		for c in coefficients_of_base_BPs:
			coef_str += str(c) + " , "
		coef_str += " ] "
		print "coefficients : ", coef_str
		
		coefficients = []
		for c in coefficients_of_base_BPs:
			coefficients.append(float(c))


		def partial_sum(spectrum, function_evaluation_mesh, local_max, positions, coefficients_of_base_BPs, plateau_shape_coefs):
			"""Calculates partial sum"""

	   		
	   		a0, a1, a2, a3 = plateau_shape_coefs

			grad = [0] * len(coefficients_of_base_BPs)

			maximum = -100
			minimum = 100
						
			from src import pyamtrack_SPC
			
			total = 0			
			for x in function_evaluation_mesh:
				dose_at_depth_x = [pyamtrack_SPC.dose_at_depth(x, local_max, positions[i], spectrum) for i in range(len(coefficients_of_base_BPs))]
			
				sum_of_BPs_for_this_x = sum([coefficients_of_base_BPs[i] * dose_at_depth_x[i] for i in range(len(coefficients_of_base_BPs)) ])

				value_at_plateau = a3 * x * x * x + a2 * x * x + a1 * x + a0

				total += (value_at_plateau - sum_of_BPs_for_this_x) ** 2

				for i in range(len(coefficients_of_base_BPs)):
					grad[i] += -2.0 * (value_at_plateau - sum_of_BPs_for_this_x) * dose_at_depth_x[i]

				if sum_of_BPs_for_this_x - value_at_plateau > maximum:
					maximum = sum_of_BPs_for_this_x - value_at_plateau
				if sum_of_BPs_for_this_x - value_at_plateau < minimum:
					minimum = sum_of_BPs_for_this_x - value_at_plateau

			return total, grad, minimum, maximum

		maximum = -100
		minimum = 100

		if parallel:
			parts = ncpus
			step = ((len(function_evaluation_mesh)) / parts) + 1

			jobs = []
			start_time = time.time()

			for index in xrange(parts):
				starti = index * step
				endi = min((index + 1) * step, len(function_evaluation_mesh))
				jobs.append(job_server.submit(partial_sum, args=(spectrum, function_evaluation_mesh[starti:endi], local_max, positions, coefficients, plateau_shape_coefs), modules=("src.pyamtrack_SPC", "math",)))

			for job in jobs:
				total_p, grad_p, min_p, max_p = job()
				total += total_p		
				for i in range(len(grad_p)):
					grad[i] += grad_p[i]
				if max_p > maximum:
					maximum = max_p
				if min_p < minimum:
					minimum = min_p

			print "Time elapsed: ", time.time() - start_time, "s"
		 	print job_server.print_stats()
			
		if serial:
			total_ser, grad_ser, minimum, maximum = partial_sum(spectrum, function_evaluation_mesh, local_max, positions, coefficients, plateau_shape_coefs)
		
		if serial and not parallel:
			total = total_ser
			grad = grad_ser
							   
		global call_count
		call_count += 1
		#if call_count % 10 == 0:
		logging.info('\t %i calls, chi^2 = %.12f, max = %f, min = %f\n',
						 call_count, total, maximum, minimum)

		if maximum < fit_interrupt_bound and (-minimum) < fit_interrupt_bound:
			raise FittingDone(coefficients=coefficients_of_base_BPs,
							  minimum=minimum,
							  maximum=maximum,
							  chi2=total)

		return total, numpy.array(grad) # TODO remove numpy array

	# choose fit algorithm
	
	
	if algorithm == 'dose':
		try:

			serial = common.check_option('serial')
			parallel = common.check_option('parallel')

			# Create jobserver
			job_server = pp.Server(secret='abcd')
			ncpus = int(common.options['ncpu'])
			job_server.set_ncpus(ncpus)

			# http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html
			opt_coeffs, chi2_at_minimum, information_dictionary = scipy.optimize.fmin_l_bfgs_b(
				func=chi2withDerivative,
				x0=coeffs,
				args=[[value_at_plateau, local_max, base_BP_positions, spectrum, plateau_shape_coefs]],
				approx_grad=False,
				bounds=bounds,
				factr=precision,
				maxfun=2000)
	
			job_server.destroy()
			
		except FittingDone as result:
			opt_coeffs = result.coefficients
			chi2_at_minimum = result.chi2
			job_server.destroy()
			logging.debug('Min when fitting: %f', result.minimum)
			logging.debug('Max when fitting: %f', result.maximum)
	elif algorithm == 'survival':
		try:
			
			serial = common.check_option('serial')
			parallel = common.check_option('parallel')

			# Create jobserver
			job_server = pp.Server(secret='abcd')
			ncpus = int(common.options['ncpu'])
			er_model = int(common.options['er_model'])
			job_server.set_ncpus(ncpus)
			
			opt_coeffs, chi2_at_minimum, information_dictionary = scipy.optimize.fmin_l_bfgs_b(
				func=chi2simple,
				x0=coeffs,
				args=[[value_at_plateau, local_max, base_BP_positions, Katz_coeff, two_beams, spectrum, er_model]],
				approx_grad=True,
				bounds=bounds,
				factr=precision,
				maxfun=2000,
				epsilon=1e-3,
				iprint=0)
			
			job_server.destroy()
						
		except FittingDone as result:
			opt_coeffs = result.coefficients
			chi2_at_minimum = result.chi2
			job_server.destroy()
			logging.debug('Min when fitting: %f', result.minimum)
			logging.debug('Max when fitting: %f', result.maximum)
	else:
		logging.error('Algorithm \'' + algorithm + '\' is not known in fit() function. Exiting.')
		exit(1)


	logging.info('Fitting done. %i calls total.', call_count)

	return opt_coeffs, chi2_at_minimum
