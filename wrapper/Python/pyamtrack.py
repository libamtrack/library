#! /usr/bin/env python

'''   pyamtrack.py
   =========

   Created on: 16.02.2010
   Author: herrmann
   
   Python class for interfacing AmTrack functions. Not all functions are fully implemented.
   
   Requirements: Python2.6

   Copyright 2006, 2010 Steffen Greilich / the libamtrack team

   This file is part of the AmTrack program (libamtrack.sourceforge.net).

   AmTrack is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   AmTrack is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with AmTrack (file: copying.txt).
   If not, see <http://www.gnu.org/licenses/>
'''

#TODO: Boolean are just available in Python 2.6, warning if using 2.5 or earlier
#TODO: SPISSto be made consistent in main ibrary

import ctypes
import sys

__status__ = 'Prototype'


class AmTrack(object):
    def __init__ (self):
        print '\npyamtrack\n-----\nAt the moment just single particle fields supported\n'
        
        self.libamtrack = ctypes.cdll.LoadLibrary("libamtrack.so")
        

    def AT_SPIFF (self, n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model, RDD_parameters, 
                  ER_model, ER_parameters, gamma_model, gamma_parameters,  N2, fluence_factor, write_output,  
                  shrink_tails, shrink_tails_under, adjust_N2, lethal_events_mode):
        ''' 
        Computes HCP response and RE/RBE using compound Poison process and successive convolutions (CPP_SC, the 'SPIFF' algorithm)

        @param  n      number of particle types in the mixed particle field  (integer)
        @param  E_MeV_u      energy of particles in the mixed particle field (float list of length 'n')
        @param  particle_no    type of the particles in the mixed particle field (integer list of length 'n')
        @see          AT_DataParticle.h for definition
        @param  fluence_cm2    fluences for the given particles, doses in Gy if negative  (float list of length 'n')
        @param  material_no    index number for detector material
        (pointer to single variable)
        @see          AT_DataMaterial.h for definition
        @param  RDD_model    index number for chosen radial dose distribution
        (pointer to single variable)
        @param  RDD_parameters    parameters for chosen radial dose distribution
        (pointer to array of size depending on chosen model)
        @see          AT_RDD.h for definition
        @param  ER_model    index number for chosen electron-range model
        (pointer to single variable)
        @param  ER_parameters    parameters for chosen electron-range model
        (pointer to array of size depending on chosen model)
        @see          AT_ElectronRange.h for definition
        @param  gamma_model    index number for chosen gamma response
        (pointer to single variable)
        @param  gamma_parameters  parameters for chosen gamma response
        (pointer to array of size depending on chosen model)
        @see          AT_GammaResponse.h for definition
        @param  N2      (algorithm specific) number of bins per factor of two in
        local dose array (pointer to single variable)
        @param  fluence_factor    factor to scale the fluences given as
        "fluence_cm2" with (pointer to single variable)
        @param  write_output    if true, a protocol is written to "SuccessiveConvolutions.txt"
        in the working directory (pointer to single variable)
        @param  shrink_tails    (algorithm specific) if true, tails of the local dose
        distribution, contributing less than "shrink_tails_under" are cut
        (pointer to single variable)
        @param  shrink_tails_under  (algorithm specific) limit for tail cutting
        in local dose distribution (pointer to single variable)
        @param  adjust_N2    (algorithm specific) if true, "N2" will be increase
        if necessary at high fluence to ensure sufficient binning resolution
        @param  lethal_events_mode (algorithm specific) if true, allows to do
        calculations for cell survival
        @param  results      pointer to array of size 10 to be allocated by the
        user which will be used to return the results
        results[0]    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D
        results[1]    d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm
        results[2]    S_HCP           (algorithm independent)  absolute particle response
        results[3]    S_gamma         (algorithm independent)  absolute gamma response
        results[4]    not used        (algorithm independent)
        results[5]    u               (algorithm specific)     mean number of tracks contributing to representative point
        results[6]    u_start         (algorithm specific)     low starting value for mean number of tracks, where linearisation is applied
        results[7]    n_convolutions  (algorithm specific)     number of convolutions performed 
        results[8]    not used        (algorithm specific)
        results[9]    not used        (algorithm specific)
        @return  none
        '''
        #conversion of Python variables to C type variables
        n_ctype =                  ctypes.byref(ctypes.c_int(n))
        if n ==1 :
            E_MeV_u_ctype =            ctypes.byref(ctypes.c_float(E_MeV_u[0]))
            particle_no_ctype =        ctypes.byref(ctypes.c_int(particle_no[0]))
            fluence_cm2_ctype =        ctypes.byref(ctypes.c_float(fluence_cm2[0]))           
        else:  # if more than one particle more sophisticated variable construction needed
            print '\nat the moment just single particle fields implemented, stay updated\n'
            sys.exit(0)
#TODO: find solution to hand over arrays to C funciton

#            n_array_c_int = ctypes.c_int * n
#            n_array_c_float = ctypes.c_float* n 
#            E_MeV_tmp =n_array_c_float()
#            particle_no_tmp = n_array_c_int()
#            fluence_cm2_tmp = n_array_c_float()
#            for i, energy  in enumerate(E_MeV_u):
#                E_MeV_tmp[i] = energy
#                particle_no_tmp[i] =  particle_no[i]
#                fluence_cm2_tmp[i] = fluence_cm2[i]
#            E_MeV_u_ctype = ctypes.pointer(E_MeV_tmp)
#            particle_no_ctype = ctypes.pointer(particle_no_tmp)
#            fluence_cm2_ctype = ctypes.pointer(fluence_cm2_tmp)
           
        material_no_ctype =        ctypes.byref(ctypes.c_int(material_no))
        RDD_model_ctype =          ctypes.byref(ctypes.c_int(RDD_model))
        RDD_parameters_ctype =     ctypes.byref(ctypes.c_float(RDD_parameters))
        ER_model_ctype =           ctypes.byref(ctypes.c_int(ER_model))
        ER_parameters_ctype =      ctypes.byref(ctypes.c_float(ER_parameters))
        gamma_model_ctype =        ctypes.byref(ctypes.c_int(gamma_model))
#TODO: IMPLEMENT gamma_parameters dimension
        gamma_parameters_ctype =   ctypes.byref(ctypes.c_float(gamma_parameters))
        N2_ctype =                 ctypes.byref(ctypes.c_int(N2))
        fluence_factor_ctype =     ctypes.byref(ctypes.c_float(fluence_factor))
        write_output_ctype =       ctypes.byref(ctypes.c_bool(write_output)) 
        shrink_tails_ctype =       ctypes.byref(ctypes.c_bool(shrink_tails))
        shrink_tails_under_ctype = ctypes.byref(ctypes.c_float(shrink_tails_under))
        adjust_N2_ctype =          ctypes.byref(ctypes.c_bool(adjust_N2))
        lethal_events_mode_ctype = ctypes.byref(ctypes.c_bool(lethal_events_mode))
        tenfloats = ctypes.c_float*10
        results = tenfloats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
        
        c_at_spiff = self.libamtrack.AT_SPIFF
        c_at_spiff.restype = ctypes.c_float *10
        c_at_spiff_output = c_at_spiff(n_ctype,  E_MeV_u_ctype,  particle_no_ctype,  fluence_cm2_ctype, material_no_ctype,
                                       RDD_model_ctype, RDD_parameters_ctype, ER_model_ctype, ER_parameters_ctype,
                                       gamma_model_ctype, gamma_parameters_ctype, N2_ctype, fluence_factor_ctype,
                                       write_output_ctype, shrink_tails_ctype, shrink_tails_under_ctype, adjust_N2_ctype,
                                       lethal_events_mode_ctype, results)
        AT_SPIFF_output = []
        for item in results:
            AT_SPIFF_output.append(item)
        return AT_SPIFF_output


    def AT_SPISS (self,  n, E_MeV_u, particle_no, fluence_cm2,  material_no, RDD_model, RDD_parameters, ER_model,
                   ER_parameters, gamma_model, gamma_parameters, n_runs, N2, fluence_factor, write_output, importance_sampling):
        ''' 
        Computes HCP response and RE/RBE using compound Poison process and
        statistical sampling (CPP_SS, the 'SPISS' algorithm)

        @param  n      number of particle types in the mixed particle field
        (pointer to single variable)
        @param  E_MeV_u      energy of particles in the mixed particle field
        (pointer to array of size n)
        @param  particle_no    type of the particles in the mixed particle field
        (pointer to array of size n)
        @see          AT_DataParticle.h for definition
        @param  fluence_cm2    fluences for the given particles, doses in Gy if
        negative (pointer to array of size n)
        @param  material_no    index number for detector material
        (pointer to single variable)
        @see          AT_DataMaterial.h for definition
        @param  RDD_model    index number for chosen radial dose distribution
        (pointer to single variable)
        @param  RDD_parameters    parameters for chosen radial dose distribution
        (pointer to array of size depending on chosen model)
        @see          AT_RDD.h for definition
        @param  ER_model    index number for chosen electron-range model
        (pointer to single variable)
        @param  ER_parameters    parameters for chosen electron-range model
        (pointer to array of size depending on chosen model)
        @see          AT_ElectronRange.h for definition
        @param  gamma_model    index number for chosen gamma response
        (pointer to single variable)
        @param  gamma_parameters  parameters for chosen gamma response
        (pointer to array of size depending on chosen model)
        @see          AT_GammaResponse.h for definition
        @param  n_runs  (algorithm specific) number of points sampled for
        local dose distribution
        @param  N2      (algorithm specific) number of bins per factor of two
        in local dose array (pointer to single variable)
        @param  fluence_factor    factor to scale the fluences given as
        "fluence_cm2" with (pointer to single variable)
        @param  write_output    if true, a protocol is written to "SuccessiveConvolutions.txt"
        in the working directory (pointer to single variable)
        @param  importance_sampling  if unequal zero importance sampling will be
        applied to the single impact local dose distribution
        @param  results      pointer to array of size 10 to be allocated by the
        user which will be used to return the results
        results[0]    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D
        results[1]    d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm
        results[2]    S_HCP           (algorithm independent)  absolute particle response
        results[3]    S_gamma         (algorithm independent)  absolute gamma response
        results[4]    not used        (algorithm independent)
        results[5]    not_used        (algorithm specific)
        results[6]    not used        (algorithm specific)
        results[7]    not_used        (algorithm specific)
        results[8]    not used        (algorithm specific)
        results[9]    not used        (algorithm specific)
        @return  none
        '''
        #conversion of Python variables to C type variables
        n_ctype =                  ctypes.byref(ctypes.c_int(n))

        if n ==1 :
            E_MeV_u_ctype =            ctypes.byref(ctypes.c_float(E_MeV_u[0]))
            particle_no_ctype =        ctypes.byref(ctypes.c_int(particle_no[0]))
            fluence_cm2_ctype =        ctypes.byref(ctypes.c_float(fluence_cm2[0]))           
        else:  # if more than one particle more sophisticated variable construction needed
            print '\nat the moment just single particle fields implemented, stay updated\n'
            sys.exit(0)
#        E_MeV_u_ctype =            ctypes.byref(ctypes.c_float(E_MeV_u))
#        particle_no_ctype =        ctypes.byref(ctypes.c_int(particle_no))
#        fluence_cm2_ctype =        ctypes.byref(ctypes.c_float(fluence_cm2))
        material_no_ctype =        ctypes.byref(ctypes.c_int(material_no))
        RDD_model_ctype =          ctypes.byref(ctypes.c_int(RDD_model))
        RDD_parameters_ctype =     ctypes.byref(ctypes.c_float(RDD_parameters))
        ER_model_ctype =           ctypes.byref(ctypes.c_int(ER_model))
        ER_parameters_ctype =      ctypes.byref(ctypes.c_float(ER_parameters))
        gamma_model_ctype =        ctypes.byref(ctypes.c_int(gamma_model))
        gamma_parameters_ctype =   ctypes.byref(ctypes.c_float(gamma_parameters))
        n_runs_ctype =             ctypes.byref(ctypes.c_long(int(n_runs)))
        N2_ctype =                 ctypes.byref(ctypes.c_int(N2))
        fluence_factor_ctype =     ctypes.byref(ctypes.c_float(fluence_factor))
        write_output_ctype =       ctypes.byref(ctypes.c_int(write_output))
        importance_sampling_ctype = ctypes.byref(ctypes.c_int(importance_sampling))
        tenfloats = ctypes.c_float*10
        results = tenfloats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
        
        c_at_spiss = self.libamtrack.AT_SPISS
        c_at_spiss.restype = ctypes.c_float *10
        
        c_at_spiss_output =  c_at_spiss(n_ctype, E_MeV_u_ctype, particle_no_ctype, fluence_cm2_ctype,
                                        material_no_ctype, RDD_model_ctype, RDD_parameters_ctype,
                                        ER_model_ctype, ER_parameters_ctype, gamma_model_ctype,
                                        gamma_parameters_ctype, n_runs_ctype, N2_ctype, fluence_factor_ctype,
                                        write_output_ctype, importance_sampling_ctype,  results)
        
        AT_SPISS_output = []
        for item in results:
            AT_SPISS_output.append(item)
        print '\nSPISS func return is nonsence, has to be fixwed im AmTrack.c,\nresults in SPISS.log\n'
        return AT_SPISS_output
        
        

    def AT_GSM (self, n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model, RDD_parameters, ER_model,
                ER_parameters, gamma_model, gamma_parameters, N_runs, N2, fluence_factor,
                write_output, nX, voxel_size_m, lethal_events_mode):
        '''
        Computes HCP response and RE/RBE using summation of tracks an a Cartesian grid (the 'GSM' algorithm)

        @param  n      number of particle types in the mixed particle field
        (pointer to single variable)
        @param  E_MeV_u      energy of particles in the mixed particle field
        (pointer to array of size n)
        @param  particle_no    type of the particles in the mixed particle field
        (pointer to array of size n)
        @see          AT_DataParticle.h for definition
        @param  fluence_cm2    fluences for the given particles, doses in Gy if negative
        (pointer to array of size n)
        @param  material_no    index number for detector material
        (pointer to single variable)
        @see          AT_DataMaterial.h for definition
        @param  RDD_model    index number for chosen radial dose distribution
        (pointer to single variable)
        @param  RDD_parameters    parameters for chosen radial dose distribution
        (pointer to array of size depending on chosen model)
        @see          AT_RDD.h for definition
        @param  ER_model    index number for chosen electron-range model
        (pointer to single variable)
        @param  ER_parameters    parameters for chosen electron-range model
        (pointer to array of size depending on chosen model)
        @see          AT_ElectronRange.h for definition
        @param  gamma_model    index number for chosen gamma response
        (pointer to single variable)
        @param  gamma_parameters  parameters for chosen gamma response
        (pointer to array of size depending on chosen model)
        @see          AT_GammaResponse.h for definition
        @param  N_runs       (algorithm specific) number of runs within which track
        positions will be resampled
        @param  N2      (algorithm specific) number of bins per factor of two in
        local dose array (pointer to single variable)
        @param  fluence_factor    factor to scale the fluences given as "fluence_cm2"
        with (pointer to single variable)
        @param  write_output    if true, a protocol is written to "SuccessiveConvolutions.txt"
        in the working directory (pointer to single variable)
        @pram  shrink_tails    (algorithm specific) if true, tails of the local dose
        distribution, contributing less than "shrink_tails_under" are cut (pointer to single variable)
        @param  shrink_tails_under  (algorithm specific) limit for tail cutting in
        local dose distribution (pointer to single variable)
        @param  adjust_N2    (algorithm specific) if true, "N2" will be increase if
        necessary at high fluence to ensure sufficient binning resolution
        @param  lethal_events_mode (algorithm specific) if true, allows to do
        calculations for cell survival
        @param  nX    (algorithm specific) number of voxels of the grid in x
        (and y as the grid is quadratic)
        @param  voxel_size_m     side length of a voxel in m
        @param  results      pointer to array of size 10 to be allocated by the
        user which will be used to return the results
        results[0]    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D
        results[1]    d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm
        results[2]    S_HCP           (algorithm independent)  absolute particle response
        results[3]    S_gamma         (algorithm independent)  absolute gamma response
        results[4]    n_particles     (algorithm independent)  average number of particle tracks on the detector grid
        results[5]    sd_efficiency   (algorithm specific)     standard deviation for results[0]
        results[6]    sd_d_check      (algorithm specific)     standard deviation for results[1]
        results[7]    sd_S_HCP        (algorithm specific)     standard deviation for results[2]
        results[8]    sd_S_gamma      (algorithm specific)     standard deviation for results[3]
        results[9]    sd_n_particles  (algorithm specific)     standard deviation for results[4]
        @return  none
        '''
        #conversion of Python variables to C type variables
        n_ctype =                  ctypes.pointer(ctypes.c_int(n))
        if n ==1 :
            E_MeV_u_ctype =            ctypes.byref(ctypes.c_float(E_MeV_u[0]))
            particle_no_ctype =        ctypes.byref(ctypes.c_int(particle_no[0]))
            fluence_cm2_ctype =        ctypes.byref(ctypes.c_float(fluence_cm2[0]))           
        else:  # if more than one particle more sophisticated variable construction needed
            print '\nat the moment just single particle fields implemented, stay updated\n'
            sys.exit(0)
        material_no_ctype =        ctypes.pointer(ctypes.c_int(material_no))
        RDD_model_ctype =          ctypes.pointer(ctypes.c_int(RDD_model))
        RDD_parameters_ctype =     ctypes.pointer(ctypes.c_float(RDD_parameters))
        ER_model_ctype =           ctypes.pointer(ctypes.c_int(ER_model))
        ER_parameters_ctype =      ctypes.pointer(ctypes.c_float(ER_parameters))
        gamma_model_ctype =        ctypes.pointer(ctypes.c_int(gamma_model))
        gamma_parameters_ctype =   ctypes.pointer(ctypes.c_float(gamma_parameters))
        N_runs_ctype = ctypes.pointer(ctypes.c_long(int(N_runs)))
        N2_ctype =                 ctypes.pointer(ctypes.c_int(N2))
        fluence_factor_ctype =     ctypes.pointer(ctypes.c_float(fluence_factor))
        write_output_ctype =       ctypes.pointer(ctypes.c_bool(write_output)) 
        nX_ctype = ctypes.pointer(ctypes.c_int(nX))
        voxel_size_m_ctype = ctypes.pointer(ctypes.c_float(voxel_size_m))
        lethal_events_mode_ctype = ctypes.pointer(ctypes.c_bool(lethal_events_mode))
        tenfloats = ctypes.c_float*10
        results = tenfloats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
        
        c_at_gsm = self.libamtrack.AT_GSM
        c_at_gsm.restype = ctypes.c_float  *10
        
        c_at_gsm_output = c_at_gsm(n_ctype,  E_MeV_u_ctype,  particle_no_ctype,
                                   fluence_cm2_ctype, material_no_ctype, RDD_model_ctype,
                                   RDD_parameters_ctype, ER_model_ctype, ER_parameters_ctype,
                                   gamma_model_ctype,  gamma_parameters_ctype, N_runs_ctype,
                                   N2_ctype, fluence_factor_ctype, write_output_ctype,  nX_ctype,
                                   voxel_size_m_ctype, lethal_events_mode_ctype,  results)

        AT_GSM_output = []
        for item in results:
            AT_GSM_output.append(item)
            
        return AT_GSM_output
        

    def AT_IGK (self):
        print 'To be included'       
        



def test_suite():
    '''test script for pyamtrack, runs all functions in a row'''
    print'\n----------\n- pyamtrack test run\n-------\n!!!!Python 2.6 has to bee installed to use pyamtrack!!!!!\n'
    # test parameters
    n = 1
    E_MeV_u = [10.0]
    particle_no = [1]
    fluence_cm2 = [-3.]
    material_no =1
    RDD_model = 3
    RDD_parameters =5e-8
    ER_model = 2
    ER_parameters = 0.0
    gamma_model = 4
    gamma_parameters = 1
    N2 = 40
    fluence_factor = 1.0
    write_output = True
    shrink_tails = True
    shrink_tails_under = 1e-30
    adjust_N2 = True
    lethal_events_mode = False
    N_runs= 1e3 
    importance_sampling = 0
    voxel_size_m = 0.001
    nX = 10


    test_obj = AmTrack()
    print test_obj.libamtrack
    print 'SPISS\nResults:\n'
    print test_obj.AT_SPISS(n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model,
                                          RDD_parameters, ER_model, ER_parameters, gamma_model, gamma_parameters, N_runs,  
                                          N2, fluence_factor, write_output, importance_sampling)


    
    print'SPIFF\nResults:\n'
    print test_obj.AT_SPIFF(n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model, RDD_parameters, 
                            ER_model, ER_parameters, gamma_model, gamma_parameters,  N2, fluence_factor, write_output,  
                            shrink_tails, shrink_tails_under, adjust_N2, lethal_events_mode)

    print'GSM\ntest skipped\n'
   # test_out= test_obj.AT_GSM(n, E_MeV_u, particle_no, fluence_cm2, material_no,
   #                          RDD_model, RDD_parameters, ER_model, ER_parameters,
   #                          gamma_model, gamma_parameters, N_runs, N2, fluence_factor,  write_output, nX, voxel_size_m, lethal_events_mode)


if __name__ == "__main__":
        test_suite()

    
