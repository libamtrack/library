#! /usr/bin/env python

'''
   pyamtrack.py

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

#TODO: SPISS to be made consistent in main library
#TODO: complete epydoc markups


import ctypes
import sys
import platform
import numpy

__status__ = 'Prototype'


class AmTrack(object):
    def __init__ (self):
        print '\npyamtrack\n-----\n'
        operating_system = platform.system()
        if operating_system == 'Windows':
            self.libamtrack = ctypes.cdll.libamtrack
            # should be working under Windows, but never has been tested
            # if you have tried it under Windows, I'D be happy to hear your
            # experience
        else:
            self.libamtrack = ctypes.cdll.LoadLibrary("libamtrack.so")
        # python version controll
        py_version = platform.python_version()
        if int(py_version[0]) < 3 and int(py_version[2]) <= 5:
            print 'ERROR:\nYou are using Python%s\nPython2.6 or later needed!\n'%py_version
            sys.exit(0)
        

    def AT_SPIFF (self, n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model, RDD_parameters, 
                  ER_model, ER_parameters, gamma_model, gamma_parameters,  N2, fluence_factor, write_output,  
                  shrink_tails, shrink_tails_under, adjust_N2, lethal_events_mode):
        ''' 
        Computes HCP response and RE/RBE using compound Poison process and successive convolutions (CPP_SC, the 'SPIFF' algorithm)

        @param  n:     number of particle types in the mixed particle field 
        @type n: integer
        @param  E_MeV_u:      energy of particles in the mixed particle field
        @type   E_MeV_u: float list of length 'n'
        @param  particle_no:    type of the particles in the mixed particle field 
        @see       AT_DataParticle.h: for definition
        @type particle_no: integer list of length 'n'
        @param  fluence_cm2:    fluences for the given particles, doses in Gy if negative
        @type  fluence_cm2:  float list of length 'n'
        @param  material_no:    index number for detector material
        @see          AT_DataMaterial.h for definition
        @type material_no: integer
        @param  RDD_model:    index number for chosen radial dose distribution
        @type RDD_model: integer
        @param  RDD_parameters:   parameters for chosen radial dose distribution
        @see          AT_RDD.h for definition
        @type RDD_parameters: list
        @param  ER_model:   index number for chosen electron-range model
        @type ER_model: integer
        @param  ER_parameters:   parameters for chosen electron-range model
        @see          AT_ElectronRange.h for definition
        @type ER_parameters: array of model depending length
        @param  gamma_model:   index number for chosen gamma response
        @type gamma_model: integer
        @param  gamma_parameters: parameters for chosen gamma response
        @see          AT_GammaResponse.h for definition
        @type gamma_parameters: array of model depending length
        @param  N2:     (algorithm specific) number of bins per factor of two in local dose array
        @type N2 : number
        @param  fluence_factor:   factor to scale the fluences given as "fluence_cm2" with
        @type fluence_factor: integer
        @param  write_output:   if true, a protocol is written to "SuccessiveConvolutions.txt"
        in the working directory
        @type write_output: boolean
        @param  shrink_tails:   (algorithm specific) if true, tails of the local dose
        distribution, contributing less than "shrink_tails_under" are cut
        @type shrink_tails: boolean
        @param  shrink_tails_under: (algorithm specific) limit for tail cutting
        in local dose distribution
        @type shrink_tails_under: float
        @param  adjust_N2:   (algorithm specific) if true, "N2" will be increase
        if necessary at high fluence to ensure sufficient binning resolution
        @type adjust_N2: boolean
        @param  lethal_events_mode:(algorithm specific) if true, allows to do
        calculations for cell survival
        @type lethal_events_mode: boolean
        @param  results:     pointer to array of size 10 to be allocated by the
        user which will be used to return the results
        @param results[0]:   efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D
        @param results[1]:   d_check        (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm
        @param results[2]:   S_HCP           (algorithm independent)  absolute particle response
        @param results[3]:   S_gamma         (algorithm independent)  absolute gamma response
        @param results[4]:   not used        (algorithm independent)
        @param results[5]:   u               (algorithm specific)     mean number of tracks contributing to representative point
        @param results[6]:   u_start         (algorithm specific)     low starting value for mean number of tracks, where linearisation is applied
        @param results[7]:   n_convolutions  (algorithm specific)     number of convolutions performed 
        @param results[8]:   not used        (algorithm specific)
        @param results[9]:   not used        (algorithm specific)
        @return:  None
        '''
        #conversion of Python variables to C type variables
        n_ctype =                  ctypes.byref(ctypes.c_int(n))
        tmp_array =         ctypes.c_float * len(E_MeV_u)
        E_MeV_u_ctype =     ctypes.byref(tmp_array(*E_MeV_u))
        tmp_array =         ctypes.c_long * len(particle_no)
        particle_no_ctype = ctypes.byref(tmp_array(*particle_no))
        tmp_array =         ctypes.c_float * len(fluence_cm2)
        fluence_cm2_ctype =  ctypes.byref(tmp_array(*fluence_cm2))
        material_no_ctype =        ctypes.byref(ctypes.c_int(material_no))
        RDD_model_ctype =          ctypes.byref(ctypes.c_int(RDD_model))
        tmp_array = ctypes.c_float* len(RDD_parameters)        
        RDD_parameters_ctype =     ctypes.byref(tmp_array(*RDD_parameters))   
        ER_model_ctype =           ctypes.byref(ctypes.c_int(ER_model))        
        tmp_array = ctypes.c_float* len(ER_parameters)        
        ER_parameters_ctype =      ctypes.byref(tmp_array(*ER_parameters))        
        gamma_model_ctype =        ctypes.byref(ctypes.c_int(gamma_model))        
        tmp_array = ctypes.c_float* len(gamma_parameters)        
        gamma_parameters_ctype =   ctypes.byref(tmp_array(*gamma_parameters))    
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
        c_at_spiff_output = c_at_spiff(n_ctype, E_MeV_u_ctype, particle_no_ctype, fluence_cm2_ctype, material_no_ctype,
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

        @param  n:     number of particle types in the mixed particle field 
        @type n: integer
        @param  E_MeV_u:      energy of particles in the mixed particle field
        @type   E_MeV_u: float list of length 'n'
        @param  particle_no:    type of the particles in the mixed particle field 
        @see       AT_DataParticle.h: for definition
        @type particle_no: integer list of length 'n'
        @param  fluence_cm2:    fluences for the given particles, doses in Gy if negative
        @type  fluence_cm2:  float list of length 'n'
        @param  material_no:    index number for detector material
        @see          AT_DataMaterial.h for definition
        @type material_no: integer
        @param  RDD_model:    index number for chosen radial dose distribution
        @type RDD_model: integer
        @param  RDD_parameters:   parameters for chosen radial dose distribution
        @see          AT_RDD.h for definition
        @type RDD_parameters: list
        @param  ER_model:   index number for chosen electron-range model
        @type ER_model: integer
        @param  ER_parameters:   parameters for chosen electron-range model
        @see          AT_ElectronRange.h for definition
        @type ER_parameters: array of model depending length
        @param  gamma_model:   index number for chosen gamma response
        @type gamma_model: integer
        @param  gamma_parameters: parameters for chosen gamma response
        @see          AT_GammaResponse.h for definition
        @type gamma_parameters: array of model depending length
        @param  n_runs: (algorithm specific) number of points sampled for
        local dose distribution
        @param  N2:     (algorithm specific) number of bins per factor of two
        in local dose array (pointer to single variable)
        @param  fluence_factor:   factor to scale the fluences given as
        "fluence_cm2" with (pointer to single variable)
        @param  write_output:   if true, a protocol is written to "SuccessiveConvolutions.txt"
        in the working directory (pointer to single variable)
        @param  importance_sampling: if unequal zero importance sampling will be
        applied to the single impact local dose distribution
        @param  results:     pointer to array of size 10 to be allocated by the
        user which will be used to return the results
        @param results[0]:   efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D
        @param results[1]:   d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm
        @param results[2]:   S_HCP           (algorithm independent)  absolute particle response
        @param results[3]:   S_gamma         (algorithm independent)  absolute gamma response
        @param results[4]:   not used        (algorithm independent)
        @param results[5]:   not_used        (algorithm specific)
        @param results[6]:   not used        (algorithm specific)
        @param results[7]:   not_used        (algorithm specific)
        @param results[8]:   not used        (algorithm specific)
        @param results[9];   not used        (algorithm specific)
        @return: None
        '''
        #conversion of Python variables to C type variables
        n_ctype =                  ctypes.byref(ctypes.c_int(n))
        tmp_array =         ctypes.c_float * len(E_MeV_u)
        E_MeV_u_ctype =     ctypes.byref(tmp_array(*E_MeV_u))
        tmp_array =         ctypes.c_long * len(particle_no)
        particle_no_ctype = ctypes.byref(tmp_array(*particle_no))
        tmp_array =         ctypes.c_float * len(fluence_cm2)
        fluence_cm2_ctype =  ctypes.byref(tmp_array(*fluence_cm2))
        material_no_ctype =        ctypes.byref(ctypes.c_int(material_no))
        RDD_model_ctype =          ctypes.byref(ctypes.c_int(RDD_model))
        tmp_array = ctypes.c_float* len(RDD_parameters)        
        RDD_parameters_ctype =     ctypes.byref(tmp_array(*RDD_parameters))   
        ER_model_ctype =           ctypes.byref(ctypes.c_int(ER_model))        
        tmp_array = ctypes.c_float* len(ER_parameters)        
        ER_parameters_ctype =      ctypes.byref(tmp_array(*ER_parameters))        
        gamma_model_ctype =        ctypes.byref(ctypes.c_int(gamma_model))        
        tmp_array = ctypes.c_float* len(gamma_parameters)        
        gamma_parameters_ctype =   ctypes.byref(tmp_array(*gamma_parameters))
        N_runs_ctype =             ctypes.byref(ctypes.c_long(int(n_runs)))
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
                                        gamma_parameters_ctype, N_runs_ctype, N2_ctype, fluence_factor_ctype,
                                        write_output_ctype, importance_sampling_ctype,  results)
        
        AT_SPISS_output = []
        for item in results:
            AT_SPISS_output.append(item)
        print '\nSPISS func return is nonsense, has to be fixed im AmTrack.c,\nresults in SPISS.log\n'
        return AT_SPISS_output
        
        

    def AT_GSM (self, n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model, RDD_parameters, ER_model,
                ER_parameters, gamma_model, gamma_parameters, N_runs, N2, fluence_factor,
                write_output, nX, voxel_size_m, lethal_events_mode):
        '''
        Computes HCP response and RE/RBE using summation of tracks an a Cartesian grid (the 'GSM' algorithm)

        @param  n:     number of particle types in the mixed particle field 
        @type n: integer
        @param  E_MeV_u:      energy of particles in the mixed particle field
        @type   E_MeV_u: float list of length 'n'
        @param  particle_no:    type of the particles in the mixed particle field 
        @see       AT_DataParticle.h: for definition
        @type particle_no: integer list of length 'n'
        @param  fluence_cm2:    fluences for the given particles, doses in Gy if negative
        @type  fluence_cm2:  float list of length 'n'
        @param  material_no:    index number for detector material
        @see          AT_DataMaterial.h for definition
        @type material_no: integer
        @param  RDD_model:    index number for chosen radial dose distribution
        @type RDD_model: integer
        @param  RDD_parameters:   parameters for chosen radial dose distribution
        @see          AT_RDD.h for definition
        @type RDD_parameters: list
        @param  ER_model:   index number for chosen electron-range model
        @type ER_model: integer
        @param  ER_parameters:   parameters for chosen electron-range model
        @see          AT_ElectronRange.h for definition
        @type ER_parameters: array of model depending length
        @param  gamma_model:   index number for chosen gamma response
        @type gamma_model: integer
        @param  gamma_parameters: parameters for chosen gamma response
        @see          AT_GammaResponse.h for definition
        @type gamma_parameters: array of model depending length
        @param  N_runs:      (algorithm specific) number of runs within which track
        positions will be resampled
        @param  N2:     (algorithm specific) number of bins per factor of two in
        local dose array (pointer to single variable)
        @param  fluence_factor:   factor to scale the fluences given as "fluence_cm2"
        with (pointer to single variable)
        @param  write_output:   if true, a protocol is written to "SuccessiveConvolutions.txt"
        in the working directory (pointer to single variable)
        @param  shrink_tails:   (algorithm specific) if true, tails of the local dose
        distribution, contributing less than "shrink_tails_under" are cut (pointer to single variable)
        @param  shrink_tails_under: (algorithm specific) limit for tail cutting in
        local dose distribution (pointer to single variable)
        @param  adjust_N2:   (algorithm specific) if true, "N2" will be increase if
        necessary at high fluence to ensure sufficient binning resolution
        @param  lethal_events_mode:(algorithm specific) if true, allows to do
        calculations for cell survival
        @param  nX:   (algorithm specific) number of voxels of the grid in x
        (and y as the grid is quadratic)
        @param  voxel_size_m:    side length of a voxel in m
        @param  results:     pointer to array of size 10 to be allocated by the
        user which will be used to return the results
        @param results[0]:   efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D
        @param results[1]:   d_check         (algorithm independent)  sanity check:  total dose (in Gy) as returned by the algorithm
        @param results[2]:   S_HCP           (algorithm independent)  absolute particle response
        @param results[3]:   S_gamma         (algorithm independent)  absolute gamma response
        @param results[4]:   n_particles     (algorithm independent)  average number of particle tracks on the detector grid
        @param results[5]:   sd_efficiency   (algorithm specific)     standard deviation for results[0]
        @param results[6]:   sd_d_check      (algorithm specific)     standard deviation for results[1]
        @param results[7]:   sd_S_HCP        (algorithm specific)     standard deviation for results[2]
        @param results[8]:   sd_S_gamma      (algorithm specific)     standard deviation for results[3]
        @param results[9]:   sd_n_particles  (algorithm specific)     standard deviation for results[4]
        @return: None
        '''
        #conversion of Python variables to C type variables
        n_ctype =                  ctypes.pointer(ctypes.c_int(n))
        tmp_array =         ctypes.c_float * len(E_MeV_u)
        E_MeV_u_ctype =     ctypes.byref(tmp_array(*E_MeV_u))
        tmp_array =         ctypes.c_long * len(particle_no)
        particle_no_ctype = ctypes.byref(tmp_array(*particle_no))
        tmp_array =         ctypes.c_float * len(fluence_cm2)
        fluence_cm2_ctype =  ctypes.byref(tmp_array(*fluence_cm2))
        material_no_ctype =        ctypes.byref(ctypes.c_int(material_no))
        RDD_model_ctype =          ctypes.byref(ctypes.c_int(RDD_model))
        tmp_array = ctypes.c_float* len(RDD_parameters)        
        RDD_parameters_ctype =     ctypes.byref(tmp_array(*RDD_parameters))   
        ER_model_ctype =           ctypes.byref(ctypes.c_int(ER_model))        
        tmp_array = ctypes.c_float* len(ER_parameters)        
        ER_parameters_ctype =      ctypes.byref(tmp_array(*ER_parameters))        
        gamma_model_ctype =        ctypes.byref(ctypes.c_int(gamma_model))        
        tmp_array = ctypes.c_float* len(gamma_parameters)        
        gamma_parameters_ctype =   ctypes.byref(tmp_array(*gamma_parameters))        
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
        

    def AT_IGK (self, n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model, RDD_parameters, ER_model,
                ER_parameters, gamma_model, gamma_parameters, saturation_cross_section_factor):
        '''
        Computes HCP response and RE/RBE using Katz\' Ion-Gamma-Kill approach
        according to Waligorski, 1988

        @param  n:     number of particle types in the mixed particle field 
        @type n: integer
        @param  E_MeV_u:      energy of particles in the mixed particle field
        @type   E_MeV_u: float list of length 'n'
        @param  particle_no:    type of the particles in the mixed particle field 
        @see       AT_DataParticle.h: for definition
        @type particle_no: integer list of length 'n'
        @param  fluence_cm2:    fluences for the given particles, doses in Gy if negative
        @type  fluence_cm2:  float list of length 'n'
        @param  material_no:    index number for detector material
        @see          AT_DataMaterial.h for definition
        @type material_no: integer
        @param  RDD_model:    index number for chosen radial dose distribution
        @type RDD_model: integer
        @param  RDD_parameters:   parameters for chosen radial dose distribution
        @see          AT_RDD.h for definition
        @type RDD_parameters: list
        @param  ER_model:   index number for chosen electron-range model
        @type ER_model: integer
        @param  ER_parameters:   parameters for chosen electron-range model
        @see          AT_ElectronRange.h for definition
        @type ER_parameters: array of model depending length
        @param  gamma_model:   index number for chosen gamma response
        @type gamma_model: integer
        @param  gamma_parameters: parameters for chosen gamma response
        @see          AT_GammaResponse.h for definition
        @type gamma_parameters: array of model depending length 
        @param  saturation_cross_section_factor:  (algorithm specific)  scaling factor for the saturation cross section
        @see          Waligorski, 1988
        @type saturation_cross_section_factor: float
        @param  results:      pointer to array of size 10 to be allocated by the user which will be used to return the results
        @param results[0]:    efficiency      (algorithm independent)  main result:   particle response at dose D / gamma response at dose D
        @param results[1]:    d_check         (algorithm independent)  not available with IGK
        @param results[2]:    S_HCP           (algorithm independent)  absolute particle response
        @param results[3]:    S_gamma         (algorithm independent)  absolute gamma response
        @param results[4]:    n_particles     (algorithm independent)  not available with IGK
        @param results[5]:    sI_cm2          (algorithm specific)     resulting ion saturation cross section in cm2
        @param results[6]:    gamma_dose_Gy   (algorithm specific)     dose contribution from gamma kills
        @param results[7]:    P_I             (algorithm specific)     ion kill probability
        @param results[8]:    P_G             (algorithm specific)     gamma kill probability
        @param results[9]:    not used        (algorithm specific)
        @return:  none
        '''
        #conversion of Python variables to C type variables
        n_ctype =                  ctypes.byref(ctypes.c_int(n))
        tmp_array =         ctypes.c_float * len(E_MeV_u)
        E_MeV_u_ctype =     ctypes.byref(tmp_array(*E_MeV_u))
        tmp_array =         ctypes.c_long * len(particle_no)
        particle_no_ctype = ctypes.byref(tmp_array(*particle_no))
        tmp_array =         ctypes.c_float * len(fluence_cm2)
        fluence_cm2_ctype =  ctypes.byref(tmp_array(*fluence_cm2))
        material_no_ctype =        ctypes.byref(ctypes.c_int(material_no))
        RDD_model_ctype =          ctypes.byref(ctypes.c_int(RDD_model))
        tmp_array = ctypes.c_float* len(RDD_parameters)        
        RDD_parameters_ctype =     ctypes.byref(tmp_array(*RDD_parameters))   
        ER_model_ctype =           ctypes.byref(ctypes.c_int(ER_model))        
        tmp_array = ctypes.c_float* len(ER_parameters)        
        ER_parameters_ctype =      ctypes.byref(tmp_array(*ER_parameters))        
        gamma_model_ctype =        ctypes.byref(ctypes.c_int(gamma_model))        
        tmp_array = ctypes.c_float* len(gamma_parameters)        
        gamma_parameters_ctype =   ctypes.byref(tmp_array(*gamma_parameters))
        saturation_cross_section_factor_ctype = ctypes.c_float(saturation_cross_section_factor)      
        tenfloats = ctypes.c_float*10
        results = tenfloats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
        
        c_at_igk = self.libamtrack.AT_IGK
        c_at_igk.restype = ctypes.c_float  *10       
        
        c_at_igk_output = c_at_igk(n_ctype,  E_MeV_u_ctype,  particle_no_ctype,
                                   fluence_cm2_ctype, material_no_ctype, RDD_model_ctype,
                                   RDD_parameters_ctype, ER_model_ctype, ER_parameters_ctype,
                                   gamma_model_ctype,  gamma_parameters_ctype, saturation_cross_section_factor_ctype,  results)
                                   
        AT_IGK_output = []
        for item in results:
            AT_IGK_output.append(item)
            
        return AT_IGK_output
   


    def AT_interparticleDistance_m(self):
        print 'To be included'

    def AT_interparticleDistance_Gy(self):
        print 'To be included'

    def AT_interparticleDistance_cm2(self):
        print 'To be included'                        




def test_suite():
    '''
    test script for pyamtrack, runs all functions in a row
    '''
    print'\n----------\n- pyamtrack test run\n----------\n\n'
    # test parameters
    #variable set
    material_no = 51# Alanine
    n = 2
    E_MeV_u = [100.0, 100.0]
    particle_no = [18,1] #C-12
    fluence_cm2 = [-10. ,-1.0]
    RDD_model = 3 # Geiss
    RDD_parameters =  [5e-8, 0]
    ER_model = 4 # Geiss
    ER_parameters = [1]
    gamma_model = 4
    gamma_parameters = [1, 10]
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
    saturation_cross_section_factor = 1.0

    test_obj = AmTrack()
    print 'SPISS\nResults:\n'
    print test_obj.AT_SPISS(n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model,
                                          RDD_parameters, ER_model, ER_parameters, gamma_model, gamma_parameters, N_runs,  
                                          N2, fluence_factor, write_output, importance_sampling)


    print'SPIFF\nResults:\n'
    print test_obj.AT_SPIFF(n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model, RDD_parameters, 
                            ER_model, ER_parameters, gamma_model, gamma_parameters,  N2, fluence_factor, write_output,  
                            shrink_tails, shrink_tails_under, adjust_N2, lethal_events_mode)

    
    print'IGK\nResults\n'
    print test_obj.AT_IGK (self, n, E_MeV_u, particle_no, fluence_cm2, material_no, RDD_model, RDD_parameters, ER_model,
                ER_parameters, gamma_model, gamma_parameters, saturation_cross_section_factor)

    print'GSM\ntest skipped\n'
    # test_out= test_obj.AT_GSM(n, E_MeV_u, particle_no, fluence_cm2, material_no,
    #                          RDD_model, RDD_parameters, ER_model, ER_parameters,
    #                          gamma_model, gamma_parameters, N_runs, N2, fluence_factor,  write_output, nX, voxel_size_m, lethal_events_mode)


if __name__ == "__main__":
        test_suite()

    
