#! /usr/bin/env python

'''
   pyam_tools.py

   =========

   Created on: 18.03.2011
   Creator: herrmann
   
   Python module holding tools to ease the use of pyamtrack,
   this file has to be adapted manually
   
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

import pyamtrack


class PyAmRun(object):
    '''
    \'Mother\'-class for AmTrack objects
    '''
    def __init__(self, name = 'amtackrun'):
        '''
        PuyAmRun object
        run object the under self.algortihm specified algorithm with the given parameters
        default \'CPPSC\'
        '''
        self.name = name
        # general parameter
        self.number_of_components = int
        self.E_MeV_u = list
        self.particle_no = list
        self.fluence_cm2 = list
        self.material_no = int
        self.rdd_model =int
        self.rdd_parameters = [0.0]*4
        self.er_model = int
        self.gamma_model = int
        self.gamma_parameters = [0.0]*9
        self.lethal_events_mode = False
        self.write_output = True
        # default algorithm CPPSC
        self.algorithm = 'CPPSC'
        
        # CPPSC specific input parameter 
        self.N2 = 120
        self.fluence_factor = 1.0
        self.shrink_tails = True
        self.shrink_tails_under = 1e-30
        self.adjust_N2 =True

        # IGK specific input parameter
        self.saturation_cross_section_factor= 1.0

        # GSM specific input parameter
        self.N_runs = 100
        self.nX = 100
        self.voxel_size_m = 1e-4
        
        #output parameter
        #general
        self.relative_efficiency = 0.0
        self.S_HCP = 0.0
        self.S_gamma = 0.0
        #CPPSC
        self.d_check = 0.0
        self.mean_number_of_tracks_contrib = 0.0
        self.start_number_of_tracks_contrib = 0.0
        self.n_convolutions = 0
        self.lower_Jensen_bound = 0.0
        self.upper_Jensen_bound = 0.0
        #IGK
        self.sI_cm2 =0.0
        self.gamma_dose_Gy =0.0
        self.P_I = 0.0
        self.P_g = 0.0
        #GSM
        self.n_particles = 0.0 
        self.sd_relative_efficiency = 0.0
        self.sd_d_check = 0.0 
        self.sd_S_HCP = 0.0 
        self.sd_S_gamma = 0.0
        self.sd_n_particles = 0.0

    def run(self):
        '''
        run specified algorithm
        '''
        if self.algorithm == 'CPPSC':
            self.runCPPSC()
        if self.algorithm == 'IGK':
            self.runIGK()
        if self.algorithm == 'GSM':
            self.runGSM()
            

    def runCPPSC(self):
        '''
        run CPPSC algorithm
        '''
        results = pyamtrack.AT_run_CPPSC_method (self.number_of_components,
                                                 self.E_MeV_u,
                                                 self.particle_no,
                                                 self.fluence_cm2,
                                                 self.material_no,
                                                 self.rdd_model,
                                                 self.rdd_parameters,
                                                 self.er_model,
                                                 self.gamma_model,
                                                 self.gamma_parameters,
                                                 self.N2,
                                                 self.fluence_factor,
                                                 self.write_output,
                                                 self.shrink_tails,
                                                 self.shrink_tails_under,
                                                 self.adjust_N2,
                                                 self.lethal_events_mode,
                                                 self.relative_efficiency,
                                                 self.d_check,
                                                 self.S_HCP,
                                                 self.S_gamma,
                                                 self.mean_number_of_tracks_contrib,
                                                 self.start_number_of_tracks_contrib,
                                                 self.n_convolutions,
                                                 self.lower_Jensen_bound,
                                                 self.upper_Jensen_bound)
        self.relative_efficiency = results[0]
        self.d_check = results[1]
        self.S_HCP = results[2]
        self.S_gamma = results[3]
        self.mean_number_of_tracks_contrib = results[4]
        self.start_number_of_tracks_contrib = results[5]
        self.n_convolutions = results[6]
        self.lower_Jensen_bound = results[7]
        self.upper_Jensen_bound  = results[8]

        
        
    def runIGK(self):
        '''
        run IGK algorithm
        '''
        results = pyamtrack.AT_run_IGK_method ( self.number_of_components,
                                                self.E_MeV_u,
                                                self.particle_no,
                                                self.fluence_cm2,
                                                self.material_no,
                                                self.rdd_model,
                                                self.rdd_parameters,
                                                self.er_model,
                                                self.gamma_model,
                                                self.gamma_parameters,
                                                self.saturation_cross_section_factor,
                                                self.write_output,
                                                self.relative_efficiency,
                                                self.S_HCP,
                                                self.S_gamma,
                                                self.sI_cm2,
                                                self.gamma_dose_Gy,
                                                self.P_I,
                                                self.P_g)
        self.relative_efficiency= results[0]
        self.S_HCP = results[1]
        self.S_gamma = results[2]
        self.sI_cm2 = results[3]
        self.gamma_dose_Gy = results[4]
        self.P_I = results[5]
        self.P_g = results[6]

    def runGSM(self):
        '''
        run GSM algorithm
        '''
        results = pyamtrack.AT_run_GSM_method ( self.number_of_components,
                                                self.E_MeV_u,
                                                self.particle_no,
                                                self.fluence_cm2,
                                                self.material_no,
                                                self.rdd_model,
                                                self.rdd_parameters,
                                                self.er_model,
                                                self.gamma_model,
                                                self.gamma_parameters,
                                                self.N_runs,
                                                self.write_output,
                                                self.nX,
                                                self.voxel_size_m,
                                                self.lethal_events_mode,
                                                self.relative_efficiency,
                                                self.d_check,
                                                self.S_HCP,
                                                self.S_gamma,
                                                self.n_particles,
                                                self.sd_relative_efficiency,
                                                self.sd_d_check,
                                                self.sd_S_HCP,
                                                self.sd_S_gamma,
                                                self.sd_n_particles)
        self.relative_efficiency = results[0]
        self.d_check = results[1]
        self.S_HCP = results[2]
        self.S_gamma = results[3]
        self.n_particles = results[4]
        self.sd_relative_efficiency = results[5]
        self.sd_d_check = results[6]
        self.sd_S_HCP = results[7]
        self.sd_S_gamma = results[8]
        self.sd_n_particles = results[9]
