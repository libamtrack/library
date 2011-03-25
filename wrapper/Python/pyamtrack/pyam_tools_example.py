#! /usr/bin/env python
'''
    pyam_tools_example.py
    =========
   Created on: 09.03.2011
   Author: herrmann
   
   Example script for using the pyam_tools.py module for libamtrack
   
   Requirements: Python2.6

   Copyright 2006, 2011 Steffen Greilich / the libamtrack team

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


import pyam_tools as pyam

# create PyAmRun object, called TEST
test = pyam.PyAmRun('TEST')

# define Attributs of the object
test.algorithm = 'IGK'
test.number_of_components = 1 
test.E_MeV_u = [10.0]   # 10 MeV
test.particle_no = [1001] # protons
test.fluence_cm2 = [-10.0] # 10 Gy
test.material_no = 5 # Alanine
test.rdd_model = 4
test.rdd_parameters[:2]= [1e-10,1e-10]
test.er_model =2
test.gamma_model = 2
test.gamma_parameters[:4]=[1,10.5e4,1,1]


# Run algorithm on object
test.run() # note test.runIGK() would call the same function
igk_re = test.relative_efficiency

# change algorithm
test.algorithm= 'CPPSC'
test.run() # note test.runCPPSC() would call the same function
cppsc_re = test.relative_efficiency

# print the results
print '### TEST ###'
print 'RE\nIGK :\t%.3f\nCPPSC :\t%.3f'%(igk_re, cppsc_re)

