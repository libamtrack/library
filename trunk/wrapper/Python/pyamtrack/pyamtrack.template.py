#! /usr/bin/env python

'''
   pyamtrack.py

   =========

   Created on: 16.02.2010
   Creator: herrmann
   
   Python class for interfacing AmTrack functions. 
   
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

import ctypes
import sys
from  platform import python_version
from platform import system as os_system
try:
    import numpy
except:
    print 'numpy library not found'


__status__ = 'Prototype'

operating_system = os_system()
if operating_system == 'Windows':
    libamtrack = ctypes.cdll.libamtrack
else:
    libamtrack = ctypes.cdll.LoadLibrary("libamtrack.so")

# python version controll
py_version = python_version()
if int(py_version[0]) < 3 and int(py_version[2]) <= 5:
    print 'ERROR:\nYou are using Python%s\nPython2.6 or later required!\n'%py_version
    sys.exit(0)

    
class AmTrack(object):
    '''
    \'Mother\' class for the AmTrack objects.
    '''
    def __init__ (self):
        print '\npyamtrack\n-----\n'

    def parameters(self):
        '''
        prints all parameters of the object
        '''
        for item in self.__dict__:
            print item, self.__dict__[item]
    
