#! /usr/bin/env python
'''
    generate_Py_wrapper.py
    =========
   Created on: 09.03.2011
   Author: herrmann
   
   Python script to generate wrapper for libamtrack 
   
   desired functions have to be given in file NAMESPACE
   
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

import py_parser

def main():
    
    print 'Reading NAMESPACE file'
    desired_func = py_parser.read_namespace()
    print 'Harvesting function definitions in header files'
    c_functions = py_parser.harvest(desired_func)
    print 'Writing \'pyamtrack.py\' wrapper file'
    py_parser.write_func_in_py(c_functions, 'pyamtrack.py')
    print 'Done'
    
    
main()
