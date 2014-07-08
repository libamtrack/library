#! /usr/bin/env python
'''
   pyamtrack_example.py
   ===========
   Created on: 11.03.2011
   Author: herrmann
   
   Example file for use of pyamtrack class.
   
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

import pyamtrack


def main():
    #example use of libamtrack function \'  AT_average_A_from_composition\'
    # first have a look at the description
    print pyamtrack.AT_average_A_from_composition.__doc__
    status_code, average_A =  pyamtrack.AT_average_A_from_composition(2, [1,16], [2./18.,16./18.], 0.0)
    
    print average_A

main()
