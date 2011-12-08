#! /usr/bin/env python
'''
    pyam_obj.py
    =========
   Created on: 09.03.2011
   Author: herrmann
   
   Python module holding pyamtrack classes for interface with libamtrack 
   
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

import string
import re

class amtrack_func:
    def __init__(self, func_name):
        '''
        @param parameter : list of amtrack_para objects
        '''
        self.name = func_name
        self.parameter = dict()
        self.comment = ''
        self.definition =''
        self.c_return = False
        
    def extract_params (self):
        '''
        analyzes the comment and searches for the input and output parameter,
        includes found parameters into self.parameter, also characteristics as pointer, array etc.
        '''
        no =0
        for lineno,  line in enumerate(self.comment):
            line_list = string.split(line)
            if len(line_list) > 1:
                if line_list[1][:6] == '@param':
                    no +=1
                    try:
                        self.parameter[line_list[2]]=amtrack_para(line_list[2], no)
                        self.parameter[line_list[2]].direction=string.split(line_list[1], '[')[1][:-1]
                        is_array = re.search('array of size', line)
                        if is_array is not None:
                            self.parameter[line_list[2]].array = True
                            self.parameter[line_list[2]].pointer = True
                            right = string.split(line,  '(')[-1]
                            middle = string.split(right, ')')[0]
                            length = string.split(middle)[-1]
                            try:
                                self.parameter[line_list[2]].array_size = lnt(length)
                            except:
                                self.parameter[line_list[2]].array_size = length
                        self.parameter[line_list[2]].works = True
                    except:
                        print self.name + ' problems occurred'
                if line_list[1] == '@return':
                    self.c_return = True
                    
        for item in self.definition:
            oneline = False
            item =string.split(item, ',')[0]
            if len(string.split(item, '(')) == 2:
                tmp_string_ls = string.split(item, '(')[-1]
                if len(string.split(item, ')')) == 2:
                    oneline =True
                if oneline :
                    tmp_name = string.split(tmp_string_ls)[-2]
                else:
                    tmp_name = string.split(tmp_string_ls)[-1]
                if tmp_name [-2:] == '[]':
                    tmp_name = tmp_name[:-2]
                    self.parameter[tmp_name].array = True
                tmp_type = string.split(tmp_string_ls)[-2]
                if oneline :
                    tmp_type = string.split(tmp_string_ls)[-3]
                if tmp_type[-1] == '*':
                    tmp_type = tmp_type[:-1]
                    self.parameter[tmp_name].pointer =True
                self.parameter[tmp_name].type=tmp_type
            elif len(string.split(item, ')')) == 2  and not oneline:
                tmp_string_ls = string.split(item, ')')
                tmp_name = string.split(tmp_string_ls[0])[-1]
                if tmp_name [-2:] == '[]':
                    tmp_name = tmp_name[:-2]
                    self.parameter[tmp_name].array = True
                tmp_type = string.split(tmp_string_ls[0])[-2] 
                if tmp_type[-1] == '*':
                    tmp_type = tmp_type[:-1]
                    self.parameter[tmp_name].pointer =True
                self.parameter[tmp_name].type=tmp_type
            else:
                tmp_name = string.split(item)[-1]
                if tmp_name [-2:] == '[]':
                    tmp_name = tmp_name[:-2]
                    self.parameter[tmp_name].array = True  
                tmp_type = string.split(item)[-2]
                if tmp_type[-1] == '*':
                    tmp_type = tmp_type[:-1]
                    self.parameter[tmp_name].pointer =True
                self.parameter[tmp_name].type=tmp_type
                    
class amtrack_para:
    def __init__(self,  name,  no , type='void',  direction='in'):
        '''
        amtrack parameter class
        
        '''
        self.name = name
        self.type = type
        self.direction = 'input'
        self.array = False
        self.array_size = 1
        self.pointer = False
        self.no = no
        self.works = False# if yes function had no problems


if __name__ == "__main__":
    print 'pyam_obj module, providing classes for pyamtrack wrapper'


