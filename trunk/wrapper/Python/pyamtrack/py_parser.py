#! /usr/bin/env python
'''
    py_wrapper.py
    =========
   Created on: 09.03.2011
   Author: herrmann
   
   Python class for interfacing AmTrack functions. 
   
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
import os

import pyam_obj as tools


DEBUG = True
DEBUG =False


def read_namespace(name_file ='NAMESPACE'):
    '''
    reads llist of functions from NAMESPACE file and 
    returns list of functions, which then should be wrapped to python
    '''
    functions = []
    in_file = open(name_file, 'r')
    content = in_file.readlines()
    for line in content:
        if line[0] != '#':
            functions.append(string.strip(line))
    in_file.close()
    return functions

def write_func_in_py(func_objects,  outfile_name = 'pyamtrack.py'):
    '''
    @param[in] func_objects  list of function objects 
    @return : string containing python code to call function
    '''
    header_string = '#! /usr/bin/env python\n#here will be a header \n\n'
    header_string +='import ctypes\n'
    header_string +='import string\n'
    header_string +='import sys\nfrom platform import python_version\nfrom platform import system as os_system\n'
    header_string +='\n__status__ = \'Prototype\'\n\n'
    header_string +='operating_system = os_system()\n'
    header_string +='if operating_system == \'Windows\':\n\tlibamtrack = ctypes.cdll.libamtrack\nelse:\n\tlibamtrack = ctypes.cdll.LoadLibrary("libamtrack.so")\n\npy_version = python_version()\n'
    header_string +='if int(py_version[0]) < 3 and int(py_version[2]) <= 5:\n\tprint \'ERROR:\\nYou are using Python%s\\nPython2.6 or later needed!\\n\'%py_version\n\tsys.exit(0)\n\n\n'
    
    outfile = open(outfile_name, 'w+')
    outfile.write(header_string)
    
    for c_func in func_objects:
        py_func_string = '\ndef ' + c_func.name + ' (' 
        i=0
        while i < len(c_func.parameter.keys()):
            for  para_name in  c_func.parameter.keys():
                if c_func.parameter[para_name].no == i+1:
                    py_func_string += ' '+ para_name
            if i < len(c_func.parameter.keys())-1:
                py_func_string +=','
            i+=1
        py_func_string+='):\n'
        #write comment, copied from source code
        py_func_string +='\t\'\'\'\n'
        for line in c_func.comment :
            line = string.replace(line ,'\'','\\\'')
            py_func_string += '\t'+ line 
        py_func_string +='\t\'\'\'\n'
        
        for i,  key in enumerate(c_func.parameter.keys()):
            py_func_string += '\t' + translate_type(c_func.parameter[key])
        
        py_func_string += '\tc_function =  libamtrack.' + c_func.name + '\n'
        py_func_string += '\tc_output = c_function('
        i=0
        # loop over positions, to match the one in the function call of libamtrack
        while i < len(c_func.parameter.keys()):
            for  para_name in  c_func.parameter.keys():
                if c_func.parameter[para_name].no == i+1:
                    py_func_string += ' c_'+ para_name
            if i < len(c_func.parameter.keys())-1:
                py_func_string +=','
            i+=1
        py_func_string +=')\n'
            
        c_output_list= []
        for para_name in c_func.parameter.keys():
            if c_func.parameter[para_name].direction == 'out':
                c_output_list.append('c_' +para_name)

        if c_output_list != []:
            py_func_string += '\treturn '
            for i , item in enumerate(c_output_list):
                py_func_string += ' ' + item+'._obj.value'
                if i<len(c_output_list)-1:
                    py_func_string += ','
        py_func_string+='\n\n'
        outfile.write(py_func_string)

    outfile.close()

def translate_type(parameter):
    '''
    takes an amtrack_para object and returns the python lines necessary to access it from python
    '''
    c_lookup = {'double': 'ctypes.c_double', 
                'long': 'ctypes.c_long', 
                'bool': 'ctypes.c_bool',
                'int' : 'ctypes.c_int'}
    
    py_line =''
    if parameter.direction == 'in':
        if not parameter.array:
            py_line += 'c_' + parameter.name + ' = '
            py_line += c_lookup[parameter.type] + '(' +parameter.name+')'
        if parameter.array:
            py_line += 'tmp_array =' +c_lookup[parameter.type] +'* len('+parameter.name + ')\n'
            py_line += '\tc_'+ parameter.name + ' = ctypes.byref(tmp_array(*' +parameter.name +'))\n'
    
    if parameter.direction == 'out':
        if not parameter.array:
            py_line += 'c_'+parameter.name +' =  ctypes.byref('+ c_lookup[parameter.type] +'(' +parameter.name+'))\n'
        if parameter.array :
            py_line += 'tmp_array = ' + c_lookup[parameter.type] + '*' +parameter.array_size + '\n'
            py_line += '\tc_' +parameter.name + '= ctypes.byref(tmp_array)\n'
        
        #print 'array, TODO ' ,  parameter.name,  parameter.array_size
    py_line += '\n'
    return py_line
    

    
def harvest(functions,  path = '../../../include/'):
    '''
    harvest, harvesting information on functions given in \"functions\" in the doc strings
    @param functions[in]: list of functions to be translated, if empty all functions will be wrapped
    @param functions[in]: path to the header directory of the libamtrack library
    @return : list of function objects
    @todo : make path independent from OS (there might be windows users )
        '''
    desired_func = []
    all  = False  
    if functions== []:
        all =True
    for file_name in os.listdir(path):
        tmp_comment = ''
        comment_start = int
        comment_stop = int
        function_end = False
        comment_end = False
        tmp_name = ''
        tmp_definition = ''
        tmp_definition_start = int
        tmp_definition_stop = int
        if file_name [-2:] == '.h':
            tmp_header = open(path+file_name, 'r')
            tmp_header_content = tmp_header.readlines()
            tmp_header.close()        
        for lineno, line in enumerate(tmp_header_content):
            
                splitted_line = string.split(line)
                if splitted_line != []:
                    if function_end:
                        function_end = False
                        #print tmp_definition_start,  len(tmp_header_content)
                        #print tmp_header_content[tmp_definition_start]
                        #print string.split(tmp_header_content[tmp_definition_start])[1]
                        #print string.split(string.split(tmp_header_content[tmp_definition_start])[1], '(') 
                        if all:
                            desired_func.append(tools.amtrack_func(string.split(string.split(tmp_header_content[tmp_definition_start] )[1] , '(') [0]))
                            desired_func[-1].comment = tmp_header_content[comment_start:comment_stop]
                            desired_func[-1].definition = tmp_header_content[tmp_definition_start:tmp_definition_stop]
                            tmp_comment = ''
                            tmp_definition =''
                        
                        elif string.split(string.split(tmp_header_content[tmp_definition_start])[1], '(') [0] in functions:              
                            desired_func.append(tools.amtrack_func(string.split(string.split(tmp_header_content[tmp_definition_start])[1], '(') [0]))
                            desired_func[-1].comment = tmp_header_content[comment_start:comment_stop]
                            desired_func[-1].definition = tmp_header_content[tmp_definition_start:tmp_definition_stop]
                            tmp_comment = ''
                            tmp_definition = ''
                            
                    if splitted_line[0] == '/**': #comment starts
                        comment_start = lineno+1
                        if comment_end:
                            comment_end = False
                    if splitted_line[-1] == '*/': #comment ends
                        comment_stop = lineno
                        comment_end =True
                    if len(string.split(line,'(') )==2  and comment_end: #function definition starts here
                        tmp_definition_start = lineno    
                        comment_end =False
#                        if tmp_line.strip()[-2:] == ');' :
#                          function_end= True
#                           tmp_definition_stop = lineno+1
                    tmp_line = line    # copy, since the following strip command would change the line
                    if tmp_line.strip()[-2:] == ');' and not function_end : #function definition ends here
                        tmp_definition_stop = lineno+1
                        function_end = True
    
    bad_list = []
    for i, item in enumerate(desired_func):
        try:
            item.extract_params()    # analyze functions
        except:
            print item.name  + ': Problems have occurred, check doc-string! Function not included in wrapper'
            bad_list.append(i)
    for  counter , bad_func in enumerate(bad_list): #Ignore non working functions
        desired_func = desired_func[:bad_func -counter] + desired_func[bad_func-counter:]
        
    
    if DEBUG:
        for func in  desired_func:
            for item in func.parameter.keys():
                print func.parameter[item].__dict__
    
    return desired_func
    

    
if __name__ == "__main__":
    print 'py_wrapper tools\nrun generate_Py_wrapper.py in order to generate the Python wrapper'

    
    
    
    
