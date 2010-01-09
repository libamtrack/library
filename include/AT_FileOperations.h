#ifndef AT_FILEOPERATIONS_H_
#define AT_FILEOPERATIONS_H_

/**
*    AT_Functions.h
*    ==============
*
*    Copyright 2006, 2009 Steffen Greilich / the libamtrack team
*
*    This file is part of the AmTrack program (libamtrack.sourceforge.net).
*
*    AmTrack is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    AmTrack is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with AmTrack (file: copying.txt).
*    If not, see <http://www.gnu.org/licenses/>
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "AT_NumericalRoutines.h"

void AT_browseInput( const char*  fileName,
    // return:
    long*  nLines,
    long*  n_gamma_parameter);

void AT_browseInputS(  char**  fileName,
    // return:
    long*  nLines,
    long*  n_gamma_parameter);


void AT_readInput( const  char*  fileName,
    const long*  nLines,
    const long*  n_gamma_parameter,
    // return:
    float*  E_MeV_u,
    long*  particle_no,
    float*  fluence_cm2,
    long*  slab_no,
    float*  parameter,
    long*  N2,
    char*  material_name,
    long*  n_slabs,
    long*  gamma_model,
    float*  gamma_parameter);

void AT_readInputS(char**  fileName,
    long*  nLines,
    long*  n_gamma_parameter,
    // return:
    float*  E_MeV_u,
    long*  particle_no,
    float*  fluence_cm2,
    long*  slab_no,
    float*  parameter,
    long*  N2,
    char**  material_name,
    long*  n_slabs,
    long*  gamma_model,
    float*  gamma_parameter);

void AT_browseSpectrum( const char* fileName,
    // return:
    long*  nLines);

void AT_browseSpectrumS(  char**  fileName,
    // return:
    long*  nLines);


void AT_readSpectrum( const  char* fileName,
    const long* nLines,
    // return:
    float*  E_MeV_u,
    long*  particle_no,
    float*  fluence_cm2,
    long*  slab_no);

void AT_readSpectrumS(  char**  fileName,
    long*  nLines,
    // return:
    float*  E_MeV_u,
    long*  particle_no,
    float*  fluence_cm2,
    long*  slab_no);


/////////////////////////////////////////////////////////////////////////////////

#endif // AT_FILEOPERATIONS_H_
