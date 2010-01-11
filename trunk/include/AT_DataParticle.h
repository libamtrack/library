#ifndef AT_DATAPARTICLE_H_
#define AT_DATAPARTICLE_H_

/**
*    AT_DataParticle.h
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
#include <stdlib.h>

#include "AT_NumericalRoutines.h"

#define PARTICLE_DATA_N    35

typedef struct {
  long    n;
  long    particle_no[PARTICLE_DATA_N];
  char*    particle_name[PARTICLE_DATA_N];
  char*    element_name[PARTICLE_DATA_N];
  long    Z[PARTICLE_DATA_N];
  long    A[PARTICLE_DATA_N];
  float    mass[PARTICLE_DATA_N];
  char*    USRTRACK_name[PARTICLE_DATA_N];
} particle_data;

static particle_data AT_Particle_Data = {
    PARTICLE_DATA_N ,
    {
        1,              2,              3,
        4,              5,
        6,              7,              8,
        9,              10,             11,
        12,             13,             14,             15,
        16,             17,             18,             19,             20,
        21,             22,             23,             24,             25,
        26,             27,             28,             29,             30,             31,             32,
        33,             34,             35
    },
    {
        "P"   ,    "D",    "T",
        "He-3", "He-4",
        "Li-6", "Li-7", "Li-8",
        "Be-7", "Be-8", "Be-9",
        "B-10", "B-11", "B-12", "B-13",
        "C-10", "C-11", "C-12", "C-13", "C-14",
        "N-13", "N-14", "N-15", "N-16", "N-17",
        "O-14", "O-15", "O-16", "O-17", "O-18", "O-19", "O-20",
        "Ne-20","Si-28","Fe-56"
    },
    {
        "H" ,    "H",    "H",
        "He",   "He",
        "Li",   "Li",   "Li",
        "Be",   "Be",   "Be",
        "B" ,    "B",    "B",    "B",
        "C" ,    "C",    "C",    "C",    "C",
        "N" ,    "N",    "N",    "N",    "N",
        "O" ,    "O",    "O",    "O",    "O",    "O",    "O",
        "Ne",   "Si",   "Fe"
    },
    {
        1,              1,              1,                                                                                      // H
        2,              2,                                                                                                      // He
        3,              3,              3,                                                                                      // Li
        4,              4,              4,                                                                                      // Be
        5,              5,              5,              5,                                                                      // B
        6,              6,              6,              6,              6,                                                      // C
        7,              7,              7,              7,              7,                                                      // N
        8,              8,              8,              8,              8,              8,              8,                      // O
        10,             14,             26
    },
    {
        1,              2,              3,                                                                                      // H
        3,              4,                                                                                                      // He
        6,              7,              8,                                                                                      // Li
        7,              8,              9,                                                                                      // Be
        10,             11,             12,             13,                                                                     // B
        10,             11,             12,             13,             14,                                                     // C
        13,             14,             15,             16,             17,                                                     // N
        14,             15,             16,             17,             18,             19,             20,                     // O
        20,             28,             56
    },
    {
        1.0078,         2.0141,         3.0161,                                                                         // H
        3.0160,         4,                                                                                                      // He
        6,              7,              8,                                                                                      // Li
        7,              8,              9,                                                                                      // Be
        10,             11,             12,             13,                                                                     // B
        10,             11,             12.0000,        13,             14,                                                     // C
        13,             14,             15,             16,             17,                                                     // N
        14,             15,             15.9949,        17,             18,             19,             20, // O
        19.992, 27.9769265,             55.9349
    },
    {
        "P       ",     "D       ",     "T       ",
        "HE3     ",     "HE4     ",
        "LI6     ",     "LI7     ",     "LI8     ",
        "BE7     ",     "BE8     ",     "BE9     ",
        "B10     ",     "B11     ",     "B12     ",     "B13     ",
        "C10     ",     "C11     ",     "C12     ",     "C13     ",     "C14     ",
        "N13     ",     "N14     ",     "N15     ",     "N16     ",     "N17     ",
        "O14     ",     "O15     ",     "O16     ",     "O17     ",     "O18     ",     "O19     ",     "O20     ",
        "NE20    ", "SI28    ", "FE56    "}
};

void AT_mass_from_particle_no(  const long*  n,
    const long*  particle_no,
    float*  mass);

void AT_A_from_particle_no(  const long*  n,
    const long*  particle_no,
    long*  A);

void AT_Z_from_particle_no(  const long*  n,
    const long*  particle_no,
    long*  Z);


void AT_Particle_Properties( const long*  n,
    const long*  particle_no,
    /* return values*/
    char**  particle_name,
    char**  USRTRACK_name,
    char**  element_name,
    long*  Z,
    long*  A,
    float*  mass);

#endif /* AT_DATAPARTICLE_H_ */
