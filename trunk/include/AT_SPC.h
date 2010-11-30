#ifndef AT_SPC_H_
#define AT_SPC_H_

/**
 * @brief LET tables and access routines
 */

/*
 *    AT_DataLET.h
 *    ==================
 *
 *    Copyright 2006, 2010 The libamtrack team
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
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

enum {
       TRPSPCBTAG_FILETYPE   =1,      /* header info */
       TRPSPCBTAG_FILEVERSION=2,
       TRPSPCBTAG_FILEDATE   =3,
       TRPSPCBTAG_TARGNAME   =4,
       TRPSPCBTAG_PROJNAME   =5,
       TRPSPCBTAG_B          =6,
       TRPSPCBTAG_P          =7,
       TRPSPCBTAG_N          =8,
       TRPSPCBTAG_NZ         =9,
       TRPSPCDTAG_Z          =10,
       TRPSPCDTAG_N          =11,
       TRPSPCDTAG_NS         =12,
       TRPSPCDTAG_S          =13,
       TRPSPCDTAG_CUM        =14,
       TRPSPCDTAG_NC         =15,
       TRPSPCDTAG_NE         =16,
       TRPSPCDTAG_E          =17,
       TRPSPCDTAG_EREF       =18,
       TRPSPCDTAG_HISTO      =19,
       TRPSPCDTAG_RUNNINGSUM =20,
};

struct STRPSPCBTAG { uint32_t ulTag; uint32_t ulLen; };

void endian_swap2(    unsigned short* x);
void endian_swap4(    unsigned int* x);
// __int64 for MSVC, "long long" for gcc
void endian_swap8(     uint64_t * x);
void readStruct(      FILE *fp,
		struct STRPSPCBTAG * str,
		bool switchEndian);

int AT_SPC_get_size(  FILE *fp);
int AT_SPC_read_data( FILE *fp,
		long depthStep[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		long particle_no[],
		double fluence_cm2[]);

#endif /* AT_SPC_H_ */
