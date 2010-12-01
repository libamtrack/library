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

/**
 * @enum spc tags
 */
enum {
       TRPSPCBTAG_FILETYPE   =1,      /** < header info */
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

/**
 * structure for spc tags
 */

struct STRPSPCBTAG { uint32_t ulTag; uint32_t ulLen; };

/**
 * Swaps endianess of two-byte types
 *
 * @param[in/out]	x	pointer to variable to be swapped
 */
void endian_swap2(    unsigned short* x);

/**
 * Swaps endianess of four-byte types
 *
 * @param[in/out]	x	pointer to variable to be swapped
 */
void endian_swap4(    unsigned int* x);

/**
 * Swaps endianess of eight-byte types (uint64_t MSVC, long long in gcc)
 *
 * @param[in/out]	x	pointer to variable to be swapped
 */
void endian_swap8(     uint64_t * x);

/**
 * Reads tag into tag structure
 *
 * @param[in/out]	fp				pointer to spc file
 * @param[in/out]	str				pointer to spc tag structure to be filled
 * @param[in]		switchEndian	if TRUE endianess will be swapped
 */
void readStruct(      FILE *fp,
		struct STRPSPCBTAG * str,
		bool switchEndian);

/**
 * Browses spc file to get total number of energy bins covered.
 * This is needed for later memory allocation when reading the
 * actual data.
 *
 * @param[in/out]	fp				pointer to spc file
 * @return							total number of bins in spc file
 */
int AT_SPC_get_size(  FILE *fp);

/**
 * Reads data from spc file into pre-allocated arrays. It will be converted
 * for direct use in libamtrack, i.e. to an R-style array of six columns
 * (each presented by a single pointer) and all of same length. That of course
 * results in redundancy in depth_step, depth_g_cm2, particle_no but enables
 * easy division into cells (i.e. passing all spectra of a specific depth and
 * particle number to another routine such as total dose). Please note that
 * the fluence IS NOT normalized to bin width but given in absolute fluence!
 *
 * @param[in/out]	fp				pointer to spc file
 * @param[in]		total_n_bins	size of pre-allocated arrays
 * @see AT_SPC_get_size
 * @param[out]		depth_step		depth step index, zero-based (array of size total_n_bins)
 * @param[out]		depth_g_cm2		depth in g/cm2 (array of size total_n_bins)
 * @param[out]		E_MeV_u			midpoints of energy bins (array of size total_n_bins)
 * @param[out]		DE_MeV_u		widths of energy bins (array of size total_n_bins)
 * @param[out]		particle_no		particle index numbers (array of size total_n_bins)
 * @param[out]      fluence_cm2		fluence values differential in energy and particle number
 */
int AT_SPC_read_data( FILE *fp,
		const long total_n_bins,
		long depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		long particle_no[],
		double fluence_cm2[]);

#endif /* AT_SPC_H_ */
