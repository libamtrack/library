#ifndef AT_SPC_H_
#define AT_SPC_H_

/**
 * @brief SPC reader
 */

/*
 *    AT_SPC.h
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
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>


/**
 * File name and path length
 */
#define FILE_NAME_NCHAR 256


/**
 * @enum spc tags
 */
enum {
       TRPSPCBTAG_FILETYPE   =1,      /** < header info *///!< TRPSPCBTAG_FILETYPE
       TRPSPCBTAG_FILEVERSION=2,                          //!< TRPSPCBTAG_FILEVERSION
       TRPSPCBTAG_FILEDATE   =3,                          //!< TRPSPCBTAG_FILEDATE
       TRPSPCBTAG_TARGNAME   =4,                          //!< TRPSPCBTAG_TARGNAME
       TRPSPCBTAG_PROJNAME   =5,                          //!< TRPSPCBTAG_PROJNAME
       TRPSPCBTAG_B          =6,                          //!< TRPSPCBTAG_B
       TRPSPCBTAG_P          =7,                          //!< TRPSPCBTAG_P
       TRPSPCBTAG_N          =8,                          //!< TRPSPCBTAG_N
       TRPSPCBTAG_NZ         =9,                          //!< TRPSPCBTAG_NZ
       TRPSPCDTAG_Z          =10,                         //!< TRPSPCDTAG_Z
       TRPSPCDTAG_N          =11,                         //!< TRPSPCDTAG_N
       TRPSPCDTAG_NS         =12,                         //!< TRPSPCDTAG_NS
       TRPSPCDTAG_S          =13,                         //!< TRPSPCDTAG_S
       TRPSPCDTAG_CUM        =14,                         //!< TRPSPCDTAG_CUM
       TRPSPCDTAG_NC         =15,                         //!< TRPSPCDTAG_NC
       TRPSPCDTAG_NE         =16,                         //!< TRPSPCDTAG_NE
       TRPSPCDTAG_E          =17,                         //!< TRPSPCDTAG_E
       TRPSPCDTAG_EREF       =18,                         //!< TRPSPCDTAG_EREF
       TRPSPCDTAG_HISTO      =19,                         //!< TRPSPCDTAG_HISTO
       TRPSPCDTAG_RUNNINGSUM =20,                         //!< TRPSPCDTAG_RUNNINGSUM
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
 * Reads data from spc file. Will reallocate the arrays to hold the
 * read content - and can therefore only called from environments
 * that are able to cope with reallocation (C, C++, but not R!).
 * The data will be converted for direct use in libamtrack, i.e. to an R-style array of six columns
 * (each presented by a single pointer) and all of same length. That of course
 * results in redundancy in depth_step, depth_g_cm2, particle_no but enables
 * easy division into cells (i.e. passing all spectra of a specific depth and
 * particle number to another routine such as total dose). Please note that
 * the fluence IS NOT normalized to bin width but given in absolute fluence!
 *
 * @param[in]	    file  	        path and name of spc file to open (incl. suffix)
 * @param[out]		depth_step		depth step index, zero-based (pointer to array of initial size 1)
 * @param[out]		depth_g_cm2		depth in g/cm2 (pointer to array of initial size 1)
 * @param[out]		E_MeV_u			midpoints of energy bins (pointer to array of initial size 1)
 * @param[out]		DE_MeV_u		widths of energy bins (pointer to array of initial size 1)
 * @param[out]		particle_no		particle index numbers (pointer to array of initial size 1)
 * @param[out]      fluence_cm2		fluence values differential in energy and particle number
 * @return                          array sizes after reallocation
 */
int AT_SPC_read( const char filename[FILE_NAME_NCHAR],
		int* depth_step[],
		double* depth_g_cm2[],
		double* E_MeV_u[],
		double* DE_MeV_u[],
		int* particle_no[],
		double* fluence_cm2[]);


/**
 * Browses spc file to get total number of energy bins covered.
 * This is needed for later memory allocation when reading the
 * actual data.
 *
 * @param[in]	    filename  	    path and name for spc file (incl. extension)
 * @return							total number of bins in spc file
 */
int AT_SPC_get_size_from_filename(const char filename[FILE_NAME_NCHAR]);


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
 * @param[in]	    filename  	    path and name for spc file (incl. extension)
 * @param[in]       n               array size, total number of bins expected
 * @see AT_SPC_get_size
 * @param[out]		depth_step		depth step index, zero-based (array of size n)
 * @param[out]		depth_g_cm2		depth in g/cm2 (array of size n)
 * @param[out]		E_MeV_u			midpoints of energy bins (array of size n)
 * @param[out]		DE_MeV_u		widths of energy bins (array of size n)
 * @param[out]		particle_no		particle index numbers (array of size n)
 * @param[out]      fluence_cm2		fluence values differential in energy and particle number
 * @return                          number of bins read. Must match the array size n
 */
int AT_SPC_read_data_from_filename( const char filename[FILE_NAME_NCHAR],
		const int n,
		int depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		int particle_no[],
		double fluence_cm2[]);


/**
 * Reads data from spc file into pre-allocated arrays. It will be converted
 * for direct use in libamtrack, i.e. to an R-style array of six columns
 * (each presented by a single pointer) and all of same length. That of course
 * results in redundancy in depth_step, depth_g_cm2, particle_no but enables
 * easy division into cells (i.e. passing all spectra of a specific depth and
 * particle number to another routine such as total dose). Please note that
 * the fluence IS NOT normalized to bin width but given in absolute fluence!
 *
 * @param[in]	    fp  	        pointer to spc file
 * @param[in]       n               array size, total number of bins expected
 * @see AT_SPC_get_size
 * @param[out]		depth_step		depth step index, zero-based (array of size total_n_bins)
 * @param[out]		depth_g_cm2		depth in g/cm2 (array of size total_n_bins)
 * @param[out]		E_MeV_u			midpoints of energy bins (array of size total_n_bins)
 * @param[out]		DE_MeV_u		widths of energy bins (array of size total_n_bins)
 * @param[out]		particle_no		particle index numbers (array of size total_n_bins)
 * @param[out]      fluence_cm2		fluence values differential in energy and particle number
 */
int AT_SPC_read_data( FILE* fp,
		const int n,
		int depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		int particle_no[],
		double fluence_cm2[]);


#endif /* AT_SPC_H_ */
