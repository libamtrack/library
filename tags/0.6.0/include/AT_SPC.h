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
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>

#include "AT_DataParticle.h"
#include "AT_DataMaterial.h"

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
 * TODO
 * @param[in]	filename  	    	path and name for spc file, incl. extension
 * @return	    number of bytes in file
 */
int AT_SPC_get_number_of_bytes_in_file( const char filename[] );


/**
 * TODO
 * @param[in]	filename  	    	path and name for spc file, incl. extension
 * @param[in]	    content_size 	number of elements in content arrays
 * @param[out]      content  	    binary content of SPC file (array of size content_size)
 * @return			status code
 */
int AT_SPC_fast_read_buffer( const char filename[],
		int content_size,
		int32_t* content);


/**
 * TODO
 * @param[in]	    content  	    table of bytes containing binary content of SPC file
 * @param[out]	    string 			string stored in first tag (item) of SPC content table
 * @param[out]      length  	    length of string
 */
void decomposeStructIntoString( const int32_t content[],
		char* string,
		int*  length );


/**
 * TODO
 * @param[in]	    content  	    table of bytes containing binary content of SPC file
 * @param[out]	    value 			floating point value (single or pointer to beginning of table) stored in first tag (item) of SPC content table
 * @param[out]      length  	    number of items under "value" pointer
 */
void decomposeStructIntoDouble( const int32_t content[],
		double* value,
		int*    length );


/**
 * TODO
 * @param[in]	    content  	    table of bytes containing binary content of SPC file
 * @param[out]	    value 			integer value (single or pointer to beginning of table) stored in first tag (item) of SPC content table
 * @param[out]      length  	    number of items under "value" pointer
 */
void decomposeStructIntoInteger( const int32_t content[],
		uint64_t* value,
		int*      length );


/**
 * Increase pointer (binary content of SPC file) to move it to next tag (item).
 * @param[in,out]	content  	    table of bytes containing binary content of SPC file
 */
int skipStruct( int32_t** content );


/**
 * TODO
 * @param[in]	content  	    table of bytes containing binary content of SPC file
 * @return
 */
int decomposeTag( const int32_t content[] );


/**
 * TODO
 * @param[in]	content  	    table of bytes containing binary content of SPC file
 * @return
 */
int decomposeLength( const int32_t content[] );


/**
 * TODO
 * @param[in]	content_size  	    size of table content_orig
 * @param[in]	content_orig  	    table of bytes containing binary content of SPC file (array of size content_size)
 * @return      number of bins
 */
int AT_SPC_decompose_size( const int content_size,
		int32_t content_orig[]);


/**
 * TODO
 * @param[in]	content_size  	    size of table content_orig
 * @param[in]	content_orig  	    table of bytes containing binary content of SPC file (array of size content_size)
 * @param[out]	E_MeV_u  	        beam energy [MeV]
 * @param[out]	peak_position_g_cm2 peak position
 * @param[out]	particle_no         projectile - particle no
 * @param[out]	material_no         target - material no
 * @param[out]	normalisation  	    normalisation
 * @param[out]  depth_steps_no      number of depth steps
 * return       status code
 */
int AT_SPC_decompose_header(		const int content_size,
		int32_t   content_orig[],
		double*   E_MeV_u,
		double*   peak_position_g_cm2,
		long*     particle_no,
		int*      material_no,
		double*   normalisation,
		int*      depth_steps_no);


/**
 * TODO
 * @param[in]	content_size  	    size of table content_orig
 * @param[in]	content_orig  	    table of bytes containing binary content of SPC file (array of size content_size)
 * @param[out]  depth_step          depth step index, zero-based
 * @param[out]  depth_g_cm2         depth in g/cm2
 * @param[out]	E_MeV_u  	        midpoints of energy bins
 * @param[out]	DE_MeV_u  	        widths of energy bins
 * @param[out]	particle_no         particle index numbers
 * @param[out]	fluence_cm2  	    fluence values differential in energy and particle number
 * return       status code
 */
int AT_SPC_decompose_data(		const int content_size,
		int32_t   content_orig[],
		int*      depth_step[],
		double*   depth_g_cm2[],
		double*   E_MeV_u[],
		double*   DE_MeV_u[],
		long*     particle_no[],
		double*   fluence_cm2[]);


/**
 * TODO
 * @param[in]	filename  	    	path and name for spc file, incl. extension
 * @return               number of bins
 */
long AT_SPC_get_number_of_bins_from_filename_fast( const char filename[] );

/**
 * Reads data from spc file into pre-allocated arrays. It will be converted
 * for direct use in libamtrack, i.e. to an R-style array of six columns
 * (each presented by a single pointer) and all of same length. That of course
 * results in redundancy in depth_step, depth_g_cm2, particle_no but enables
 * easy division into cells (i.e. passing all spectra of a specific depth and
 * particle number to another routine such as total dose). Please note that
 * the fluence IS NOT normalized to bin width but given in absolute fluence!
 *
 * @param[in]	filename  	    	path and name for spc file, incl. extension
 * @param[out]	E_MeV_u			primary beam energy in MeV/u
 * @param[out]	peak_position_g_cm2	position of peak in g/cm2
 * @param[out]	particle_no         projectile - particle no
 * @param[out]	material_no         target - material no
 * @param[out]	normalisation		normalisation
 * @param[out]  depth_steps_no      number of depth steps
 * @return               status code
 */
int AT_SPC_read_header_from_filename_fast( const char filename[],
		double*   E_MeV_u,
		double*   peak_position_g_cm2,
		long*     particle_no,
		int*      material_no,
		double*   normalisation,
		int*      depth_steps_no);

/**
 * Reads data from spc file into pre-allocated arrays. It will be converted
 * for direct use in libamtrack, i.e. to an R-style array of six columns
 * (each presented by a single pointer) and all of same length. That of course
 * results in redundancy in depth_step, depth_g_cm2, particle_no but enables
 * easy division into cells (i.e. passing all spectra of a specific depth and
 * particle number to another routine such as total dose). Please note that
 * the fluence IS NOT normalized to bin width but given in absolute fluence!
 *
 * @param[in]	    filename  	    path and name for spc file, incl. extension
 * @param[in]       n               array size, total number of bins expected
 * @see AT_SPC_get_size
 * @param[out]		depth_step		depth step index, zero-based (array of size n)
 * @param[out]		depth_g_cm2		depth in g/cm2 (array of size n)
 * @param[out]		E_MeV_u			midpoints of energy bins (array of size n)
 * @param[out]		DE_MeV_u		widths of energy bins (array of size n)
 * @param[out]		particle_no		particle index numbers (array of size n)
 * @param[out]      fluence_cm2		fluence values differential in energy and particle number (array of size n)
 * @return                          number of bins read. Must match the array size n
 */
int AT_SPC_read_data_from_filename_fast( const char filename[],
		int    n,
		int    depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		long   particle_no[],
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
 * @param[in]	filename  	    	path and name for spc file, incl. extension
 * @param[in]   n                   array size, total number of bins expected
 * @see AT_SPC_get_size
 * @param[out]	E_MeV_u_initial		primary beam energy in MeV/u
 * @param[out]	peak_position_g_cm2	position of peak in g/cm2
 * @param[out]	particle_no_initial projectile - particle no
 * @param[out]	material_no         target - material no
 * @param[out]	normalisation		normalisation
 * @param[out]  depth_steps_no      number of depth steps
 * @param[out]	depth_step		    depth step index, zero-based (array of size n)
 * @param[out]	depth_g_cm2		    depth in g/cm2 (array of size n)
 * @param[out]	E_MeV_u			    midpoints of energy bins (array of size n)
 * @param[out]	DE_MeV_u		    widths of energy bins (array of size n)
 * @param[out]	particle_no		    particle index numbers (array of size n)
 * @param[out]  fluence_cm2		    fluence values differential in energy and particle number (array of size n)
 * @return                          number of bins read. Must match the array size n
 */
int AT_SPC_read_from_filename_fast( const char filename[],
		int 	  n,
		double*   E_MeV_u_initial,
		double*   peak_position_g_cm2,
		long*     particle_no_initial,
		int*      material_no,
		double*   normalisation,
		int*      depth_steps_no,
		int       depth_step[],
		double    depth_g_cm2[],
		double    E_MeV_u[],
		double    DE_MeV_u[],
		long      particle_no[],
		double    fluence_cm2[]);
	

struct spc_pair {
	double range_gcm2;
	char filename[2048];
};


/**
 * TODO
 * @param[in]   a  TODO
 * @param[in]   b  TODO
 * @return  TODO
 */
int compare_SPC_Pairs (const void *a,
		const void *b);


/**
 * TODO
 *
 * @param[in]   path            	path to spc file dir (array of size 2048)
 * @param[in]   range_g_cm2         range in g/cm2
 * @return                          number of bins for given range
 */
long AT_SPC_number_of_bins_at_range( const char path[],
        double range_g_cm2);


/**
 * TODO
 *
 * @param[in]   path            	path to spc file dir (array of size 2048)
 * @param[in]   range_g_cm2         range in g/cm2
 * @param[in]   n                   array size, total number of bins expected
 * @param[out]	depth_step		    depth step index, zero-based (array of size n)
 * @param[out]	depth_g_cm2		    depth in g/cm2 (array of size n)
 * @param[out]	E_MeV_u			    midpoints of energy bins (array of size n)
 * @param[out]	DE_MeV_u		    widths of energy bins (array of size n)
 * @param[out]	particle_no		    particle index numbers (array of size n)
 * @param[out]  fluence_cm2		    fluence values differential in energy and particle number (array of size n)
 * @return                          status code
 */
int AT_SPC_spectrum_at_range( const char path[],
		double range_g_cm2,
		int n,
		int    depth_step[],
		double depth_g_cm2[],
		double E_MeV_u[],
		double DE_MeV_u[],
		long   particle_no[],
		double fluence_cm2[]);


#endif /* AT_SPC_H_ */
