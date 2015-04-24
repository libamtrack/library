#ifndef AT_STOPPINGPOWERDATA_H_
#define AT_STOPPINGPOWERDATA_H_

/**
 * @brief Stopping power
 */

/*
 *    AT_DataStoppingPower.h
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

#include "AT_StoppingPowerDataBethe.h"
#include "AT_StoppingPowerDataPSTAR.h"
#include "AT_StoppingPowerDataICRU.h"
#include "AT_StoppingPowerDataFromFile.h"


#define STOPPING_POWER_SOURCE_N					4
#define STOPPING_POWER_SOURCE_NAME_LENGTH    	255

/**
 * @enum stoppingPowerSource_no 	Stopping-power data-source ids
 */
enum stoppingPowerSource_no{
	FromFile             = 0, /**< Data read from file */
	Bethe                = 1, /**< Analytical Bethe fomula */
	PSTAR                = 2, /**< PSTAR data from NIST */
	ICRU                 = 3, /**< ICRU 49 and 73 data*/
};

typedef int (*stopping_power_function)(const long,
		const double[],
		const long[],
		const long,
		const char[],
		double[]);

/**
 * Structure to connect stopping power number and corresponding function
 *
 * @struct AT_stopping_power_functions_struct
 */
typedef struct {
	const char stopping_power_source_name[STOPPING_POWER_SOURCE_N][STOPPING_POWER_SOURCE_NAME_LENGTH];
	stopping_power_function function[STOPPING_POWER_SOURCE_N];
} AT_stopping_power_functions_struct;


extern AT_stopping_power_functions_struct AT_stopping_power_functions;

/**
 * TODO
 * @param[in]  source_no
 * @param[out] source_name
 * @return statuc code
 */
int AT_stopping_power_source_model_name_from_number( const long source_no,
		char* source_name);

/**
 * TODO
 * @param[in] source_name
 * @return    source number
 */
long AT_stopping_power_source_model_number_from_name( const char* source_name);

#endif /* AT_STOPPINGPOWERDATA_H_ */
