#ifndef AT_ERROR_H_
#define AT_ERROR_H_

/**
 * @file
 * @brief error codes
 */

/*
 *    AT_Error.h
 *    ==============
 *
 *    Created on: 03.05.2010
 *    Author: greilich
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


/**
 * TODO
 */
enum AT_error_no{
  AT_Success                            = 0,
  AT_Material_Already_Established,
  AT_Energy_Outside_Range,
  AT_Particle_Not_Defined,
  AT_No_PSTAR_Data,
  AT_Unknown_LET_Data_Source
};


/**
 * TODO
 */
enum AT_energy_ranges{
  AT_energy_range_for_PSTAR_data,
  AT_energy_range_for_PowerLaw_data,
  AT_energy_range_for_Katz_method,
  AT_energy_range_for_SPIFF_method
};


/**
 * TODO
 */
typedef struct {
  const long    error_no;
  const char*   error_msg;
} AT_error_msg;


/**
 * TODO
 */
static const AT_error_msg error_messages[6] = {
    {AT_Success,                        "Success"},
    {AT_Material_Already_Established,   "The material has already been established."},
    {AT_Energy_Outside_Range,           "The energies are outside the range of the selected purpose."},
    {AT_Particle_Not_Defined,           "The particle definition is not correct."},
    {AT_No_PSTAR_Data,                  "No PSTAR data for this material (yet)."},
    {AT_Unknown_LET_Data_Source,        "The data source for LET data is not available."}
};


/**
 * TODO
 * @param error_no
 * @return
 */
char* AT_get_error_msg(const long error_no);


/**
 * TODO
 * @param n
 * @param E_MeV_u
 * @param purpose_energy_range
 * @return
 */
long AT_check_E_MeV_u(   const long n,
                        const double E_MeV_u[],
                        const long purpose_energy_range);


/**
 * TODO
 * @param n
 * @param particle_no
 * @return
 */
long AT_check_particle_no(   const long n,
                             const long particle_no[]);

#endif /* AT_ERROR_H_ */
