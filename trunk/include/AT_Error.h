/*
 * AT_Error.h
 *
 *  Created on: 02.05.2010
 *      Author: Steffen
 */

#ifndef AT_ERROR_H_
#define AT_ERROR_H_

enum AT_error_no{
  AT_Success                            = 0,
  AT_Material_Already_Established,
  AT_Energy_Outside_Range,
  AT_Particle_Not_Defined,
  AT_No_PSTAR_Data,
  AT_Unknown_LET_Data_Source
};

typedef struct {
  const long    error_no;
  const char*   error_msg;
} AT_error_msg;

static const AT_error_msg error_messages[6] = {
    {AT_Success,                        "Success"},
    {AT_Material_Already_Established,   "The material has already been established."},
    {AT_Energy_Outside_Range,           "The energies are outside the range of libamtrack (1 - 250 MeV/u)."},
    {AT_Particle_Not_Defined,           "The particle definition is not correct."},
    {AT_No_PSTAR_Data,                  "No PSTAR data for this material (yet)."},
    {AT_Unknown_LET_Data_Source,        "The data source for LET data is not available."}
};

char* AT_get_error_msg(const long error_no);

long AT_check_E_MeV_u(   const long n,
                        const double E_MeV_u[]);
long AT_check_particle_no(   const long n,
                             const long particle_no[]);
#endif

/* AT_ERROR_H_ */
