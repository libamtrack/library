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
  AT_Material_Already_Established       = 1,
  AT_Energy_Outside_Range               = 2
};

typedef struct {
  const long    error_no;
  const char*   error_msg;
} AT_error_msg;

static const AT_error_msg error_messages[2] = {
    {AT_Success,                        "Success"},
    {AT_Material_Already_Established,   "The material has already been estalished."}
};

char* AT_get_error_msg(const long error_no);

long AT_check_E_MeV_u(   const long n,
                        const double E_MeV_u[]);
#endif

/* AT_ERROR_H_ */
