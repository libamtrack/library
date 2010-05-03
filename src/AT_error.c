/*
 * AT_error.c
 *
 *  Created on: 03.05.2010
 *      Author: Steffen
 */

#include "AT_error.h"

char* AT_get_error_msg(const long error_no){
  return "Look up error_no in error_messages and return String.";
}

long AT_check_E_MeV_u(   const long n,
                        const double E_MeV_u[])
{
  long i;
  for (i = 0; i < n; i++)
    {
      if ((E_MeV_u[i] < 1) | (E_MeV_u[i] > 250)){
        return AT_Energy_Outside_Range;
      }
    }
  return AT_Success;
}

long AT_check_particle_no(   const long n,
                             const long particle_no[])
{
  long i;
  for (i = 0; i < n; i++)
    {
      if (      (particle_no[i] < 0) |
                (AT_Z_from_particle_no_single(particle_no[i]) < 1) |
                (AT_Z_from_particle_no_single(particle_no[i]) > 98) |
                (AT_A_from_particle_no_single(particle_no[i]) < 1) |
                (AT_A_from_particle_no_single(particle_no[i]) > 251))
      {
        return AT_Energy_Outside_Range;
      }
    }
  return AT_Success;
}

