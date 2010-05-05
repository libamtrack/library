/*
 * AT_error.c
 *
 *  Created on: 03.05.2010
 *      Author: Steffen
 */

#include "AT_Error.h"

char* AT_get_error_msg(const long error_no){
  return "Look up error_no in error_messages and return String.";
}

long AT_check_E_MeV_u(   const long n,
                        const double E_MeV_u[],
                        const long purpose_energy_range)
{
  double E_max_MeV_u;
  double E_min_MeV_u;

  switch(purpose_energy_range){
  case AT_energy_range_for_PSTAR_data:       {     E_min_MeV_u     = 0.001;
                      E_max_MeV_u     = 1.0e4;
                      break;}
  case AT_energy_range_for_PowerLaw_data:    {     E_min_MeV_u     = 1.0;
                      E_max_MeV_u     = 250.0;
                      break;}
  case AT_energy_range_for_Katz_method:      {     E_min_MeV_u     = 0.1;
                      E_max_MeV_u     = 1000;
                      break;}
  case AT_energy_range_for_SPIFF_method:     {     E_min_MeV_u     = 3;
                      E_max_MeV_u     = 300.0;
                      break;}
  default:                              {     E_min_MeV_u     = 1.0;
                      E_max_MeV_u     = 250.0;
                      break;}
  };


  long i;
  for (i = 0; i < n; i++)
    {
      if ((E_MeV_u[i] < E_min_MeV_u) | (E_MeV_u[i] > E_max_MeV_u)){
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
        return AT_Particle_Not_Defined;
      }
    }
  return AT_Success;
}

