/**
 * @file
 * @brief LET tables
 */

/*
 *    AT_DataLET.c
 *    ==============
 *
 *    Created on: 09.01.2010
 *    Creator: kongruencja
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

#include "AT_DataLET.h"
#include "AT_DataParticle.h"
#include <math.h>

void get_table_value(
    const long    n,
    const double  x[],
    const long    subset_no,
    const double  x_table[],
    const double  y_table[],
    double        y[])
{
  // first: find those PSTAR entries that match the material number
  bool*    matches    =  (bool*)calloc(AT_PSTAR_Data.n, sizeof(bool));
  is_element_int(    subset_no,
      AT_PSTAR_Data.material_no,
      AT_PSTAR_Data.n,
      matches);

  long    n_matches  = 0;
  long    i;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      n_matches++;
    }
  }

  // allocate vectors for extracted LET entries
  double*  x_c  =  (double*)calloc(n_matches, sizeof(double));
  double*  y_c  =  (double*)calloc(n_matches, sizeof(double));

  // and get the values
  long     j  = 0;
  for (i = 0; i < AT_PSTAR_Data.n; i++){
    if (matches[i]){
      x_c[j]  = x_table[i];
      y_c[j]  = y_table[i];
      j++;
    }
  }
  long  n_pol      = 4 + 1;
  for (i = 0; i < n; i++){
    // Get proton-LET for scaled energy from table E, L using 4th degree polynomial (n_pol - 1 = 2) interpolation
    double  err_y_tmp  = 0.0;    // dummy
    interp(    x_c,
        y_c,
        n_matches,
        n_pol,
        x[i],
        &y[i],
        &err_y_tmp);
  }

  free(x_c);
  free(y_c);
  free(matches);
}


double AT_LET_MeV_cm2_g_single(  const double  E_MeV_u,
    const long    particle_no,
    const long    material_no){

  // get LET for proton of same energy / nucleon
  const long number_of_particles  =  1;
  double LET_MeV_cm2_g            =  0.0;
  get_table_value(number_of_particles, &E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, &LET_MeV_cm2_g);

  double Zeff_ion    =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u, particle_no);

  double Zeff_proton =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u, PARTICLE_PROTON_NUMBER);

  // scale proton LET by ratio of effective Z
  if( particle_no != PARTICLE_PROTON_NUMBER){ // for particles other than proton scale LET by (Zeff_ion / Zeff_proton)^2
      LET_MeV_cm2_g *=   gsl_pow_2(Zeff_ion / Zeff_proton);
  }

  return LET_MeV_cm2_g;
}


void AT_LET_MeV_cm2_g(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    double        LET_MeV_cm2_g[])
{
  // get LET for proton of same energy / nucleon
  get_table_value(number_of_particles, E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, LET_MeV_cm2_g);

  // get effective charge for all given particles and energies
  double*  Zeff_ion  =  (double*)calloc( number_of_particles, sizeof(double));
  AT_effective_charge_from_E_MeV_u(  number_of_particles,
      E_MeV_u,
      particle_no,
      Zeff_ion);

  // get effective charge for protons of same energy
  long   i;
  double*  Zeff_proton           =  (double*)calloc(number_of_particles, sizeof(double));
  long*    particle_no_proton    =  (long*)calloc(number_of_particles, sizeof(long));
  for (i = 0; i < number_of_particles; i++){
    particle_no_proton[i] = PARTICLE_PROTON_NUMBER;
  }

  AT_effective_charge_from_E_MeV_u(  number_of_particles,
      E_MeV_u,
      particle_no_proton,
      Zeff_proton);

  // scale proton LET by ratio of effective Z
  for (i = 0; i < number_of_particles; i++){
    if( particle_no[i] != PARTICLE_PROTON_NUMBER){ // for particles other than proton scale LET by (Zeff_ion / Zeff_proton)^2
      LET_MeV_cm2_g[i] *=   gsl_pow_2(Zeff_ion[i] / Zeff_proton[i]);
    }
  }

  free(particle_no_proton);
  free(Zeff_proton);

  free(Zeff_ion);
}


void AT_LET_keV_um(  const long  number_of_particles,
    const double E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    double       LET_keV_um[])
{
  // Get material density
  double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);

  // Get mass-norm. LET
  AT_LET_MeV_cm2_g(  number_of_particles,
      E_MeV_u,
      particle_no,
      material_no,
      LET_keV_um);

  long  i;
  for (i = 0; i < number_of_particles; i++){
    LET_keV_um[i]  *=  material_density_g_cm3 * 0.1;
  }
}


void AT_CSDA_range_g_cm2(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    double        CSDA_range_g_cm2[])
{
  get_table_value(number_of_particles, E_MeV_u, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.range_cdsa_g_cm2, CSDA_range_g_cm2);

  // Conversion CSDA_proton => CSDA_ion
  long*  Z  =  (long*)calloc(number_of_particles, sizeof(long));
  long*  A  =  (long*)calloc(number_of_particles, sizeof(long));
  AT_Z_from_particle_no(        number_of_particles,
                                particle_no,
                                Z);
  AT_A_from_particle_no(        number_of_particles,
                                particle_no,
                                A);
  long i = 0;
  for (i = 0; i < number_of_particles; i++){
    if (particle_no[i] != 1){
      CSDA_range_g_cm2[i]  *=   (double)(A[i])/(double)((Z[i])*(Z[i]));
    }
  }

  free(Z);
  free(A);
}


void AT_CSDA_range_m(  const long  number_of_particles,
    const double  E_MeV_u[],
    const long   particle_no[],
    const long   material_no,
    double        CSDA_range_m[])
{
  // Get material density
  double material_density_g_cm3 = AT_density_g_cm3_from_material_no(material_no);

  // Get mass-norm. CSDA range
  AT_CSDA_range_g_cm2(  number_of_particles,
      E_MeV_u,
      particle_no,
      material_no,
      CSDA_range_m);

  long  i;
  for (i = 0; i < number_of_particles; i++){
    CSDA_range_m[i]  /=  material_density_g_cm3 * 100.0;
  }

}


void AT_E_MeV_from_CDSA_range(  const long  number_of_particles,
    const double  CSDA_range_g_cm2[],
    const long   particle_no[],
    const long   material_no,
    double        E_MeV[])
{
  // scaled energies
  double*  sE  =  (double*)calloc(number_of_particles, sizeof(double));

  get_table_value(number_of_particles, CSDA_range_g_cm2, material_no, AT_PSTAR_Data.range_cdsa_g_cm2, AT_PSTAR_Data.kin_E_MeV, sE);

  free( sE );
}


void AT_E_MeV_from_LET(  const long  number_of_particles,
    const double  LET_MeV_cm2_g[],
    const long   particle_no[],
    const long   material_no,
    double        E_MeV[])
{
  long  i;

  long*    Z             =  (long*)calloc(number_of_particles, sizeof(long));
  double*  Z_eff         =  (double*)calloc(number_of_particles, sizeof(double));
  long*    A             =  (long*)calloc(number_of_particles, sizeof(long));
  AT_Z_from_particle_no(        number_of_particles,
      particle_no,
      Z);
  AT_A_from_particle_no(        number_of_particles,
      particle_no,
      A);

  // 1. find approx. energy to compute effective charge
  double*       LET_MeV_cm2_g_copy = (double*)calloc(number_of_particles,sizeof(double));
  for(i = 0; i < number_of_particles; i++){
    LET_MeV_cm2_g_copy[i] = LET_MeV_cm2_g[i] / ((double)Z[i] * (double)Z[i] );
  }

  // 2. Do 10 rounds of approximation
  long j;
  for (j = 0; j < 10; j++){
    get_table_value(        number_of_particles,
                          LET_MeV_cm2_g_copy,
                          material_no,
                          AT_PSTAR_Data.stp_pow_el_MeV_cm2_g,
                          AT_PSTAR_Data.kin_E_MeV,
                          E_MeV);

    // 2a. Compute effective charge from approx. energies
    AT_effective_charge_from_E_MeV_u(  number_of_particles,
        E_MeV,
        particle_no,
        Z_eff);

    // 2b. find more accurate energy by using the computed effective charge
    for(i = 0; i < number_of_particles; i++){
      LET_MeV_cm2_g_copy[i] = LET_MeV_cm2_g[i] / ((double)Z_eff[i] * (double)Z_eff[i] );
    }
  }

  // 3. find eventual energy
  get_table_value(        number_of_particles,
                        LET_MeV_cm2_g_copy,
                        material_no,
                        AT_PSTAR_Data.stp_pow_el_MeV_cm2_g,
                        AT_PSTAR_Data.kin_E_MeV,
                        E_MeV);


  free( Z );
  free( A );
  free( Z_eff );
  free( LET_MeV_cm2_g_copy );
}


/////////////////////////////////////////////////////////
/* TEST FUNCTIONS FOR NEW MATERIAL / LET DATA HANDLING */
long AT_new_LET_MeV_cm2_g(  const long  number_of_particles,
    const double        E_MeV_u[],
    const long          particle_no[],
    AT_single_material_data_struct         material,
    double              LET_MeV_cm2_g[]){

    /* check input */
    int error;
    int purpose_energy_range = (material.LET_data_source == PSTAR) ? AT_energy_range_for_PSTAR_data : AT_energy_range_for_PowerLaw_data;
    error = AT_check_energy_range_single_field(   number_of_particles,
                                E_MeV_u,
                                purpose_energy_range);
    if(error != AT_Success)     return error;


    error = AT_check_particle_no_single_field(   number_of_particles,
                                particle_no);
    if(error != AT_Success)     return error;

    /* establish material if not yet done */
    int material_return_code = AT_check_material(&material);
    if((material_return_code != AT_Success) & (material_return_code != AT_Material_Already_Established)){
      return (material_return_code);
    }

    /* Loop over single particle function */
    long i;
    for (i = 0; i < number_of_particles; i++){
      LET_MeV_cm2_g[i]  =       AT_new_LET_MeV_cm2_g_single(  E_MeV_u[i],
          particle_no[i],
          material);
    }

    /* if material was established here, release memory */
    if(material_return_code == AT_Success){
      AT_free_material(&material);
    }

    return(AT_Success);
}

double AT_new_LET_MeV_cm2_g_single(  const double        E_MeV_u,
    const long          particle_no,
    AT_single_material_data_struct         material){

    /* establish material if not yet done */
    int material_return_code = AT_check_material(&material);
    if((material_return_code != AT_Success) & (material_return_code != AT_Material_Already_Established)){
      return (material_return_code);
    }

    /* TODO: Check if dedicated LET data is available for chosen particle */
    double LET_MeV_cm2_g = get_table_value_new( E_MeV_u,
        material.LET_data.LET_data_single[0].n,
        material.LET_data.LET_data_single[0].kin_E_MeV,
        material.LET_data.LET_data_single[0].stp_pow_el_MeV_cm2_g);

    /* If no dedicated table: get LET for proton of same energy / nucleon and scale by effective Z */
    double Zeff_ion    =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u, particle_no);
    double Zeff_proton =  AT_effective_charge_from_E_MeV_u_single(E_MeV_u, PARTICLE_PROTON_NUMBER);
    // scale proton LET by ratio of effective Z
    if( particle_no != PARTICLE_PROTON_NUMBER){ // for particles other than proton scale LET by (Zeff_ion / Zeff_proton)^2
        LET_MeV_cm2_g *=   gsl_pow_2(Zeff_ion / Zeff_proton);
    }

    /* if material was established here, release memory */
    if(material_return_code == AT_Success){
      AT_free_material(&material);
    }

    return LET_MeV_cm2_g;
}

double AT_CDSA_range_g_cm2_from_power_law_single(  const double E_MeV_u,
     const long particle_no,
     const double p_MeV,
     const double alpha_g_cm2_MeV)
{
  double Z_eff          = AT_effective_charge_from_E_MeV_u_single(  E_MeV_u,
      particle_no);
  double A              = (double)AT_A_from_particle_no_single( particle_no);
  return (A / (Z_eff * Z_eff) * alpha_g_cm2_MeV * pow(E_MeV_u, p_MeV));
}


double AT_LET_MeV_cm2_g_from_power_law_single(  const double E_MeV_u,
     const long particle_no,
     const double p_MeV,
     const double alpha_g_cm2_MeV)
{
  /* Get CDSA range for energy first */
  double CDSA_range_g_cm2 = AT_CDSA_range_g_cm2_from_power_law_single( E_MeV_u,
      particle_no,
      p_MeV,
      alpha_g_cm2_MeV);

  return 1.0 / (p_MeV * pow(alpha_g_cm2_MeV, 1.0 / p_MeV)) * pow(CDSA_range_g_cm2, 1.0 / p_MeV - 1.0);
}

#define POWERLAW_ENERGY_STEPS   100
#define POWERLAW_ENERGY_MIN_MEV 1
#define POWERLAW_ENERGY_MAX_MEV 250


int AT_establish_LET_data( AT_single_material_data_struct* material){

  if (material->LET_data_source == PowerLaw){
    material->LET_data.n                 = 1;
    material->LET_data.LET_data_single   = (AT_LET_data_single*)malloc(sizeof(AT_LET_data_single));

    material->LET_data.LET_data_single[0].n                      =       POWERLAW_ENERGY_STEPS;
    material->LET_data.LET_data_single[0].particle_no            =       1001;                           // build proton table
    material->LET_data.LET_data_single[0].kin_E_MeV              =       (double*)malloc(material->LET_data.LET_data_single[0].n * sizeof(double));
    material->LET_data.LET_data_single[0].stp_pow_el_MeV_cm2_g   =       (double*)malloc(material->LET_data.LET_data_single[0].n * sizeof(double));
    material->LET_data.LET_data_single[0].range_cdsa_g_cm2       =       (double*)malloc(material->LET_data.LET_data_single[0].n * sizeof(double));

    long i;
    double E_step_size = (log(POWERLAW_ENERGY_MAX_MEV) - log(POWERLAW_ENERGY_MIN_MEV)) / (POWERLAW_ENERGY_STEPS - 1);
    for (i = 0; i < material->LET_data.LET_data_single[0].n; i++){
      material->LET_data.LET_data_single[0].kin_E_MeV[i]                 = POWERLAW_ENERGY_MIN_MEV * exp(i * E_step_size);
      material->LET_data.LET_data_single[0].stp_pow_el_MeV_cm2_g[i]      = AT_LET_MeV_cm2_g_from_power_law_single( material->LET_data.LET_data_single[0].kin_E_MeV[i],
                                                                                                                  material->LET_data.LET_data_single[0].particle_no,
                                                                                                                  material->p_MeV,
                                                                                                                  material->alpha_g_cm2_MeV);
      material->LET_data.LET_data_single[0].range_cdsa_g_cm2[i]          = AT_CDSA_range_g_cm2_from_power_law_single(    material->LET_data.LET_data_single[0].kin_E_MeV[i],
                                                                                                                        material->LET_data.LET_data_single[0].particle_no,
                                                                                                                        material->p_MeV,
                                                                                                                        material->alpha_g_cm2_MeV);

    }
    return AT_Success;
  }

    if (material->LET_data_source == PSTAR){
       /* TODO: In case of non-predefined material: Read-In data from external PSTAR file */
       if (material->material_no == User_Defined_Material){
          return AT_No_PSTAR_Data;
        }

        /* find pre-defined material by AT_material_no */
        bool* matches           = (bool*)malloc(AT_PSTAR_Data.n * sizeof(bool));
        is_element_int( material->material_no,
                        AT_PSTAR_Data.material_no,
                        AT_PSTAR_Data.n,
                        matches);

        long i;
        long n_matches = 0;
        for (i = 0; i < AT_PSTAR_Data.n; i++){
          if (matches[i]){
            n_matches++;
          }
        }

        /* alloc data arrays and copy data */
        material->LET_data.n                 = 1;
        material->LET_data.LET_data_single   = (AT_LET_data_single*)malloc(sizeof(AT_LET_data_single));

        material->LET_data.LET_data_single[0].n                      =       n_matches;
        material->LET_data.LET_data_single[0].particle_no            =       1001;                           // build proton table
        material->LET_data.LET_data_single[0].kin_E_MeV              =       (double*)malloc(material->LET_data.LET_data_single[0].n * sizeof(double));
        material->LET_data.LET_data_single[0].stp_pow_el_MeV_cm2_g   =       (double*)malloc(material->LET_data.LET_data_single[0].n * sizeof(double));
        material->LET_data.LET_data_single[0].range_cdsa_g_cm2       =       (double*)malloc(material->LET_data.LET_data_single[0].n * sizeof(double));

        long j = 0;
        for (i = 0; i < material->LET_data.LET_data_single[0].n; i++){
          if(matches[i]){
            material->LET_data.LET_data_single[0].kin_E_MeV[j]             =       AT_PSTAR_Data.kin_E_MeV[i];
            material->LET_data.LET_data_single[0].stp_pow_el_MeV_cm2_g[j]  =       AT_PSTAR_Data.stp_pow_el_MeV_cm2_g[i];
            material->LET_data.LET_data_single[0].range_cdsa_g_cm2[j]      =       AT_PSTAR_Data.range_cdsa_g_cm2[i];
            j++;
          }
        }
        free(matches);

        return (AT_Success);
    }
    return AT_Unknown_LET_Data_Source;
}


double get_table_value_new( const double  x,
    const long    n,
    const double  x_table[],
    const double  y_table[])
{
  double  y = 0.0;
  double  err_y_tmp  = 0.0;    // dummy
  // Get tablulated value using 4th degree polynomial (n_pol - 1 = 2) interpolation
  long  n_pol      = 4 + 1;
  interp(    x_table,
      y_table,
      n,
      n_pol,
      x,
      &y,
      &err_y_tmp);

  return y;
}

/* END OF TEST FUNCTIONS FOR NEW MATERIAL / LET DATA HANDLING */
////////////////////////////////////////////////////////////////

