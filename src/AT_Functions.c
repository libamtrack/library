/**
*    AT_FileOperations.c
*    ===================
*
*    Created on: 28.07.2009
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

#include "AT_Functions.h"

///////////////////////////////////////////////////////////////////////
// Routines to access PARTICLE data
///////////////////////////////////////////////////////////////////////
void AT_Particle_Properties(  long*  particle_no,
    /* return values*/
    char**  particle_name,
    char**  USRTRACK_name,
    char**  element_name,
    long*  Z,
    long*  A,
    float*  mass)
{
  long i = (*particle_no) - 1;

  strcpy(*particle_name, AT_Particle_Data.particle_name[i]);
  strcpy(*USRTRACK_name, AT_Particle_Data.USRTRACK_name[i]);
  strcpy(*element_name,  AT_Particle_Data.element_name[i]);
  *Z    = AT_Particle_Data.Z[i];
  *A    = AT_Particle_Data.A[i];
  *mass = AT_Particle_Data.mass[i];
}

///////////////////////////////////////////////////////////////////////
// Routines to access MATERIAL data
///////////////////////////////////////////////////////////////////////
void AT_getMaterialData(    long*  n,
    long*  material_no,
    float*  density_g_cm3,
    float*  electron_density_m3,
    float*  I_eV,
    float*  alpha_g_cm2_MeV,
    float*  p_MeV,
    float*  m_g_cm2)
{
  long*  match  =  (long*)calloc(*n, sizeof(long));
  pmatchi(  material_no,
      n,
      AT_Material_Data.material_no,
      &AT_Material_Data.n,
      match);

  long i;
  for(i = 0; i < *n; i++){
    density_g_cm3[i]      = AT_Material_Data.density_g_cm3[match[i]];
    electron_density_m3[i]= AT_Material_Data.electron_density_m3[match[i]];
    I_eV[i]               = AT_Material_Data.I_eV[match[i]];
    alpha_g_cm2_MeV[i]    = AT_Material_Data.alpha_g_cm2_MeV[match[i]];
    p_MeV[i]              = AT_Material_Data.p_MeV[match[i]];
    m_g_cm2[i]            = AT_Material_Data.m_g_cm2[match[i]];
  }

  free(match);
}


#define matchIt      long  n_mat  = 1;                  \
		pmatchc(  &material_name,                \
				&n_mat,                    \
				AT_Material_Data.material_name,      \
				&AT_Material_Data.n,            \
				&match);

void AT_density_g_cm3( long*  n,
    char*  material_name,
    float*  density_g_cm3)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    density_g_cm3[i]    = AT_Material_Data.density_g_cm3[match];
  }
}

void AT_density_g_cm3S(    char**  material_name,
    float*  density_g_cm3){
  long  n;
  n    = 1;
  AT_density_g_cm3(  &n,
      *material_name,
      density_g_cm3);
}

void AT_electron_density_m3(  long*  n,
    char*  material_name,
    float*  electron_density_m3)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    electron_density_m3[i]  = AT_Material_Data.electron_density_m3[match];
  }
}

void AT_electron_density_m3S(  char**  material_name,
    float*  electron_density_m3){
  long  n;
  n    = 1;
  AT_electron_density_m3(  &n,
      *material_name,
      electron_density_m3);
}

void AT_alpha_g_cm2_MeV(    long*  n,
    char*  material_name,
    float*  alpha_g_cm2_MeV)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    alpha_g_cm2_MeV[i]    = AT_Material_Data.alpha_g_cm2_MeV[match];
  }
}

void AT_p_MeV(          long*  n,
    char*  material_name,
    float*  p_MeV)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    p_MeV[i]        = AT_Material_Data.p_MeV[match];
  }
}

void AT_m_g_cm2(        long*  n,
    char*  material_name,
    float*  m_g_cm2)
{
  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    m_g_cm2[i]        = AT_Material_Data.m_g_cm2[match];
  }
}



///////////////////////////////////////////////////////////////////////
// Routines to access PSTAR data
///////////////////////////////////////////////////////////////////////

void AT_LET_MeV_cm2_g(  long*  n,
    float*  E_MeV_u,
    long*  particle_no,
    long*  material_no,
    float*  LET_MeV_cm2_g)
{

#ifdef _DEBUG
  indnt_init();
  indnt_inc();
  fprintf(debf,"%sbegin AT_LET_MeV_cm2_g\n",isp);

  fprintf(debf,"%sn = %ld, material_no = %ld\n", isp, *n, *material_no);
  long ii;
  for( ii = 0 ; ii < *n ; ii++){
    fprintf(debf,"%sE_MeV_u[%ld]=%e\n", isp, ii , E_MeV_u[ii]);
    fprintf(debf,"%sparticle_no[%ld]=%ld\n", isp, ii , particle_no[ii]);
  }
#endif

  // get scaled energies for all given particles and energies
  float*  sE  =  (float*)calloc(*n, sizeof(float));
  AT_scaled_energy(  n,
      E_MeV_u,
      particle_no,
      sE);

#ifdef _DEBUG
  fprintf(debf,"%sE[0]=%e\n", isp, sE[0]);
#endif

  // get effective charge for all given particles and energies
  float*  eC  =  (float*)calloc(*n, sizeof(float));
  AT_effective_charge_from_particle_no(  n,
      E_MeV_u,
      particle_no,
      eC);

  getPSTARvalue(n, sE, material_no, AT_PSTAR_Data.kin_E_MeV, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, LET_MeV_cm2_g);
  long   i;
  for (i = 0; i < *n; i++){
    LET_MeV_cm2_g[i] *=   eC[i] * eC[i];
  }

  free(eC);
  free(sE);

#ifdef _DEBUG
  fprintf(debf,"%send AT_LET_MeV_cm2_g\n", isp);
  indnt_dec();
#endif

}


void AT_LET_keV_um(  long*  n,
    float*  E_MeV_u,
    long*  particle_no,
    long*  material_no,
    float*  LET_keV_um)
{
  long  match, n_tmp = 1;
  pmatchi(  material_no,
      &n_tmp,
      AT_Material_Data.material_no,
      &AT_Material_Data.n,
      &match);

  // Get mass-norm. LET
  AT_LET_MeV_cm2_g(  n,
      E_MeV_u,
      particle_no,
      material_no,
      LET_keV_um);

  long  i;
  for (i = 0; i < *n; i++){
    LET_keV_um[i]  *=  AT_Material_Data.density_g_cm3[match] * 0.1f;
  }

}

void AT_CSDA_range_g_cm2(  long*  n,
    float*  E_MeV_u,
    long*  particle_no,
    long*  material_no,
    float*  CSDA_range_g_cm2)
{
  getPSTARvalue(n, E_MeV_u, material_no, AT_PSTAR_Data.range_cdsa_g_cm2, AT_PSTAR_Data.kin_E_MeV, CSDA_range_g_cm2);
}

void AT_E_MeV_from_CDSA_range(  long*  n,
    float*  CSDA_range_g_cm2,
    long*  particle_no,
    long*  material_no,
    float*  E_MeV)
{
  // scaled energies
  float*  sE  =  (float*)calloc(*n, sizeof(float));

  getPSTARvalue(n, CSDA_range_g_cm2, material_no, AT_PSTAR_Data.range_cdsa_g_cm2, AT_PSTAR_Data.kin_E_MeV, sE);

  // scale energy back
  AT_E_MeV_u_from_scaled_energy(n , sE, particle_no, E_MeV);

  free( sE );
}

void AT_E_MeV_from_LET(  long*  n,
    float*  LET_MeV_cm2_g,
    long*  particle_no,
    long*  material_no,
    float*  E_MeV)
{

#ifdef _R
  int n_int = (int)(*n);
  *n = (long)n_int;

  int material_no_int  = (int)(*material_no);
  *material_no = (long)material_no_int;

  int particle_no_int = (int)(*particle_no);
  *particle_no = (long)particle_no_int;
#endif

  //printf("n = %ld\n", *n);
  //printf("particle_no = %ld\n", *particle_no);
  //printf("material_no = %ld\n", *material_no);
  //printf("AT_E_MeV_from_LET LET = %g\n", *LET_MeV_cm2_g);

  // scaled energies
  float*  sE  =  (float*)calloc(*n, sizeof(float));

  //TODO add effective charge correction !!

  // Conversion LETion => LETproton
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  float*  charge  =  (float*)calloc(*n, sizeof(float));
  pmatchi(  particle_no,
      n,
      AT_Particle_Data.particle_no,
      &AT_Particle_Data.n,
      matches);

  //printf("AT_E_MeV_from_LET LET ion = %g\n", LET_MeV_cm2_g[0]);

  // loop over n to find charge for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    charge[i]  = AT_Particle_Data.Z[matches[i]];
    LET_MeV_cm2_g[i] /= (charge[i]*charge[i]);
  }

  //printf("AT_E_MeV_from_LET LET proton = %g\n", LET_MeV_cm2_g[0]);

  getPSTARvalue(n, LET_MeV_cm2_g, material_no, AT_PSTAR_Data.stp_pow_el_MeV_cm2_g, AT_PSTAR_Data.kin_E_MeV, sE);


  // scale energy back
  AT_E_MeV_u_from_scaled_energy(n , sE, particle_no, E_MeV);

  free( sE );
}

