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

  strcpy(*particle_name,     AT_Particle_Data.particle_name[i]);
  strcpy(*USRTRACK_name,     AT_Particle_Data.USRTRACK_name[i]);
  strcpy(*element_name,    AT_Particle_Data.element_name[i]);
  *Z              = AT_Particle_Data.Z[i];
  *A              = AT_Particle_Data.A[i];
  *mass            = AT_Particle_Data.mass[i];
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
    density_g_cm3[i]    = AT_Material_Data.density_g_cm3[match[i]];
    electron_density_m3[i]  = AT_Material_Data.electron_density_m3[match[i]];
    I_eV[i]          = AT_Material_Data.I_eV[match[i]];
    alpha_g_cm2_MeV[i]    = AT_Material_Data.alpha_g_cm2_MeV[match[i]];
    p_MeV[i]        = AT_Material_Data.p_MeV[match[i]];
    m_g_cm2[i]        = AT_Material_Data.m_g_cm2[match[i]];
  }

  free(match);
}


#define matchIt      long  n_mat  = 1;                  \
            pmatchc(  &material_name,                \
                  &n_mat,                    \
                  AT_Material_Data.material_name,      \
                  &AT_Material_Data.n,            \
                  &match);

void AT_density_g_cm3(      long*  n,
                char*  material_name,
                float*  density_g_cm3)
{  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    density_g_cm3[i]    = AT_Material_Data.density_g_cm3[match];}}

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
{  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    electron_density_m3[i]  = AT_Material_Data.electron_density_m3[match];}}

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
{  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    alpha_g_cm2_MeV[i]    = AT_Material_Data.alpha_g_cm2_MeV[match];}}

void AT_p_MeV(          long*  n,
                char*  material_name,
                float*  p_MeV)
{  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    p_MeV[i]        = AT_Material_Data.p_MeV[match];}}

void AT_m_g_cm2(        long*  n,
                char*  material_name,
                float*  m_g_cm2)
{  long  match; matchIt;
  long  i;
  for(i = 0; i < *n; i++){
    m_g_cm2[i]        = AT_Material_Data.m_g_cm2[match];}}



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


//////////////////////////////////////////////////
// MISC ROUTINES
//////////////////////////////////////////////////

void AT_beta_from_particle_no(  long*  n,
                  float*  E_MeV_u,
                  long*  particle_no,
                  float*  beta)
{
  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  float*  mass  =  (float*)calloc(*n, sizeof(float));

  pmatchi(  particle_no,
        n,
        AT_Particle_Data.particle_no,
        &AT_Particle_Data.n,
        matches);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    mass[i]  = AT_Particle_Data.mass[matches[i]];}

  AT_beta_from_mass(  n,
            E_MeV_u,
            mass,
            beta);
  free(mass);
  free(matches);
}

void AT_beta_from_mass(  long*  n,
              float*  E_MeV_u,
              float*  mass,
              float*  beta)
{
  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){

    float  E_MeV    =  E_MeV_u[i] * mass[i];

    // Get rest energy
    float  E0_MeV    =  (float)proton_mass_MeV_c2 * mass[i];

    // Return relativistic speed
    beta[i]        =  (float)sqrt(1 - 1/((1 + E_MeV/E0_MeV)*(1 + E_MeV/E0_MeV)));
  }
}

void AT_E_from_beta_and_mass(  long*  n,
                float*  beta,
                float*  mass,
                float*  E_MeV_u)
{
  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){


    // Get rest energy
    float  E0_MeV    =  (float)proton_mass_MeV_c2 * mass[i];

    E_MeV_u[i]      =  E0_MeV * (1.0f / (1 - beta[i]*beta[i]) - 1);

    E_MeV_u[i]      /=  mass[i];
  }
}

void AT_E_from_beta_and_particle_no(  long*  n,
                      float*  beta,
                      long*  particle_no,
                      float*  E_MeV_u)
{
  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  float*  mass  =  (float*)calloc(*n, sizeof(float));

  pmatchi(  particle_no,
        n,
        AT_Particle_Data.particle_no,
        &AT_Particle_Data.n,
        matches);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    mass[i]  = AT_Particle_Data.mass[matches[i]];}

  AT_E_from_beta_and_mass(  n,
                beta,
                mass,
                E_MeV_u);
  free(mass);
  free(matches);
}


void AT_effective_charge_from_particle_no(  long*  n,
                        float*  E_MeV_u,
                        long*  particle_no,
                        float*  effective_charge)
{
  // get relativistic speeds for all given particles and energies
  float*  beta  =  (float*)calloc(*n, sizeof(float));
  long*  Z    =  (long*)calloc(*n, sizeof(long));

  AT_beta_from_particle_no(  n,
                  E_MeV_u,
                  particle_no,
                  beta);

  // find look-up indices for Z's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  pmatchi(  particle_no,
        n,
        AT_Particle_Data.particle_no,
        &AT_Particle_Data.n,
        matches);

  long i;
  for (i = 0; i < *n; i++){
    Z[i]  =  AT_Particle_Data.Z[matches[i]];
  }

  AT_effective_charge_from_beta(  n,
                  beta,
                  Z,
                  effective_charge);

  free(beta);
  free(Z);
  free(matches);
}

void AT_effective_charge_from_beta(  long*  n,
                    float*  beta,
                    long*  Z,
                    float*  effective_charge)
{
  // loop over n
  long  i;
  for (i = 0; i < *n; i++){
    // Return effective charge according to Barkas-Bethe-approximation (but not for protons!)
//    if (Z[i]!=1){
      effective_charge[i]  = (float)(Z[i]) * (1 - (float)exp(-125.0f * beta[i] / (pow(Z[i], 2.0f/3.0f))));//}
//    else{
//      effective_charge[i]  = (float)(Z[i]);}
  }
}



void AT_scaled_energy(  long*  n,
            float*  E_MeV_u,
            long*  particle_no,
            float*  scaled_energy)
{
//  printf("begin AT_scaled_energy\n");
//  printf("n = %ld\n", *n);
//  long ii;
//  for( ii = 0 ; ii < *n ; ii++){
//    printf("E_MeV_u[%ld]=%e\n", ii , E_MeV_u[ii]);
//    printf("particle_no[%ld]=%ld\n", ii , particle_no[ii]);
//  }

  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  pmatchi(  particle_no,
        n,
        AT_Particle_Data.particle_no,
        &AT_Particle_Data.n,
        matches);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){

    //printf("matches[0] = %ld\n", matches[i]);

    // total kinetic energy
    float  E_MeV    =  E_MeV_u[i] * AT_Particle_Data.A[matches[i]];

    float mass_fraction = AT_Particle_Data.mass[matches[i]] / AT_Particle_Data.mass[0];

    //printf("mass_fraction = %g\n", mass_fraction );

    // Return mass-scaled energy
    scaled_energy[i]  =  E_MeV / mass_fraction ;
  }

  free(matches);

  //printf("end AT_scaled_energy\n");

}

void AT_E_MeV_u_from_scaled_energy(  long*  n,
            float*  scaled_energy,
            long*  particle_no,
            float*  E_MeV_u)
{
  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));
  pmatchi(  particle_no,
        n,
        AT_Particle_Data.particle_no,
        &AT_Particle_Data.n,
        matches);

  // loop over n to find beta for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){

    //printf("matches[0] = %ld\n", matches[i]);

    float mass_fraction = AT_Particle_Data.mass[matches[i]] / AT_Particle_Data.mass[0];

    float  E_MeV  = scaled_energy[i] * mass_fraction;

    // Return energy per nucleon
    E_MeV_u[i]  =  E_MeV / AT_Particle_Data.A[matches[i]];

    printf(" AT_E_MeV_u_from_scaled_energy E_MeV_u = %f\n", E_MeV_u[i]);

  }

  free(matches);

}


void AT_max_E_transfer_MeV(  long*  n,
                float*  E_MeV_u,
                long*  particle_no,
                float*  max_E_transfer_MeV)
{
  /* if E_MeV_u < 0:    use non-relativistic formula
   * if E_MeV_u > 0:    use relativistic formula
   */
  int*  relativistic    =  (int*)calloc(*n, sizeof(int));
  long  i;
  for (i = 0; i < *n; i++){
    if(E_MeV_u[i] >= 0){
      relativistic[i]      = 1;
    }else{
      relativistic[i]      = 0;
      E_MeV_u[i]        *= -1.0f;
    }
  }

  // get relativistic speeds for all given particles and energies
  float*  beta  =  (float*)calloc(*n, sizeof(float));
  AT_beta_from_particle_no(  n,
                  E_MeV_u,
                  particle_no,
                  beta);

  for (i = 0; i < *n; i++){
    if(relativistic[i] == 0){
      max_E_transfer_MeV[i]  =  4.0f * electron_mass_MeV_c2 / proton_mass_MeV_c2* E_MeV_u[i];
    }else{
      max_E_transfer_MeV[i]  =  2.0f * electron_mass_MeV_c2 * beta[i] * beta[i] / (1.0f - beta[i] * beta[i]);
    }
  }

  free(relativistic);
  free(beta);
}

// ELECTRON RANGE
#ifdef _R
void AT_max_electron_range_mS(  int*  n,
                float*  E_MeV_u,
                int*  particle_no,
                int*  material_no,
                int*   er_model,
                float*  max_electron_range_m)
{
  long n_long = (long)(*n);
  long material_no_long = (long)(*material_no);
  long er_model_long = (long)(*er_model);

  long i;
  long * particle_no_long = (long*)calloc(*n,sizeof(long));
  for(i = 0 ; i < *n ; i++){
    particle_no_long[i] = (long)particle_no[i];
//    printf("output particle_no[%ld]=%ld\n", i , particle_no_long[i]);
  }

  AT_max_electron_range_m( &n_long,
                E_MeV_u,
                particle_no_long,
                &material_no_long,
                &er_model_long,
                max_electron_range_m);

  free(particle_no_long);

}
#endif

#ifdef _S
void AT_max_electron_range_mS(  long*  n,
                float*  E_MeV_u,
                long*  particle_no,
                long*  material_no,
                long*   er_model,
                float*  max_electron_range_m)
{

  AT_max_electron_range_m( n,
                E_MeV_u,
                particle_no,
                material_no,
                er_model,
                max_electron_range_m);

}
#endif

void AT_max_electron_range_m(  long*  n,
                float*  E_MeV_u,
                long*  particle_no,
                long*  material_no,
                long*   er_model,
                float*  max_electron_range_m)
{
#ifdef _DEBUG
  indnt_init();
  indnt_inc();
  fprintf(debf,"%sbegin AT_max_electron_range_m\n",isp);
  fprintf(debf,"%sn = %ld, material_no = %ld, er_model = %ld\n", isp, *n, *material_no, *er_model);
  long ii;
  for( ii = 0 ; ii < *n ; ii++){
    fprintf(debf,"%sE_MeV_u[%ld]=%e\n", isp, ii , E_MeV_u[ii]);
    fprintf(debf,"%sparticle_no[%ld]=%ld\n", isp, ii , particle_no[ii]);
  }
#endif

  // Get density matching to material_name (only 1 name therefore n_mat = 1)
  //long  n_mat  = 1;
  //long  match;
  //pmatchi(  material_no,
  //      &n_mat,
  //      AT_Material_Data.material_no,
  //      &AT_Material_Data.n,
  //      &match);

  long*  matches  =  (long*)calloc(*n, sizeof(long));
  float*  mass  =  (float*)calloc(*n, sizeof(float));

  pmatchi(  particle_no,
        n,
        AT_Particle_Data.particle_no,
        &AT_Particle_Data.n,
        matches);

  long  i;
  for (i = 0; i < *n; i++){
    float tmpE  = E_MeV_u[i];
    if (tmpE < 0) {tmpE *= -1.0f;}  // E can be set neg. if non-PSTAR are given --> use pos. value

    mass[i]  = AT_Particle_Data.mass[matches[i]];
    float E_div_E0 = E_MeV_u[i] / (mass[i]*proton_mass_MeV_c2);
    float w_keV;
    if( *er_model == ER_ButtsKatz ){
      w_keV = 2 * electron_mass_MeV_c2 * ( E_div_E0*E_div_E0 + 2*E_div_E0) * 1e3;
      max_electron_range_m[i] = 1e-5 * w_keV;
    }
    if( *er_model == ER_Waligorski ){
      double alpha = 1.667;
      w_keV = 2 * electron_mass_MeV_c2 * ( E_div_E0*E_div_E0 + 2*E_div_E0) * 1e3;
      if( w_keV < 1. ) alpha = 1.079;
      max_electron_range_m[i] =  6* 1e-6 * (float)pow( w_keV, alpha );
    }
    if( *er_model == ER_Geiss ){
      max_electron_range_m[i] = 4e-5 * (float)pow(tmpE, 1.5);
    }
    if( *er_model == ER_Scholz ){
      max_electron_range_m[i] = 5e-5 * (float)pow(tmpE, 1.7);
    }

    // Scale maximum el. range with material density relative to water (1/rho)
    max_electron_range_m[i]    /= AT_Material_Data.density_g_cm3[*material_no];

    // covert cm to m
    max_electron_range_m[i]    /= 1e2;  // cm to m

  }
#ifdef _DEBUG
  for( ii = 0 ; ii < *n ; ii++){
    fprintf(debf,"%srange[%ld]=%e\n", isp, ii , max_electron_range_m[ii]);
  }
  fprintf(debf,"%send AT_max_electron_range_m\n",isp);
  indnt_dec();
#endif
  free(matches);
  free(mass);


}

void AT_Bohr_Energy_Straggling_g_cm2(  long*  n,
                    char**  material_name,
                    float*  dsE2dz)
{
  // Get Bohr's energy spread (Wilson, 1947, Phys Rev 71, 385)
  long  i;
  for (i = 0; i < *n; i++){
    float  electron_density_m3;
    long  n_dummy = 1;
    AT_electron_density_m3(  &n_dummy,
                  material_name[i],
                  &electron_density_m3);

    double  tmp          =  e_C * e_C * e_C * e_C * electron_density_m3;
    tmp              /=  4.0 * pi * e0_F_m * e0_F_m;
    tmp              /=  MeV_to_J * MeV_to_J * m_to_cm;
    dsE2dz[i]          =  (float)tmp;
  }
}

void AT_Z_from_particle_no(  long*  n,
                long*  particle_no,
                long*  Z)
{
  // find look-up indices for A's for particle numbers in particle data
  long*  matches  =  (long*)calloc(*n, sizeof(long));

  pmatchi(  particle_no,
        n,
        AT_Particle_Data.particle_no,
        &AT_Particle_Data.n,
        matches);

  // loop over n to find Z for all given particles and energies
  long  i;
  for(i = 0; i < *n; i++){
    Z[i]  = AT_Particle_Data.Z[matches[i]];}

  free(matches);
}
