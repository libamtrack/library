/**
 * @file
 * @brief Radial Dose Distribution models
 */

/*
 *    AT_RDD.c
 *    ========
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

#include "AT_RDD.h"


void getRDDName(const long* RDD_no, char* RDD_name){
  strcpy(RDD_name,"*** invalid choice ***");
  long i;
  for (i = 0; i < RDD_DATA_N; i++){
    if (AT_RDD_Data.RDD_no[i] == *RDD_no){
      strcpy(RDD_name, AT_RDD_Data.RDD_name[i]);
    }
  }
}

void getRDDNo(const char* RDD_name, long* RDD_no){
  *RDD_no = 0;
  long i;
  for (i = 0; i < RDD_DATA_N; i++){
    if (strcmp(RDD_name, AT_RDD_Data.RDD_name[i]) == 0){
      *RDD_no = AT_RDD_Data.RDD_no[i];
      break;
    }
  }
}


inline float AT_RDD_Katz_point_kernel(const float* x, const float* alpha){
  return (1.0f/((*x)*(*x)) )*pow(1.0f - (*x), 1.0f / (*alpha));
}

void AT_RDD_Katz_point_kernelS(int *n, float* x, float* alpha, float* f){
  int i;
  for( i = 0 ; i < *n ; i++){
    f[i] = AT_RDD_Katz_point_kernel(&(x[i]),alpha);
  };
}

inline float AT_RDD_Katz_point_coeff_Gy(const float* C_J_m,const float* Z_eff, const float* beta, const float* alpha, const float* density_kg_m3, const float* r_max_m){
  return (*C_J_m) * (*Z_eff)*(*Z_eff) / (2.0f * M_PI * (*beta)*(*beta) * (*alpha) * (*density_kg_m3) * (*r_max_m)* (*r_max_m));
}

inline float AT_RDD_Katz_point_Gy(const float* r_m, const float* alpha, const float* r_max_m, const float* Katz_point_coeff_Gy){
  float x = (*r_m)/(*r_max_m);
  return (*Katz_point_coeff_Gy) * AT_RDD_Katz_point_kernel(&x,alpha);
}

void AT_RDD_Katz_point_GyS(int *n, float* r_m, float* alpha, float* r_max_m,float* Katz_point_coeff_Gy, float * D){
  int i;
  for( i = 0 ; i < *n ; i++){
    D[i] = AT_RDD_Katz_point_Gy(&(r_m[i]),alpha,r_max_m,Katz_point_coeff_Gy);
  };
}

inline float AT_RDD_Katz_dEdx_kernel(const float* x, const float* alpha){
  return (1.0f/(*x) )*pow(1.0f - (*x), 1.0f / (*alpha));
}

void AT_RDD_Katz_dEdx_kernelS(int *n, float* x, float* alpha, float* f){
  int i;
  for( i = 0 ; i < *n ; i++){
    f[i] = AT_RDD_Katz_dEdx_kernel(&(x[i]),alpha);
  };
}

double AT_RDD_Katz_dEdx_integrand(double x, void * params){
  float alpha = ((float*)params)[0];
  float f_x = (float)(x);
  return (double)AT_RDD_Katz_dEdx_kernel( &f_x, &alpha);
}

inline float AT_RDD_Katz_dEdx_coeff_J_m(const float* r_max_m, const float* density_kg_m3, const float* Katz_point_coeff_Gy){
  return 2 * M_PI * (*density_kg_m3) * (*r_max_m)*(*r_max_m) * (*Katz_point_coeff_Gy);
}

float AT_RDD_Katz_dEdx_J_m(  const float* alpha,
    const float* r_min_m,
    const float* r_max_m,
    const float* Katz_dEdx_coeff_J_m){
  double dEdx_integral = 0.0;
  if( (*r_min_m) < (*r_max_m)){
    dEdx_integral = (*alpha/(1.+(*alpha))) * pow( 1.-(*r_min_m)/(*r_max_m) , 1. +
        1./(*alpha) ) *
        gsl_sf_hyperg_2F1(1.,1.+1./(*alpha),2.+1./(*alpha),1.-(*r_min_m)/(*r_max_m));
  }
  return (*Katz_dEdx_coeff_J_m)*(float)dEdx_integral;
}

inline float AT_RDD_Katz_site_Gy(const float* r_m, const float* alpha, const float* r_min_m, const float* r_max_m, const float* LET_J_m, const float* density_kg_m3, const float* Katz_dEdx_J_m, const float* Katz_point_coeff_Gy){
  if( (*r_m) < (*r_min_m) ){
    return (1.0f / ((*density_kg_m3) * M_PI * (*r_min_m)*(*r_min_m)))*((*LET_J_m) - (*Katz_dEdx_J_m));
  } else {
    return AT_RDD_Katz_point_Gy(r_m, alpha, r_max_m, Katz_point_coeff_Gy);
  }
  return 0.0f;
}

void AT_RDD_Katz_site_GyS(int *n, float* r_m, float* alpha, float* r_min_m, float* r_max_m, float* LET_J_m, float* density_kg_m3, float* Katz_dEdx_J_m, float* Katz_point_coeff_Gy, float * D_Gy){
  int i;
  for( i = 0 ; i < *n ; i++){
    D_Gy[i] = AT_RDD_Katz_site_Gy(&(r_m[i]),alpha,r_min_m,r_max_m,LET_J_m,density_kg_m3,Katz_dEdx_J_m,Katz_point_coeff_Gy);
  };
}

float
geometryFunctionPhi (const float* r0_m, const float* a0_m, const float* r_m)
{
  float res = 0.;
  double factor = 0.;
  gsl_complex carg, cres;
  if ((*r_m) <= fabs ((*r0_m) - (*a0_m)))
  {
    if ((*r0_m) >= (*a0_m))
      res = 0.0f;
    else
      res = (float)M_PI;
  }
  else
  {
    factor = gsl_pow_2 (*a0_m) - gsl_pow_2 ((*r0_m) - (*r_m));
    factor /= gsl_pow_2 ((*r_m) + (*r0_m)) - gsl_pow_2 (*a0_m);
    GSL_SET_COMPLEX (&carg, sqrt (factor), 0.);
    cres = gsl_complex_arctan (carg);
    res = 2.0f * (float)(GSL_REAL (cres));
  }
  return res;
}

inline float AT_RDD_Katz_ext_kernel_Gy(const float* t_m, const float *r_m, const float* a0_m, const float* alpha, const float* r_min_m, const float* r_max_m, const float* Katz_point_coeff_Gy){
  if( (*t_m) < (*r_min_m ) )
    return 0.0;
  if( (*t_m) >= (*a0_m) + (*r_m) )
    return 0.0;
  if( ((*r_m) >= (*a0_m)) && ((*t_m) <= (*r_m) - (*a0_m)))
    return 0.0;
  else
    return  (1.0f/ (M_PI * (*a0_m)*(*a0_m))) * AT_RDD_Katz_point_Gy(t_m,alpha,r_max_m,Katz_point_coeff_Gy) *  geometryFunctionPhi(r_m,a0_m,t_m) * (*t_m);
}

void AT_RDD_Katz_ext_kernel_GyS(int *n, float* t_m, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy){
  int i;
  for( i = 0 ; i < *n ; i++){
    D_Gy[i] = AT_RDD_Katz_ext_kernel_Gy(  &(t_m[i]), r_m,a0_m, alpha, r_min_m, r_max_m, Katz_point_coeff_Gy);
  };
}

double AT_RDD_Katz_ext_integrand_Gy(  double t_m, void * params){
  float r_m           = ((float*)params)[0];
  float a0_m           = ((float*)params)[1];
  float alpha         = ((float*)params)[2];
  float r_min_m         = ((float*)params)[3];
  float r_max_m         = ((float*)params)[4];
  float Katz_point_coeff_Gy   = ((float*)params)[5];
  float f_t_m         = (float)(t_m);
  return (double)AT_RDD_Katz_ext_kernel_Gy(  &f_t_m, &r_m, &a0_m, &alpha,
      &r_min_m, &r_max_m, &Katz_point_coeff_Gy);
}

inline float AT_RDD_Katz_ext_Gy(  const float *r_m, const float* a0_m, const float* alpha,
    const float* r_min_m, const float* r_max_m, const float* Katz_point_coeff_Gy){
  double int_lim_m = 0.0f;
  if( (*r_m) > (*a0_m) ){
    int_lim_m = GSL_MAX((*r_m) - (*a0_m),*r_min_m);
  }
    gsl_set_error_handler_off();

  double ext_integral_Gy;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (10000);
  gsl_function F;
  F.function = &AT_RDD_Katz_ext_integrand_Gy;
  float params[] = {*r_m,*a0_m,*alpha,*r_min_m,*r_max_m,*Katz_point_coeff_Gy};
  F.params = params;
  int status = gsl_integration_qags (&F, int_lim_m, (*r_m)+(*a0_m), 1e-9, 1e-4, 10000, w1, &ext_integral_Gy, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
    ext_integral_Gy = -1.0f;
  }
  gsl_integration_workspace_free (w1);

  return ext_integral_Gy;
}

void AT_RDD_Katz_ext_GyS(int *n, float *r_m, float* a0_m, float* alpha, float* r_min_m, float* r_max_m, float* Katz_point_coeff_Gy, float * D_Gy){
  int i;
  for( i = 0 ; i < *n ; i++){
    D_Gy[i] = AT_RDD_Katz_ext_Gy(&(r_m[i]), a0_m, alpha, r_min_m, r_max_m, Katz_point_coeff_Gy);
  };
}


inline float AT_RDD_Cucinotta_C(float Z, float beta){
  return Z*beta;
}

void AT_RDD_f1_parameters(  /* radiation field parameters */
    const float*  E_MeV_u,
    const long*  particle_no,
    /* detector parameters */
    const long*  material_no,
    /* radial dose distribution model */
    const long*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const long*  er_model,
    const float*  er_parameter,
    /* calculated parameters */
    float * f1_parameters)
{
  long  n_tmp    = 1;
  float beta     = 0.0f;
  long  Z        = 0;
  float Z_eff    = 0.0f;

  float LET_MeV_cm2_g = 0.0f;
  float max_electron_range_m = 0.0f;
  float r_min_m = 0.0f;
  float single_impact_fluence_cm2 = 0.0f;
  float single_impact_dose_Gy = 0.0f;
  float norm_constant_Gy = 0.0f;
  float d_min_Gy = 0.0f;
  float d_max_Gy = 0.0f;
  float dEdx_MeV_cm2_g = 0.0f;

  AT_beta_from_particle_no(  &n_tmp,
                  E_MeV_u,
                  particle_no,
                  &beta);
  AT_Z_from_particle_no(  &n_tmp,
                particle_no,
                &Z);
  AT_effective_charge_from_beta(  &n_tmp,
                  &beta,
                  &Z,
                  &Z_eff);

  // Get density
  float   density_g_cm3, density_kg_m3, electron_density_m3;
  AT_getMaterialData(  &n_tmp,
              material_no,
              &density_g_cm3,
              &electron_density_m3,
              NULL, NULL, NULL, NULL, NULL, NULL);
  density_kg_m3      =  density_g_cm3 * 1000.0f;

  ///////////////////////////////////////////////////////////////
  // PARAMETER 0: Get the LET (same for all models)
  AT_LET_MeV_cm2_g(  &n_tmp,
            E_MeV_u,
            particle_no,
            material_no,
            &LET_MeV_cm2_g);

  /////////////////////////////////////////////////////////////////////////////
  // PARAMETER 2: Get the maximum electron range (same for all RDD models)
  AT_max_electron_range_m(  &n_tmp,
                E_MeV_u,
                particle_no,
                material_no,
                er_model,
                &max_electron_range_m);


  /////////////////////////////////////////////////////////////////////////////
  // PARAMETER 6: Get the single impact fluence (same for all RDD models)
  single_impact_fluence_cm2  = M_1_PI / gsl_pow_2( max_electron_range_m * m_to_cm ); // pi * r_max_m^2 = Track area -> single_impact_fluence [1/cm2]


  /////////////////////////////////////////////////////////////////////////////
  // MODEL SPECIFIC PARAMETERS
  if( *rdd_model == RDD_Test){
    r_min_m = 0.0f;                                                                             // r_min_m
    norm_constant_Gy = single_impact_fluence_cm2 * LET_MeV_cm2_g * MeV_to_J * 1000.0f;          // LET  / track area = Norm.constant k
    single_impact_dose_Gy = LET_MeV_cm2_g * f1_parameters[5] * MeV_g_to_J_kg;                   // single_impact_dose = LET * fluence;
    d_min_Gy = norm_constant_Gy;                                                                // d_min_Gy = k
    d_max_Gy = norm_constant_Gy;                                                                // d_max_Gy = k
    dEdx_MeV_cm2_g = LET_MeV_cm2_g;                                                             // dEdx = LET

    // Save parameters to f1_parameters table
    f1_parameters[1]  = r_min_m;
  }

  if( *rdd_model == RDD_KatzPoint || *rdd_model == RDD_Site || *rdd_model == RDD_ExtTarget || *rdd_model == RDD_Edmund){
    //////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////
    float alpha        =  1.667f;
    float wmax_MeV     =  0.0f;
    AT_max_E_transfer_MeV(&n_tmp,E_MeV_u,particle_no,&wmax_MeV);
    if(wmax_MeV <= 1e-3){  // if wmax < 1keV
      alpha           = 1.079f;
    }
    //////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////

    r_min_m          =  rdd_parameter[0];
    if (max_electron_range_m <= r_min_m){
      r_min_m = max_electron_range_m;
    }                  // If r.max < a0 or r.min, r.min = r.max, not a0

    float  N_el_cm3    =  electron_density_m3 / (1e6);
    float  C_J_cm      =  2.0f * M_PI * N_el_cm3 * gsl_pow_4(e_esu) / (electron_mass_g * gsl_pow_2(c_cm_s)) * 1e-7;   // energy constant [J/cm] not [erg/cm] hence 10^-7
    float  C_J_m       =  C_J_cm * 100.0f;
    norm_constant_Gy   =  C_J_m;                      // Norm.constant k


    // Get dEdx by simple integration from r_min_m (in case of RDD_Site = a0) to r_max_m
    float  dEdx_J_m    =  0.0f;
    float  LET_J_m     =  LET_MeV_cm2_g * 1.602e-13 * density_g_cm3 * 100.0f;

    float Katz_point_coeff_Gy = AT_RDD_Katz_point_coeff_Gy(&C_J_m,&Z_eff,&beta,&alpha,&density_kg_m3,&max_electron_range_m);
    float Katz_dEdx_coeff_J_m = AT_RDD_Katz_dEdx_coeff_J_m(&max_electron_range_m,&density_kg_m3,&Katz_point_coeff_Gy);
    dEdx_J_m                  = AT_RDD_Katz_dEdx_J_m(&alpha,&r_min_m,&max_electron_range_m,&Katz_dEdx_coeff_J_m);

    dEdx_MeV_cm2_g            =  dEdx_J_m / 100.0f / density_g_cm3 / MeV_to_J;

    if( *rdd_model != RDD_ExtTarget )
      d_min_Gy    =  rdd_parameter[1];
    else
      d_min_Gy    =  rdd_parameter[2];

    // d_max_Gy
    if(*rdd_model == RDD_Site){
      float tmp = 0.0f;
      d_max_Gy    =  AT_RDD_Katz_site_Gy(&tmp,&alpha,&r_min_m,&max_electron_range_m,&LET_J_m,&density_kg_m3,&Katz_dEdx_coeff_J_m,&Katz_point_coeff_Gy);
    }else if(  *rdd_model == RDD_Edmund){
      float tmp =  LET_MeV_cm2_g - dEdx_MeV_cm2_g;  // LET - dEdx (MeV_g_cm2)
      if(tmp > 0){
        tmp *=  density_g_cm3; // MeV / cm
        tmp *=  MeV_to_J;      // J / cm
        tmp *=  100.0f;        // J / m
        float core_kg_m  = M_PI * gsl_pow_2(rdd_parameter[0]) * density_kg_m3;  // kg / m
        d_max_Gy = tmp / core_kg_m;                                             // J / kg = Gy
      }else{
        d_max_Gy = 1e-11;
      }
    }else if(*rdd_model == RDD_KatzPoint){
      d_max_Gy    =  AT_RDD_Katz_point_Gy(&r_min_m,&alpha,&max_electron_range_m,&Katz_point_coeff_Gy);
    } else { // RDD_ExtTarget
      d_max_Gy    =  AT_RDD_Katz_ext_Gy(&r_min_m,&(rdd_parameter[1]),&alpha,&r_min_m,&max_electron_range_m,&Katz_point_coeff_Gy);
    }

    // single_impact_dose
    if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund){
      single_impact_dose_Gy    =  LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;        // LET * fluence
    }else{
      single_impact_dose_Gy    =  dEdx_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;       // dEdx * fluence
    }

    // Save parameters to f1_parameters table
    f1_parameters[1]  = r_min_m;
  }

  if( *rdd_model == RDD_Geiss){
    float a0_m                 = 0.0;
    a0_m                       =  rdd_parameter[0];                                          // "r_min_m" = a0
    if (max_electron_range_m <= a0_m){
      a0_m = max_electron_range_m;}      // If r.max < a0, r.min = r.max, not a0
    // Normalization to match with LET
    float  tmp                 = (float)(0.5f + log(max_electron_range_m / a0_m));
    tmp                       *= 2.0f * M_PI * gsl_pow_2(a0_m * m_to_cm);
    norm_constant_Gy           = LET_MeV_cm2_g * MeV_g_to_J_kg / tmp;                         // k = LET / tmp
    single_impact_dose_Gy      = LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;   // single_impact_dose = LET * single_impact_fluence
    d_max_Gy                   = norm_constant_Gy;                                            // d_max_Gy = k
    d_min_Gy                   = norm_constant_Gy * gsl_pow_2( a0_m / max_electron_range_m);  // d_min_Gy
    dEdx_MeV_cm2_g             = LET_MeV_cm2_g;                                               // dEdx = LET

    // Save parameters to f1_parameters table
    f1_parameters[1]  = a0_m;
  }

  if( *rdd_model == RDD_Cucinotta){ // TODO to be implemented
    r_min_m  = 0.0f;
    norm_constant_Gy  = 0.0f;  // k = ????
    single_impact_dose_Gy  = LET_MeV_cm2_g * MeV_g_to_J_kg * single_impact_fluence_cm2;  // single_impact_dose = LET * single_impact_fluence
    d_max_Gy  = 0.0f;                                      // d_max_Gy = k
    d_min_Gy  = 0.0f;
    dEdx_MeV_cm2_g  = LET_MeV_cm2_g;                       // dEdx = LET

    // Save parameters to f1_parameters table
    f1_parameters[1]  = r_min_m;
  }

  // write data to output table (apart from f1_parameters[0] which sometimes is
  // r_min_m and sometimes a0_m )
  f1_parameters[0] = LET_MeV_cm2_g;
  f1_parameters[2] = max_electron_range_m;
  f1_parameters[3]  = d_min_Gy;
  f1_parameters[4]  = d_max_Gy;
  f1_parameters[5]  = norm_constant_Gy;
  f1_parameters[6]  = single_impact_fluence_cm2;
  f1_parameters[7]  = single_impact_dose_Gy;
  f1_parameters[8]  = LET_MeV_cm2_g;

}

void AT_D_RDD_Gy  ( const  long*  n,
    const float*  r_m,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const long*  particle_no,
    /* detector parameters */
    const long*  material_no,
    /* radial dose distribution model */
    const long*  rdd_model,
    const float*  rdd_parameter,
    /* electron range model */
    const long*  er_model,
    const float*  er_parameter,
    float*  D_RDD_Gy)
{
  long     i;
  long    n_tmp      = 1;
  // Get f1 parameters
  long     n_f1_parameters = 9;
  float*   f1_parameters   = (float*)calloc(n_f1_parameters, sizeof(float));
  AT_RDD_f1_parameters(      /* radiation field parameters */
      E_MeV_u,
      particle_no,
      /* detector parameters */
      material_no,
      /* radial dose distribution model */
      rdd_model,
      rdd_parameter,
      /* electron range model */
      er_model,
      er_parameter,
      /* calculated parameters */
      f1_parameters);

  // TODO assign f1_parameters to some more human-readable names

  // Get material data
  float     density_g_cm3;
  AT_getMaterialData(  &n_tmp,
        material_no,
        &density_g_cm3,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  float    density_kg_m3  =  density_g_cm3 * 1000.0f;


  if( *rdd_model == RDD_Test){
    // Loop over all r_m given
    for (i = 0; i < *n; i++){
      D_RDD_Gy[i]    =  0.0f;
      if (r_m[i] < f1_parameters[2]){
        D_RDD_Gy[i]    =  f1_parameters[5];
      } else {
        D_RDD_Gy[i]    =  0;
      }
    }
  }// end RDD_Test


  if( *rdd_model == RDD_KatzPoint || *rdd_model == RDD_Site || *rdd_model == RDD_ExtTarget || *rdd_model == RDD_Edmund){
    // Get beta, Z and Zeff
    float  beta    = 0.0f;
    long   Z       = 0;
    float  Z_eff   = 0.0f;
    AT_beta_from_particle_no(  &n_tmp,
                  E_MeV_u,
                  particle_no,
                  &beta);
    AT_Z_from_particle_no(  &n_tmp,
                  particle_no,
                  &Z);
    AT_effective_charge_from_beta(  &n_tmp,
                  &beta,
                  &Z,
                  &Z_eff);

    //////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////
    float alpha        =  1.667f;
    float wmax_MeV     =  0.0f;
    AT_max_E_transfer_MeV(&n_tmp,E_MeV_u,particle_no,&wmax_MeV);
    if(wmax_MeV <= 1e-3){  // if wmax < 1keV
      alpha           = 1.079f;
    }
    //////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////

    // Loop over all r_m given
    float Katz_point_coeff_Gy = AT_RDD_Katz_point_coeff_Gy(&(f1_parameters[5]),&Z_eff,&beta,&alpha,&density_kg_m3,&(f1_parameters[2]));
    float a0         = rdd_parameter[0];
    for (i = 0; i < *n; i++){
      D_RDD_Gy[i]        =  0.0f;                        // r < r_min_m (for RDD_KatzPoint) or r > r_max_m --> D = 0
      if (r_m[i] >= f1_parameters[1] && r_m[i] <= f1_parameters[2] && *rdd_model != RDD_ExtTarget){          // in between r_min and r_max --> D = KatzPoint
        //D_RDD_Gy[i]        = f1_parameters[5] * Z_eff*Z_eff / (2.0f * pi * r_m[i]*r_m[i] * beta*beta * alpha * density_kg_m3) * pow(1.0f - r_m[i] / f1_parameters[2], 1.0f / alpha);
        D_RDD_Gy[i]        = AT_RDD_Katz_point_Gy(&(r_m[i]),&alpha,&(f1_parameters[2]),&Katz_point_coeff_Gy);
        D_RDD_Gy[i]        = fmaxf(D_RDD_Gy[i], f1_parameters[3]);          // Cut-off low doses
      }
      if (r_m[i] <= f1_parameters[1] && *rdd_model == RDD_Site){            // r < r_min_m (for RDD_Site) --> D = d_max_Gy
        D_RDD_Gy[i]        = f1_parameters[4];
      }
      if (r_m[i] < f1_parameters[1] && *rdd_model == RDD_Edmund){            // r < r_min_m (for RDD_Site) --> D = d_max_Gy
        D_RDD_Gy[i]        = f1_parameters[4];
      }
      if (r_m[i] <= f1_parameters[2] && *rdd_model == RDD_ExtTarget){
        D_RDD_Gy[i]        = AT_RDD_Katz_ext_Gy(&(r_m[i]),&a0,&alpha,&(f1_parameters[1]),&(f1_parameters[2]),&Katz_point_coeff_Gy);
      }
    }
  }// end RDD_KatzPoint & RDD_Site

  if( *rdd_model == RDD_Geiss){
    // Loop over all r_m given
    for (i = 0; i < *n; i++){
      D_RDD_Gy[i]    =  0.0f;
      if (r_m[i] < f1_parameters[1]){
        D_RDD_Gy[i]    =  f1_parameters[5];
      }
      if ((f1_parameters[1] <= r_m[i]) && (r_m[i] <= f1_parameters[2])){
        D_RDD_Gy[i]    =  f1_parameters[5] * (f1_parameters[1] / r_m[i]) * (f1_parameters[1] / r_m[i]);
      }
    }
  }// end RDD_Geiss

  if( *rdd_model == RDD_Cucinotta){
    // Loop over all r_m given
    for (i = 0; i < *n; i++){
      D_RDD_Gy[i]    =  0.0f;
      if (r_m[i] < f1_parameters[1]){
//        D_RDD_Gy[i]    =  f1_parameters[5];
        D_RDD_Gy[i]    =  0.0f;
      }
      if ((f1_parameters[1] <= r_m[i]) && (r_m[i] <= f1_parameters[2])){
//        D_RDD_Gy[i]    =  f1_parameters[5] * (f1_parameters[1] / r_m[i]) * (f1_parameters[1] / r_m[i]);
        D_RDD_Gy[i]    =  0.0f;
      }
    }
  }// end RDD_Cucinotta

  free(f1_parameters);
}

float AT_D_RDD_Gy_solver( const float r , void * params ){
  AT_D_RDD_Gy_parameters* params_struct = (AT_D_RDD_Gy_parameters*)(params);
  *((*params_struct).n) = 1;
  params_struct->r_m = (float*)calloc(1,sizeof(float));
  (params_struct->r_m)[0] = r;
  AT_D_RDD_Gy(  params_struct->n,
          params_struct->r_m,
          params_struct->E_MeV_u,
          params_struct->particle_no,
          params_struct->material_no,
          params_struct->rdd_model,
          params_struct->rdd_parameter,
          params_struct->er_model,
          params_struct->er_parameter,
          params_struct->D_RDD_Gy);
  free(params_struct->r_m);
  return (params_struct->D_RDD_Gy)[0]-(params_struct->D0);
}


void AT_r_RDD_m  ( const  long*  n,
    const float*  D_RDD_Gy,
    /* radiation field parameters */
    const float*  E_MeV_u,
    const long*  particle_no,
    /* detector parameters */
    const long*  material_no,
    /* radial dose distribution model */
    const long*  rdd_model,       /* */
    const float*  rdd_parameter,   /* parameters: LEM: E_MeV_u, particle_no, material_name, a0 */
    /* electron range model */
    const long*  er_model,
    const float*  er_parameter,
    float*  r_RDD_m)
{
  long     i;
  long     n_tmp      = 1;
  // Get f1 parameters
  long     n_f1_parameters = 9;
  float*   f1_parameters   = (float*)calloc(n_f1_parameters, sizeof(float));
  AT_RDD_f1_parameters(      /* radiation field parameters */
    E_MeV_u,
    particle_no,
    /* detector parameters */
    material_no,
    /* radial dose distribution model */
    rdd_model,
    rdd_parameter,
    /* electron range model */
    er_model,
    er_parameter,
    /* calculated parameters */
    f1_parameters);
  // Get material data
  float     density_g_cm3;
  AT_getMaterialData(  &n_tmp,
      material_no,
      &density_g_cm3,
      NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  float  density_kg_m3  =  density_g_cm3 * 1000.0f;

  if( *rdd_model == RDD_Test){
    // Loop over all doses given
    for (i = 0; i < *n; i++){
      r_RDD_m[i]    =  0.0f;
      if (D_RDD_Gy[i] > 0.0f){
        r_RDD_m[i]    =  f1_parameters[2];
      }
    }
  }// end RDD_Test

  if( *rdd_model == RDD_Geiss){
    // Loop over all doses given
    for (i = 0; i < *n; i++){
      // if D is outside the definition, return -1
      r_RDD_m[i]    =  -1.0f;
      if ((f1_parameters[3] <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= f1_parameters[4])){
        r_RDD_m[i]    =  f1_parameters[1] * (float)sqrt(f1_parameters[5] / D_RDD_Gy[i]);}
    }
  }// end RDD_Geiss

  if( (*rdd_model == RDD_KatzPoint) || (*rdd_model == RDD_Site) || (*rdd_model == RDD_ExtTarget) || (*rdd_model == RDD_Edmund)){
    AT_D_RDD_Gy_parameters* params = (AT_D_RDD_Gy_parameters*)calloc(1,sizeof(AT_D_RDD_Gy_parameters));
    params->n             = (long*)calloc(1,sizeof(long));
    params->E_MeV_u       = E_MeV_u;
    params->particle_no   = particle_no;
    params->material_no   = material_no;
    params->rdd_model     = rdd_model;
    params->rdd_parameter = rdd_parameter;
    params->er_model      = er_model;
    params->er_parameter  = er_parameter;
    params->D_RDD_Gy      = (float*)calloc(1,sizeof(float));
    // Loop over all doses given
    float dev = 0.01f;
    float critical_d_Gy    = 0.0f;
    float critical_r_m     = 0.0f;
    float inv2_d_Gy        = 0.0f;
    float inv2_r_m         = 0.0f;
    if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund){
      critical_r_m    = rdd_parameter[0] * (1.0f + 1e-6f);
      inv2_r_m      = fmaxf(rdd_parameter[0], f1_parameters[2] * dev);
    }
    if(*rdd_model == RDD_ExtTarget){
      critical_r_m    = rdd_parameter[1];
      inv2_r_m      = fmaxf(rdd_parameter[0], f1_parameters[2] * dev);
    }
    if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund || *rdd_model == RDD_ExtTarget){
      AT_D_RDD_Gy  (  &n_tmp,                    // Use D(r) to find dose at jump of D_Site
                &critical_r_m,
                /* radiation field parameters */
                E_MeV_u,
                particle_no,
                /* detector parameters */
                material_no,
                /* radial dose distribution model */
                rdd_model,
                rdd_parameter,
                /* electron range model */
                er_model,
                er_parameter,
                &critical_d_Gy);
      AT_D_RDD_Gy  (  &n_tmp,                    // Use D(r) to find dose at jump of D_Site
                &inv2_r_m,
                /* radiation field parameters */
                E_MeV_u,
                particle_no,
                /* detector parameters */
                material_no,
                /* radial dose distribution model */
                rdd_model,
                rdd_parameter,
                /* electron range model */
                er_model,
                er_parameter,
                &inv2_d_Gy);
    }

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
    // Get beta, Z and Zeff
    float  beta   = 0.0f;
    long  Z       = 0;
    float  Z_eff  = 0.0f;
    AT_beta_from_particle_no(  &n_tmp,
                  E_MeV_u,
                  particle_no,
                  &beta);
    AT_Z_from_particle_no(  &n_tmp,
                  particle_no,
                  &Z);
    AT_effective_charge_from_beta(  &n_tmp,
                  &beta,
                  &Z,
                  &Z_eff);

    //////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////
    float alpha        =  1.667f;
    float wmax_MeV     =  0.0f;
    AT_max_E_transfer_MeV(&n_tmp,E_MeV_u,particle_no,&wmax_MeV);
    if(wmax_MeV <= 1e-3){  // if wmax < 1keV
      alpha           = 1.079f;
    }

    //////////////////////// PRELIMINARY: alpha only according to Katz E-R model ////////////////////////////

    float Katz_point_coeff_Gy = AT_RDD_Katz_point_coeff_Gy(&(f1_parameters[5]),&Z_eff,&beta,&alpha,&density_kg_m3,&(f1_parameters[2]));

    for (i = 0; i < *n; i++){
      params->D0 = D_RDD_Gy[i];
      // if D is outside the definition, return -1
      r_RDD_m[i]        =  -1.0f;
      float  solver_accuracy  =  1e-13f;
      if ((f1_parameters[3] <= D_RDD_Gy[i]) && (D_RDD_Gy[i] <= f1_parameters[4])){
        if(D_RDD_Gy[i] >= critical_d_Gy){
          if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund){
            r_RDD_m[i] = rdd_parameter[0];}
          if(*rdd_model == RDD_ExtTarget){
            r_RDD_m[i] = rdd_parameter[1];}
        }
        if(*rdd_model == RDD_Site || *rdd_model == RDD_Edmund || *rdd_model == RDD_ExtTarget){
          if(D_RDD_Gy[i] < critical_d_Gy && D_RDD_Gy[i] >= inv2_d_Gy){
            r_RDD_m[i] = sqrt(Katz_point_coeff_Gy / D_RDD_Gy[i]) * f1_parameters[2];}
          if(D_RDD_Gy[i] < inv2_d_Gy){
            r_RDD_m[i] = zriddr(AT_D_RDD_Gy_solver,
                    (void*)params,
                    f1_parameters[1],
                    f1_parameters[2],
                    solver_accuracy);
          }
        }else{
          r_RDD_m[i] = zriddr(AT_D_RDD_Gy_solver,
                  (void*)params,
                  f1_parameters[1],
                  f1_parameters[2],
                  solver_accuracy);

        }
      }
    }
    free(params->n);
    free(params->D_RDD_Gy);
    free(params);
  }// end RDD_KatzPoint & RDD_Site

  free(f1_parameters);

}


double AT_D_RDD_Gy_int( double  r_m,
    void*   params){
  float D_Gy;
  long  n_tmp              = 1;
  AT_P_RDD_parameters* par = (AT_P_RDD_parameters*)params;
  float fr_m               = (float)r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  &n_tmp,
      &fr_m,
      par->E_MeV_u,
      par->particle_no,
      par->material_no,
      par->rdd_model,
      par->rdd_parameters,
      par->er_model,
      par->er_parameters,
      &D_Gy);

  return (2.0 * M_PI * r_m * (double)D_Gy);
}

double AT_sI_int( double  r_m,
    void*   params){
  double  P = AT_P_RDD(r_m, params);
  return (r_m * P);
}

double AT_P_RDD( double  r_m,
    void*   params)
{
  float  D_Gy;
  long   n_tmp              = 1;
  AT_P_RDD_parameters* par  = (AT_P_RDD_parameters*)params;
  float  fr_m               = (float)r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  &n_tmp,
      &fr_m,
      par->E_MeV_u,
      par->particle_no,
      par->material_no,
      par->rdd_model,
      par->rdd_parameters,
      par->er_model,
      par->er_parameters,
      &D_Gy);

  long gamma_model = GR_GeneralTarget;
  float P;
  AT_gamma_response(  &n_tmp,
      &D_Gy,
      &gamma_model,
      par->gamma_parameters,
      // return
      &P);
  return ((double)P);
}
