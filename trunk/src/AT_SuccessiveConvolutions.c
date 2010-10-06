/**
 * @file
 * @brief ...
 */

/*
 *    AT_SuccessiveConvolutions.c
 *    ===========================
 *
 *    Created on: 28.07.2009
 *    Creator: greilich
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

#include "AT_SuccessiveConvolutions.h"

#define MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT  0.002
#define DEBUG_INTERVALS                      8
#define DEBUG_MEAN                           10.0
#define DEBUG_SIGMA                          1.0


long  AT_n_bins_for_singe_impact_local_dose_distrib(
    const long    n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    N2)
{
  /* get lowest and highest dose */
  double d_max_Gy    =  0.0;
  double d_min_Gy    =  0.0;

  // TODO think if d_min calculations can be done in smarter way. LET is only needed for Geiss RDD

  long  i;
  for (i = 0; i < n; i++){

    double max_electron_range_m = AT_max_electron_range_m( E_MeV_u[i], (int)material_no, (int)er_model);
    double LET_MeV_cm2_g        = AT_LET_MeV_cm2_g_single(E_MeV_u[i], particle_no[i], material_no);
    double norm_constant_Gy     = AT_RDD_precalculated_constant_Gy(max_electron_range_m, LET_MeV_cm2_g, E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model);
    double current_d_min_Gy     = AT_RDD_d_min_Gy( E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model, norm_constant_Gy);
    double current_d_max_Gy     = AT_RDD_d_max_Gy( E_MeV_u[i], particle_no[i], material_no, rdd_model, rdd_parameter, er_model);

    if(i == 0){
      d_min_Gy      =  current_d_min_Gy;
      d_max_Gy      =  current_d_max_Gy;
    }
    else{
      d_min_Gy      =  GSL_MIN(d_min_Gy, current_d_min_Gy);
      d_max_Gy      =  GSL_MAX(d_max_Gy, current_d_max_Gy);
    }
  }

  long n_bins_for_singe_impact_local_dose_ditrib = 0;
  // get number of bins needed to span that dose range
  if( (d_min_Gy > 0) && (d_max_Gy >0) ){
    double tmp        =  log10(d_max_Gy/d_min_Gy) / log10(2.0) * ((double)N2);
    n_bins_for_singe_impact_local_dose_ditrib        =  (long)(floor(tmp) + 1.0);
  } else {
    printf("AT_n_bins_for_singe_impact_local_dose_ditrib: problem in evaluating n_bins_for_singe_impact_local_dose_ditrib: d_min = %g [Gy], d_max = %g [Gy] \n", d_min_Gy, d_max_Gy);
    exit(EXIT_FAILURE);
  }
  return n_bins_for_singe_impact_local_dose_ditrib;
}


void  AT_single_impact_local_dose_distrib(
    const long    n,
    const double  E_MeV_u[],
    const long    particle_no[],
    const double  fluence_cm2_or_dose_Gy[],
    const long    material_no,
    const long    rdd_model,
    const double  rdd_parameter[],
    const long    er_model,
    const long    N2,
    const long    n_bins_f1,
    const double  f1_parameters[],
    double        f1_d_Gy[],
    double        f1_dd_Gy[],
    double        f1[])
{
  double*  fluence_cm2    =  (double*)calloc(n, sizeof(double));

  long i;
  if(fluence_cm2_or_dose_Gy[0] < 0){
    double*  dose_Gy        =  (double*)calloc(n, sizeof(double));
    for (i = 0; i < n; i++){
      dose_Gy[i] = -1.0 * fluence_cm2_or_dose_Gy[i];
    }
    AT_fluence_cm2(  n,
        E_MeV_u,
        particle_no,
        dose_Gy,
        material_no,
        fluence_cm2);
    free( dose_Gy );
  }else{
    for (i = 0; i < n; i++){
      fluence_cm2[i] = fluence_cm2_or_dose_Gy[i];
    }
  }
  double*  norm_fluence                                 =  (double*)calloc(n, sizeof(double));

  // Normalize fluence vector
  AT_normalize(    n,
                fluence_cm2,
                norm_fluence);

  free( fluence_cm2 );

  if(n_bins_f1 > 0){
    double  d_min      =  f1_parameters[0*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3];
    double  d_max      =  f1_parameters[0*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4];

    for (i = 1; i < n; i++){
      d_min          =  GSL_MIN(f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3], d_min);
      d_max          =  GSL_MAX(f1_parameters[i*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4], d_max);
    }

    double  U        =  (log(2.0) / (double)N2);

    double*  d_df_low      =  (double*)calloc(n_bins_f1, sizeof(double));
    double*  d_df_mid      =  (double*)calloc(n_bins_f1, sizeof(double));
    double*  d_df_high     =  (double*)calloc(n_bins_f1, sizeof(double));
    double*  dd_df         =  (double*)calloc(n_bins_f1, sizeof(double));

    for (i = 0; i < n_bins_f1; i++){
      // TODO: check if n.bins sufficient

      d_df_low[i]          =   d_min * exp((double)i * U);
      d_df_high[i]         =   d_min * exp(((double)i + 1.0) * U);
      d_df_mid[i]          =   d_min * exp(((double)i + 0.5) * U);
      dd_df[i]             =   d_df_high[i] - d_df_low[i];              // OBS: using Kellerer's bin-width = mid.point * U is not entirely correct

      f1[i]            =  0.0;
    }

    long n_bins_used = 1;

    // loop over all particles and energies, compute contribution to f1
    long   k;
    for (k = 0; k < n; k++){

      double  d_min_k        =  f1_parameters[k*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 3];
      double  d_max_k        =  f1_parameters[k*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 4];

      // find first and last bin to fit this particle's contribution into the all-over f1-frame
      long  i_low, i_high;
      i_low  = locate(d_df_low, n_bins_f1, d_min_k);
      i_high = locate(d_df_low, n_bins_f1, d_max_k);
      i_low            -=  1;
      i_high           -=  1;

      long  n_bins_df      =  i_high - i_low + 1;  // changed from + 2

      if (n_bins_df > 1){
        double*  d_low       =  (double*)calloc(n_bins_df, sizeof(double));
        double*  d_mid       =  (double*)calloc(n_bins_df, sizeof(double));
        double*  d_high      =  (double*)calloc(n_bins_df, sizeof(double));
        double*  dd          =  (double*)calloc(n_bins_df, sizeof(double));
        double*  r           =  (double*)calloc(n_bins_df, sizeof(double));
        double*  F1_1        =  (double*)calloc(n_bins_df, sizeof(double));
        double*  f1_k        =  (double*)calloc(n_bins_df - 1, sizeof(double));

        // extract the corresponding part from the all-over frame
        for (i = 0; i < n_bins_df; i++){
          d_low[i]           =   d_df_low[i_low + i];
          d_high[i]          =   d_df_high[i_low + i];
          d_mid[i]           =   d_df_mid[i_low + i];
          dd[i]              =   d_high[i] - d_low[i];
        };

        // and adjust the edges
        d_low[0]                =  d_min_k;
        d_low[n_bins_df-1]      =  d_max_k;

        d_mid[0]                =  sqrt(d_low[0] * d_high[0]);
        d_mid[n_bins_df-1-1]    =  sqrt(d_low[n_bins_df - 1] * d_high[n_bins_df - 1]);
        d_mid[n_bins_df-1]      =  0.0;

        d_high[n_bins_df-1-1]   =  d_max_k;
        d_high[n_bins_df-1]     =  0.0;  //TODO ??

        dd[n_bins_df-1]         =  0.0;

        // now compute r, F1, and f1, this could be any RDD if implemented
        int inverse_RDD_status_code = AT_r_RDD_m  (  n_bins_df,
            d_low,
            E_MeV_u[k],
            particle_no[k],
            material_no,
            rdd_model,
            rdd_parameter,
            er_model,
            r);

        if( inverse_RDD_status_code != 0 ){
          printf("Problem in evaluating inverse RDD in AT_SC_get_f1, probably wrong combination of ER and RDD used\n");
          char rdd_model_name[100];
          AT_RDD_name_from_number(rdd_model, rdd_model_name);
          char er_model_name[100];
          getERName( er_model, er_model_name);
          printf("rdd_model: %ld (%s), er_model: %ld (%s)\n", rdd_model, rdd_model_name, er_model, er_model_name);
          exit(EXIT_FAILURE);
        }

        for (i = 0; i < n_bins_df; i++){
          F1_1[i]            = gsl_pow_2(r[i] / f1_parameters[k*AT_SC_F1_PARAMETERS_SINGLE_LENGTH + 2]);        // F1 - 1 instead of F1 to avoid numeric cut-off problems
        }

        F1_1[n_bins_df-1]    =  0.0;

        // now compute f1 as the derivative of F1
        for (i = 0; i < (n_bins_df - 1); i++){
          f1_k[i]          =  -1.0 * (F1_1[i + 1] - F1_1[i]) / dd[i];
        }

        // adjust the density in first and last bin, because upper limit is not d.max.Gy and lower not d.min.Gy
        f1_k[0]              *=  dd[0] / dd_df[i_low];
        f1_k[n_bins_df-2]    *=  dd[n_bins_df - 2] / dd_df[i_high - 1];

        // and paste f1 for this energy /particle into the over all data frame according to rel. fluence
        for (i = 0; i < (n_bins_df - 1); i++){
          f1[i_low + i]      +=  norm_fluence[k] * f1_k[i];
        }

        free(d_low);
        free(d_mid);
        free(d_high);
        free(dd);
        free(r);
        free(F1_1);
        free(f1_k);
      }
      else{ // n_bins_df == 1
        f1[i_low ]        +=  norm_fluence[k] * 1.0 / dd_df[i_low];
      }

      // remember highest bin used
      n_bins_used          =  GSL_MAX(n_bins_used, i_high);
    }

    // copy back for the dose axis
    for (i = 0; i < n_bins_f1; i++){
      f1_d_Gy[i]    =  d_df_mid[i];
      f1_dd_Gy[i]   =  dd_df[i];
    }

    free(d_df_low);
    free(d_df_mid);
    free(d_df_high);
    free(dd_df);

    // normalize f1 (should be ok anyway but there could be small round-off errors)
    double  f1_norm    =  0.0;
    for (i = 0; i < n_bins_f1; i++){
      f1_norm    +=    f1[i] * f1_dd_Gy[i];
    }
    for (i = 0; i < n_bins_f1; i++){
      f1[i]    /=    f1_norm;
    }

  } // if(f1_d_Gy != NULL)

  free( norm_fluence );

}


void  AT_n_bins_for_low_fluence_local_dose_distribution( const double  u,
    const double   fluence_factor,
    const long     N2,
    const long     n_bins_f1,
    const double   f1_d_Gy[],
    const double   f1_dd_Gy[],
    const double   f1[],
    long*          n_bins_f,
    double*        u_start,
    long*          n_convolutions)
{
  // Get expectation value of dose from f1
  double  d_f1_Gy    =  0.0;

  long   i;
  for (i = 0; i < n_bins_f1; i++){
    d_f1_Gy    +=  f1_d_Gy[i] * f1_dd_Gy[i] * f1[i];
  }

  // The dose set by the input data is therefore
  double  d_f_Gy    = u * fluence_factor * d_f1_Gy;

  // How many convolution are necessary starting from a small mean
  // impact number that allows linear approximation (e.g. 0.002)
  *n_convolutions    = 0;

  *u_start      =    d_f_Gy  / d_f1_Gy;
  while(*u_start  > MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT){
    *u_start      =    0.5 * (*u_start);
    (*n_convolutions)++;
  }

  // Every convolution will add a factor of two, so the array for f has to
  // be expanded from f1 by N2 * n_convolutions
  *n_bins_f      =  (*n_convolutions + 1) * N2;
  *n_bins_f      +=  n_bins_f1;
}


void  AT_low_fluence_local_dose_distribution(  const long     n_bins_f1,
    const long     N2,
    const double   f1_d_Gy[],
    const double   f1_dd_Gy[],
    const double   f1[],
    const long     n_bins_f,
    double         f_d_Gy[],
    double         f_dd_Gy[],
    double         f_start[])
{
  // temporary arrays
  double*  d_low         =  (double*)calloc(n_bins_f, sizeof(double));
  double*  d_high        =  (double*)calloc(n_bins_f, sizeof(double));

  double  U              =  log(2.0) / (double)N2;

  long   i;
  for (i = 0; i < n_bins_f; i++){
    d_low[i]            =   f1_d_Gy[0] * exp(((double)i - 0.5)* U);
    d_high[i]           =   f1_d_Gy[0] * exp(((double)i + 0.5)* U);
    if (i < n_bins_f1){
      f_d_Gy[i]         =  f1_d_Gy[i];
      f_dd_Gy[i]        =  f1_dd_Gy[i];
      f_start[i]        =  f1[i];
    }else{
      f_d_Gy[i]         =  sqrt(d_low[i] * d_high[i]);
      f_dd_Gy[i]        =  d_high[i] - d_low[i];
      f_start[i]        =  0.0;
    }
  }

  free(d_low);
  free(d_high);
}


AT_aKList  AT_SC_NORMAL(AT_aKList theKList){

  double  Y        =  theKList.CM1 * 2;
  double  Z        =  theKList.CM2 * 2;
  double  CM0      =  theKList.H0;
  theKList.CM1     =  0;

  long    N        =  theKList.MIH - theKList.MIE;
  long     L;
  for (L = 1; L <= theKList.LEH; L++){
    long    LE     =  L + N;
    double  S      =  theKList.H[L-1] * theKList.DE[LE-1];

    CM0            =  CM0 + S;
    theKList.CM1   =  theKList.CM1 + S * theKList.E[LE-1];
  }

  double  TT       =  (1.0 - theKList.H0) / (CM0 - theKList.H0);
  theKList.CM1     =  theKList.CM1 * TT;
  double  R        =  theKList.CM1 * theKList.CM1;
  theKList.CM2     =  R * theKList.H0;
  theKList.CM3     =  -1.0 * theKList.CM1 * R * theKList.H0;
  theKList.CM4     =  R * R * theKList.H0;

  for (L = 1; L <= theKList.LEH; L++){
    long   LE        =    L + N;
    double  EC       =    theKList.E[LE-1] - theKList.CM1;
    double  E2       =    EC * EC;
    theKList.H[L-1]  =    theKList.H[L-1] * TT;
    double  S        =    theKList.H[L-1] * theKList.DE[LE-1] * E2;
    theKList.CM2     =    theKList.CM2 + S;
    theKList.CM3     =    theKList.CM3 + S * EC;
    theKList.CM4     =    theKList.CM4 + S * E2;
  }

  theKList.X      =  theKList.X * CM0;
  Y          =  theKList.CM1 / Y;
  Z          =  theKList.CM2 / Z;

  return(theKList);
}


AT_aKList  AT_SC_OUTPUT(AT_aKList theKList){

  //double*  SD            =  (double*)calloc(theKList.array_size, sizeof(double));

  double  B             =  theKList.CM2 / (theKList.CM1 * theKList.CM1);
  double  C             =  theKList.CM3 / sqrt(theKList.CM2 * theKList.CM2 * theKList.CM2);
  double  D             =  theKList.CM4 / (theKList.CM2 * theKList.CM2);
  double  S1            =  theKList.CN * theKList.D1;
  double  S2            =  theKList.CN * theKList.D2 / (S1 * S1);
  double  S3            =  theKList.D3 / sqrt(theKList.CN * theKList.D2 * theKList.D2 * theKList.D2);
  double  S4            =  theKList.D4 / (theKList.D2 * theKList.D2 * theKList.CN) + 3;

  if(theKList.N1 <= 0){
    S2                =  B;
    S3                =  C;
    S4                =  D;
  }

  if(theKList.write_output){
    long   L;
    for (L = 1; L <= theKList.array_size; L++){
      fprintf(  theKList.output_file,
          "%ld;%e;%e;%e;%e;%e;%ld;%ld;%ld;%ld;%ld;%ld\n",
          L,
          theKList.E[L-1],
          theKList.DE[L-1],
          theKList.H[L-1],
          theKList.H0,
          theKList.F[L-1],
          theKList.N1,
          theKList.MIF,
          theKList.LEF,
          theKList.MIH,
          theKList.LEH,
          theKList.MIE);
    }
  }

  return(theKList);
}


AT_aKList  AT_SC_INTERP(AT_aKList theKList){

  theKList.A[0]             =  theKList.F[1] - theKList.F[0];
  theKList.BI[0]            =  0.0;
  theKList.F[theKList.LEF]  =  0.0;

  long   K;
  for (K = 1; K <= theKList.N2; K++){
    long L           =  theKList.LEF + K;
    theKList.A[L-1]  =  0.0;
    theKList.BI[L-1] =  0.0;
  }

  long   L;
  for (L = 2; L <= theKList.LEF; L++){
    theKList.A[L -1]    =  0.5 * (theKList.F[L] - theKList.F[L - 1 -1]);
    theKList.BI[L -1]   =  theKList.A[L-1] + theKList.F[L - 1 -1] - theKList.F[L -1];
  }

  return(theKList);
}


AT_aKList  AT_SC_RESET(AT_aKList theKList){

  if (theKList.N2 <= 256){
    if(theKList.LEF <= 64){
      double S              =  log(2.0);
      double TT             =  (double)theKList.N2;
      //      theKList.N2            =  theKList.N2 * 2;
      theKList.N2          +=  (long)(0.1 + exp((double)((long)(log(TT) / S - 0.99)) * S));
      TT                   =  TT / (double)theKList.N2;
      theKList.U           =  S / (double)theKList.N2;

      theKList            =  AT_SC_INTERP(theKList);
      theKList.F[theKList.LEF]  =  0;
      long N              =  theKList.MIF;
      theKList.MIF          =  (long)((double)theKList.MIF / TT) + 1;    /////////////////////
      theKList.LEF          =  (long)((double)theKList.LEF / TT) - 1;    // added (SG) : -1 //
      /////////////////////
      long   K;
      for (K = 1; K <= theKList.LEF; K++){
        long   L              =  theKList.LEF - K + 1;
        double  FLF           =  (double)(L + theKList.MIF) * TT - (double)N;
        long   LFF            =  (long)(FLF + 0.5);
        double S              =  FLF - (double)LFF;

        ////////////////////////////////////////////////////////////////////////////////
        // Replaced Kellerer's original quadratic interpolation by a slower           //
        // but more correct approach. The original produced wrong interpolation &     //
        // negative values of F on the new N2's grid at very steep, irregular parts   //
        // of F  and eventually systematic deviation of moments.                       //
        ////////////////////////////////////////////////////////////////////////////////

        // original:  theKList.F[L -1]        =  theKList.F[LFF -1] + S * (theKList.A[LFF - 1] + S * theKList.BI[LFF - 1]);

        theKList.F[L -1]        =  theKList.F[LFF -1];
        if((S < 0 ) && (LFF >= 2)){
          theKList.F[L -1]      =  pow(theKList.F[LFF - 1 -1], -1.0 * S) * pow(theKList.F[LFF -1], 1.0 + S);
        }
        if((S > 0 ) && (LFF <= theKList.LEF - 1)){
          theKList.F[L -1]      =  pow(theKList.F[LFF -1], 1.0 - S) * pow(theKList.F[LFF], S);
        }
      }

      long   L;
      for (L = theKList.N2; L <= theKList.array_size; L++){;
      double S           =  (double)(L - theKList.N2) * theKList.U;
      double tmp         =  -1.0 * log(1.0 - 0.5 * exp(-S)) / theKList.U;
      theKList.DI[L -1]  =  tmp - (double)theKList.N2;    // type casts necessary to prevent round of errors (that will eventually yield negative H-values in AT_SC_FOLD
      }

      theKList.MIE          =  theKList.MIF;

      long   J;
      for (J = 1; J <= theKList.array_size; J++){
        double S            =  (double)(J + theKList.MIE);
        theKList.E[J -1]    =  exp(S * theKList.U) * theKList.E0;
        ///////////////////////////////////////////////////////////////////////////
        // addition SG: not to use Kellerer's formula for new DE's, as it is     //
        // not exact (but deviation are small)                                   //
        ///////////////////////////////////////////////////////////////////////////
        double* high_E      =  (double*)calloc(theKList.array_size, sizeof(double));
        S                   =  (double)(J + theKList.MIE + 1);
        high_E[J - 1]       =  exp(S * theKList.U) * theKList.E0;
        theKList.DE[J -1]   =  high_E[J -1] - theKList.E[J -1];
        free(high_E);
      }
    }else{
      return(theKList);
    }
  }else{
    theKList.MIE          =  theKList.MIF;

    long   J;
    for (J = 1; J <= theKList.array_size; J++){
      double S             =  (double)(J + theKList.MIE);
      theKList.E[J -1]     =  exp(S * theKList.U) * theKList.E0;
      ///////////////////////////////////////////////////////////////////////////
      // addition SG: not to use Kellerer's formula for new DE's, as it is     //
      // not exact (but deviation are small)                                   //
      ///////////////////////////////////////////////////////////////////////////
      double* high_E       =  (double*)calloc(theKList.array_size, sizeof(double));
      S                    =  (double)(J + theKList.MIE + 1);
      high_E[J - 1]        =  exp(S * theKList.U) * theKList.E0;
      theKList.DE[J -1]    =  high_E[J -1] - theKList.E[J -1];
      free(high_E);
    }
  }

  return(theKList);
}


AT_aKList  AT_SC_ZERO(AT_aKList theKList){

  theKList.X          =  0;
  long N              =  theKList.MIH - theKList.MIE;

  long   L;
  for (L = 1; L <= theKList.LEH; L++){
    long K            =  L + N;
    theKList.X        =  theKList.X + theKList.H[L -1] * theKList.DE[K -1];
  }

  double S            =  (1.0 - theKList.F0) * (1.0 - theKList.F0) / theKList.X;
  theKList.X          =  2.0 / S;

  for (L = 1; L <= theKList.LEH; L++){;
  theKList.H[L -1]    =  theKList.H[L -1] * S;
  }

  N                   =  theKList.MIH - theKList.MIF;
  theKList.MIH        =  theKList.MIF;
  theKList.LEH        =  theKList.LEH + N;

  long   LL;
  for (LL = 1; LL <= theKList.LEH; LL++){
    long L            =  theKList.LEH + 1 - LL;
    long K            =  L + N;
    theKList.H[K -1]  =  theKList.H[L -1];
  }

  for (L = 1; L <= N; L++){
    theKList.H[L -1]  =  0.0;
  }

  S                   =  theKList.F0 * 2.0;

  for (L = 1; L <= theKList.LEF; L++){
    theKList.H[L -1]  =  theKList.H[L -1] + theKList.F[L -1] * S;
  }

  return(theKList);
}


AT_aKList  AT_SC_SHRINK(AT_aKList theKList){

  double  EX          =  theKList.shrink_tails_under;
  double  S           =  0.0;
  long  N             =  theKList.MIH - theKList.MIE;

  long   L;
  for (L = 1; L <= theKList.LEH; L++){
    long K            =  L + N;
    S                 =  S + theKList.H[L -1] * theKList.DE[K -1];
    if(S > 1000.0 * EX){
      theKList.MIH    =  theKList.MIH + L - 1;
      break;}
  }

  long    M           =  L - 1;
  S                   =  0;

  long   K;
  for (K = 1; K <= theKList.LEH; K++){
    L                 =  theKList.LEH + 1 - K;
    long KK           =  L + N;
    S                 =  S + theKList.H[L - 1] * theKList.DE[KK - 1];
    if(S > EX){
      break;
    }
  }

  theKList.LEH        =  L - M;
  for (L = 1; L <= theKList.LEH; L++){
    K                 =  L + M;
    theKList.H[L -1]  =  theKList.H[K -1];
  }

  K                   =  theKList.LEH + 1;
  long  KK            =  theKList.LEH + M;
  for (L = K; L <= KK; L++){
    theKList.H[L -1]  =  0;
  }

  return(theKList);

}


AT_aKList AT_SC_FOLD(AT_aKList theKList){
  double*  FDE        =  (double*)calloc(theKList.array_size, sizeof(double));

  if((theKList.CN >= 10.0) && (theKList.adjust_N2 == true)){
    theKList          =  AT_SC_RESET(theKList);
  }

  theKList.H0         =  theKList.F0 * theKList.F0;
  theKList.MIH        =  theKList.MIF + theKList.N2;
  theKList.LEH        =  theKList.LEF;
  long  K             =  theKList.LEF + 1;
  long KK             =  K + theKList.N2;

  long   L;
  for (L = K; L <= KK; L++){
    theKList.F[L -1]  =  0;
  }

  theKList            =  AT_SC_INTERP(theKList);
  long N              =  theKList.MIF - theKList.MIE;

  for (L = 1; L <= theKList.LEH; L++){
    K                 =  L + N;
    FDE[L -1]         =  theKList.F[L -1] * theKList.DE[K -1];
  }

  long   LH;
  for (LH = 1; LH <= theKList.LEH; LH++){
    double   HLH      =  0;
    long   LL         =  LH + theKList.N2;
    long   LF;
    for (LF = 1; LF <= LH; LF++){
      K               =  LL - LF;
      double FLF      =  (double)LH - theKList.DI[K -1];
      long LFF        =  (long)(FLF + 0.5);
      double S        =  FLF - (double)LFF;
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Modification SG: if Kellerer's quadratic interpolation fails, use simple estimate
      double tmp      =  theKList.F[LFF -1] + S * (theKList.A[LFF -1] + S * theKList.BI[LFF -1]);
      if (tmp <0){
        tmp = 0.0;        // Very crude - better to replace by interpolation as done in RESET
      }                   // which is time-consuming, however.
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      HLH             =  HLH + FDE[LF -1] * tmp;
    }
    theKList.H[LH -1] =  HLH - FDE[LH -1] * theKList.F[LH -1] * 0.5;
  }

  free(FDE);

  if (theKList.F0 < 1e-10){
    theKList.X        =  2.0;
  }else{
    theKList          =  AT_SC_ZERO(theKList);
  }
  return(theKList);
}


void   AT_SuccessiveConvolutions( const double  u,
    const long    n_bins_f,
    long*         N2,
    long*         n_bins_f_used,
    double        f_d_Gy[],
    double        f_dd_Gy[],
    double        f[],
    double*       f0,
    double        fdd[],
    double        dfdd[],
    double*       d,
    const bool    write_output,
    const bool    shrink_tails,
    const double  shrink_tails_under,
    const bool    adjust_N2)
{
  //////////////////////////////////////////
  // Init KList structure
  // (Constructor)
  //////////////////////////////////////////
  AT_aKList KList;

  KList.array_size  = n_bins_f;
  KList.N2          = *N2;
  KList.U           = log(2.0) / (double)KList.N2;

  //////////////////////////////////

  KList.write_output    =  write_output;
  KList.output_file     =  NULL;
  if(KList.write_output){
    KList.output_file    =  fopen("SuccessiveConvolutions.log","w");
    if (KList.output_file == NULL) return;                      // File error
    fprintf(KList.output_file, "i;E;DE;H;H0;F;convolution.no;MIF;LEF;MIH;LEH;MIE\n");
   }

  //////////////////////////////////

  KList.shrink_tails        =  shrink_tails;
  KList.shrink_tails_under  =  shrink_tails_under;
  KList.adjust_N2           =  adjust_N2;

  //////////////////////////////////

  KList.F        = (double*)calloc(KList.array_size, sizeof(double));
  KList.H        = (double*)calloc(KList.array_size, sizeof(double));
  KList.E        = (double*)calloc(KList.array_size, sizeof(double));
  KList.DE       = (double*)calloc(KList.array_size, sizeof(double));
  KList.DI       = (double*)calloc(KList.array_size, sizeof(double));
  KList.A        = (double*)calloc(KList.array_size, sizeof(double));
  KList.BI       = (double*)calloc(KList.array_size, sizeof(double));

  // Some other initializations
  KList.MIH      = 0;
  KList.MIE      = 0;
  KList.N1       = 0;
  KList.H0       = 0;
  KList.X        = 1;
  KList.CN       = 1;
  KList.CM1      = 1;
  KList.CM2      = 1;

  // Added by Leszek
  KList.MIF      = 0;
  KList.LEF      = 0;

  // Copy input data
  KList.E0      = f_d_Gy[0] * exp(-1.0 * KList.U);

  long   L;
  for (L = 1; L <= KList.array_size; L++){
    KList.E[L -1]      = f_d_Gy[L -1];
    KList.DE[L -1]     = f_dd_Gy[L -1];
    KList.H[L -1]      = f[L -1];
  }

  KList.LEH            = *n_bins_f_used;

  ///////////////////////////////////////
  // Fill array for auxilary function that enables easy index operations
  for  (L = KList.N2; L <= KList.array_size; L++){
    double S          =  (double)(L - KList.N2) * KList.U;
    double tmp        =  -1.0 * log(1.0 - 0.5 * exp(-S)) / KList.U;
    KList.DI[L -1]    =  tmp - (double)KList.N2;
  }    // type casts necessary to prevent round of errors (that will eventually yield negative H-values in AT_SC_FOLD

  ///////////////////////////////////////
  // Normalize distribution
  ///////////////////////////////////////
  KList  = AT_SC_NORMAL(KList);

  if(KList.write_output){
    long i;
    for (i = 0; i < KList.array_size; i++){
        fprintf(  KList.output_file,
            "%ld;%e;%e;%e;%e;%e;%ld;%ld;%ld;%ld;%ld;%ld\n",
            i,
            KList.E[L-1],
            KList.DE[L-1],
            KList.H[L-1],
            KList.H0,
            KList.F[L-1],
            KList.N1,
            KList.MIF,
            KList.LEF,
            KList.MIH,
            KList.LEH,
            KList.MIE);
    }
  }

  ///////////////////////////////////////
  // Get moments of single impact f1
  ///////////////////////////////////////

  KList.D1    =    KList.CM1;
  double  S   =    KList.D1 * KList.D1;
  KList.D2    =    KList.CM2 + S;
  KList.D3    =    KList.CM3 + 3.0 * KList.CM2 * KList.D1 + S * KList.D1;
  KList.D4    =    KList.CM4 + 4.0 * KList.CM3 * KList.D1 + 6.0 * S * KList.CM2 + S * S;

  double  S2        =    KList.D2 / KList.D1;
  double  S3        =    KList.D3 / KList.D1;
  double  S4        =    KList.D4 / KList.D1;
  S                 =    S3 / sqrt(gsl_pow_3(S2));
  double  TT        =    S4  / gsl_pow_2(S2);

  ///////////////////////////////////////
  // AT_SC_SHRINK distribution
  ///////////////////////////////////////
  if(KList.shrink_tails){
    KList  = AT_SC_SHRINK(KList);
  }

  ///////////////////////////////////////
  // AT_SC_OUTPUT
  ///////////////////////////////////////
  KList  = AT_SC_OUTPUT(KList);

  ///////////////////////////////////////
  // Get approximation for small hit
  // numbers
  ///////////////////////////////////////

  KList.FINAL        = u * KList.D1;  // Final mean impact number

  long n_convolutions    = 0;
  KList.CN               =    KList.FINAL  / KList.D1; // i.e. CN = u , TODO remove FINAL
  while(KList.CN > MEAN_HIT_NUMBER_LINEAR_APPROX_LIMIT){
    KList.CN             =    0.5 * KList.CN;
    n_convolutions++;
  }

  KList.H0         =    1.0 - KList.CN;

  for (L = 1; L <= KList.LEH; L++){
    KList.H[L -1]  =  KList.H[L -1] * KList.CN;
  }

  ///////////////////////////////////////
  // Convolution loop
  ///////////////////////////////////////
  long   j;
  for(j = 0; j < n_convolutions; j++){
    KList.N1        =  KList.N1 + 1;
    KList.CN        =  KList.CN * 2.0;

    for (L = 1; L <= KList.LEH; L++){
      KList.F[L -1]      =  KList.H[L -1];
    }

    KList.F0         =  KList.H0;
    KList.LEF        =  KList.LEH;
    KList.MIF        =  KList.MIH;
    KList            =  AT_SC_FOLD(KList);
    if(KList.shrink_tails){
      KList          =  AT_SC_SHRINK(KList);
    }
    KList            =  AT_SC_NORMAL(KList);
    KList            =  AT_SC_OUTPUT(KList);
  }


  //////////////////////////////////////////
  // Copy results back to input structure
  // and adjust according to MIH, MIE
  //////////////////////////////////////////

  *d    = 0.0;

  for (L = 1; L <= KList.array_size; L++){
    f_d_Gy[L -1]      =  0.0;
    f_dd_Gy[L -1]     =  0.0;
    f[L -1]           =  0.0;
    fdd[L -1]         =  0.0;
    dfdd[L -1]        =  0.0;
  }

  long  N        = KList.MIH - KList.MIE;
  for (L = 1; L <= KList.LEH; L++){
    long LE          =  L + N;
    f_d_Gy[L -1]     =  KList.E[LE -1];
    f_dd_Gy[L -1]    =  KList.DE[LE -1];
    f[L -1]          =  KList.H[L-1];
    fdd[L -1]        =  f[L -1] * f_dd_Gy[L -1];
    dfdd[L -1]       =  fdd[L -1] * f_d_Gy[L -1];
    *d              +=  dfdd[L -1];
  }

  *n_bins_f_used = KList.LEH;

  *f0            = KList.H0;
  *N2            = KList.N2;      // could have been changed by RESET --> report back

  //////////////////////////////////////////
  // Free allocated KList structures
  //////////////////////////////////////////
  free(KList.F);
  free(KList.H);
  free(KList.E);
  free(KList.DE);
  free(KList.DI);
  free(KList.A);
  free(KList.BI);

  if(KList.write_output){
    // Close file
    fclose(KList.output_file);
  }
}
