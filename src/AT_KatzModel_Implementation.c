/**
 * @brief TODO
 */

/*
 *    AT_KatzModel.c
 *    ===========================
 *
 *    Created on: 01.03.2010
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

#include "AT_KatzModel_Implementation.h"


double AT_KatzModel_KatzExtTarget_inactivation_probability(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const long    er_model,
    const double  alpha,
    const double  Katz_plateau_Gy,
    const double  Katz_point_coeff_Gy,
    const double  D0_Gy,
    const double  c,
    const double  m){

  if( (r_m >= 0.0) && (r_m <= a0_m + max_electron_range_m ) ){
    double        D_Gy = AT_RDD_ExtendedTarget_KatzPoint_Gy( r_m, a0_m, er_model, KatzPoint_r_min_m, max_electron_range_m, alpha, Katz_plateau_Gy, Katz_point_coeff_Gy);

    const long    gamma_model                =  GR_GeneralTarget;
    const double  gamma_parameters[5]        =  {1.0, D0_Gy, c, m, 0.0};
    double        inactivation_probability;
    const long    n_tmp                      =  1;
    AT_gamma_response(  n_tmp,
        &D_Gy,
        gamma_model,
        gamma_parameters,
        false,
        &inactivation_probability);

    return inactivation_probability;
  } else {
    return 0.0;
  }
}


double AT_KatzModel_CucinottaExtTarget_inactivation_probability(
    const double  r_m,
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy,
    const double  KatzPoint_point_coeff_Gy,
    const double  D0_Gy,
    const double  c,
    const double  m){

  if( (r_m >= 0.0) && (r_m <= a0_m + max_electron_range_m ) ){
    double        D_Gy = AT_RDD_ExtendedTarget_CucinottaPoint_Gy( r_m, a0_m, KatzPoint_r_min_m, max_electron_range_m, beta, KatzPoint_point_coeff_Gy, C_norm, Cucinotta_plateau_Gy);

    const long    gamma_model                =  GR_GeneralTarget;
    const double  gamma_parameters[5]        =  {1.0, D0_Gy, c, m, 0.0};
    double        inactivation_probability;
    const long    n_tmp                      =  1;
    AT_gamma_response(  n_tmp,
        &D_Gy,
        gamma_model,
        gamma_parameters,
        false,
        &inactivation_probability);

    return inactivation_probability;
  } else {
    return 0.0;
  }
}


int AT_KatzModel_inactivation_probability(
    const long    n,
    const double  r_m[],
    const double  E_MeV_u,
    const long    particle_no,
    const long    rdd_model,
    const double  rdd_parameters[],
    const long    er_model,
    const double  gamma_parameters[],
    const long    stop_power_source,    
    double        inactivation_probability[]){

  const double max_electron_range_m  =  AT_max_electron_range_m( E_MeV_u, Water_Liquid, (int)er_model);

  const double a0_m                  =  AT_RDD_a0_m( max_electron_range_m, rdd_model, rdd_parameters );

  const double KatzPoint_r_min_m     =  AT_RDD_r_min_m( max_electron_range_m, rdd_model, rdd_parameters );

  const double Katz_point_coeff_Gy   =  AT_RDD_Katz_coeff_Gy_general( E_MeV_u, particle_no, Water_Liquid, er_model);

  const double r_max_m               =  GSL_MIN(a0_m, max_electron_range_m);

  const double D0_Gy  =  gamma_parameters[1];
  const double c      =  gamma_parameters[2];
  const double m      =  gamma_parameters[3];

  if( rdd_model == RDD_KatzExtTarget ){
    double Katz_plateau_Gy  =  0.0;
    double alpha                       =  0.0;
    if( (er_model == ER_Waligorski) || (er_model == ER_Edmund) ){ // "new" Katz RDD
      alpha            =  AT_ER_PowerLaw_alpha(E_MeV_u);
      Katz_plateau_Gy  =  AT_RDD_Katz_PowerLawER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, alpha, Katz_point_coeff_Gy );
    } else if (er_model == ER_ButtsKatz){ // "old" Katz RDD
      Katz_plateau_Gy  =  AT_RDD_Katz_LinearER_Daverage_Gy( KatzPoint_r_min_m, r_max_m, max_electron_range_m, Katz_point_coeff_Gy );
    }

    long i;
    for( i = 0 ; i < n ; i++){
      inactivation_probability[i] = AT_KatzModel_KatzExtTarget_inactivation_probability(
              r_m[i],
              a0_m,
              KatzPoint_r_min_m,
              max_electron_range_m,
              er_model,
              alpha,
              Katz_plateau_Gy,
              Katz_point_coeff_Gy,
              D0_Gy,
              c,
              m);
    }
  }

    if (rdd_model == RDD_CucinottaExtTarget) {
        double LET_MeV_cm2_g;
        AT_Mass_Stopping_Power_with_no(stop_power_source,
                                       1,
                                       &E_MeV_u,
                                       &particle_no,
                                       Water_Liquid,
                                       &LET_MeV_cm2_g);
        const double density_g_cm3 = AT_density_g_cm3_from_material_no(Water_Liquid);
        const double density_kg_m3 = density_g_cm3 * 1000.0;
        const double LET_J_m = LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
        const double beta = AT_beta_from_E_single(E_MeV_u);
        const double C_norm = AT_RDD_Cucinotta_Cnorm(KatzPoint_r_min_m, max_electron_range_m, beta, density_kg_m3,
                                                     LET_J_m, Katz_point_coeff_Gy);
        double Cucinotta_plateau_Gy = AT_RDD_Cucinotta_Ddelta_average_Gy(KatzPoint_r_min_m, r_max_m,
                                                                         max_electron_range_m, beta,
                                                                         Katz_point_coeff_Gy);
        Cucinotta_plateau_Gy += C_norm *
                                AT_RDD_Cucinotta_Dexc_average_Gy(KatzPoint_r_min_m, r_max_m, max_electron_range_m, beta,
                                                                 Katz_point_coeff_Gy);

        long i;
        for (i = 0; i < n; i++) {
            inactivation_probability[i] = AT_KatzModel_CucinottaExtTarget_inactivation_probability(
                    r_m[i],
                    a0_m,
                    KatzPoint_r_min_m,
                    max_electron_range_m,
                    beta,
                    C_norm,
                    Cucinotta_plateau_Gy,
                    Katz_point_coeff_Gy,
                    D0_Gy,
                    c,
                    m);
        }
    }


  return EXIT_SUCCESS;
}


double AT_KatzModel_KatzExtTarget_inactivation_cross_section_integrand_m(
    double t_m,
    void* params){

  assert( params != NULL );
  AT_KatzModel_KatzExtTarget_inactivation_probability_parameters * inact_prob_parameters = (AT_KatzModel_KatzExtTarget_inactivation_probability_parameters*)params;

  const double  inactivation_probability   =  AT_KatzModel_KatzExtTarget_inactivation_probability( t_m,
      inact_prob_parameters->a0_m,
      inact_prob_parameters->KatzPoint_r_min_m,
      inact_prob_parameters->max_electron_range_m,
      inact_prob_parameters->er_model,
      inact_prob_parameters->alpha,
      inact_prob_parameters->Katz_plateau_Gy,
      inact_prob_parameters->Katz_point_coeff_Gy,
      inact_prob_parameters->D0_Gy,
      inact_prob_parameters->c,
      inact_prob_parameters->m);

  return inactivation_probability * t_m;
}


double AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2(
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const long    er_model,
    const double  alpha,
    const double  Katz_plateau_Gy,
    const double  Katz_point_coeff_Gy,
    const double  D0_Gy,
    const double  c,
    const double  m){

  double low_lim_m = 0.01*a0_m;
  gsl_set_error_handler_off();

  double integral_m2;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &AT_KatzModel_KatzExtTarget_inactivation_cross_section_integrand_m;

  AT_KatzModel_KatzExtTarget_inactivation_probability_parameters inact_prob_parameters;

  inact_prob_parameters.a0_m                       =  a0_m;
  inact_prob_parameters.KatzPoint_r_min_m          =  KatzPoint_r_min_m;
  inact_prob_parameters.max_electron_range_m       =  max_electron_range_m;
  inact_prob_parameters.er_model                   =  er_model;
  inact_prob_parameters.alpha                      =  alpha;
  inact_prob_parameters.Katz_plateau_Gy            =  Katz_plateau_Gy;
  inact_prob_parameters.Katz_point_coeff_Gy        =  Katz_point_coeff_Gy;
  inact_prob_parameters.D0_Gy  =  D0_Gy;
  inact_prob_parameters.c      =  c;
  inact_prob_parameters.m      =  m;

  F.params = (void*)(&inact_prob_parameters);
  int status = gsl_integration_qag (&F, low_lim_m, max_electron_range_m + a0_m, 0, 1e-4, 1000, GSL_INTEG_GAUSS21, w1, &integral_m2, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
#ifndef NDEBUG
	  fprintf(stderr,"Error in AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2: er_model = %ld, integration from %g to %g [m] + %g [m]\n", er_model, low_lim_m, max_electron_range_m, a0_m);
#endif
	  integral_m2 = 0.0;
  }
  gsl_integration_workspace_free (w1);

  const long    gamma_model                =  GR_GeneralTarget;
  const double  gamma_parameters[5]        =  {1.0, D0_Gy, c, m, 0.0};
  double        inactivation_probability_plateau;
  const long    n_tmp                      =  1;
  AT_gamma_response(  n_tmp,
      &Katz_plateau_Gy,
      gamma_model,
      gamma_parameters,
      false,
      &inactivation_probability_plateau);

  return 2.0 * M_PI * integral_m2 + M_PI * gsl_pow_2(0.01*a0_m) * inactivation_probability_plateau;
}


double AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_integrand_m(
    double t_m,
    void* params){

  assert( params != NULL );
  AT_KatzModel_CucinottaExtTarget_inactivation_probability_parameters * inact_prob_parameters = (AT_KatzModel_CucinottaExtTarget_inactivation_probability_parameters*)params;

  const double  inactivation_probability   =  AT_KatzModel_CucinottaExtTarget_inactivation_probability( t_m,
      inact_prob_parameters->a0_m,
      inact_prob_parameters->KatzPoint_r_min_m,
      inact_prob_parameters->max_electron_range_m,
      inact_prob_parameters->beta,
      inact_prob_parameters->C_norm,
      inact_prob_parameters->Cucinotta_plateau_Gy,
      inact_prob_parameters->KatzPoint_coeff_Gy,
      inact_prob_parameters->D0_Gy,
      inact_prob_parameters->c,
      inact_prob_parameters->m);

  return inactivation_probability * t_m;
}


double AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2(
    const double  a0_m,
    const double  KatzPoint_r_min_m,
    const double  max_electron_range_m,
    const double  beta,
    const double  C_norm,
    const double  Cucinotta_plateau_Gy,
    const double  KatzPoint_point_coeff_Gy,
    const double  D0_Gy,
    const double  c,
    const double  m){

  double low_lim_m = 0.01*a0_m;
  gsl_set_error_handler_off();

  double integral_m2;
  double error;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_integrand_m;

  AT_KatzModel_CucinottaExtTarget_inactivation_probability_parameters inact_prob_parameters;

  inact_prob_parameters.a0_m                    =  a0_m;
  inact_prob_parameters.KatzPoint_r_min_m       =  KatzPoint_r_min_m;
  inact_prob_parameters.max_electron_range_m    =  max_electron_range_m;
  inact_prob_parameters.beta                    =  beta;
  inact_prob_parameters.C_norm                  =  C_norm;
  inact_prob_parameters.Cucinotta_plateau_Gy    =  Cucinotta_plateau_Gy;
  inact_prob_parameters.KatzPoint_coeff_Gy      =  KatzPoint_point_coeff_Gy;
  inact_prob_parameters.D0_Gy                   =  D0_Gy;
  inact_prob_parameters.c                       =  c;
  inact_prob_parameters.m                       =  m;

  F.params = (void*)(&inact_prob_parameters);
  int status = gsl_integration_qag (&F, low_lim_m, max_electron_range_m + a0_m, 0, 1e-4, 1000, GSL_INTEG_GAUSS21, w1, &integral_m2, &error);
  if (status == GSL_EROUND || status == GSL_ESING){
#ifndef NDEBUG
    fprintf(stderr,"Error in AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2: integration from %g to %g [m] + %g [m]\n", low_lim_m, max_electron_range_m, a0_m);
#endif
    integral_m2 = 0.0;
  }
  gsl_integration_workspace_free (w1);

  const long    gamma_model                =  GR_GeneralTarget;
  const double  gamma_parameters[5]        =  {1.0, D0_Gy, c, m, 0.0};
  double        inactivation_probability_plateau;
  const long    n_tmp                      =  1;
  AT_gamma_response(  n_tmp,
      &Cucinotta_plateau_Gy,
      gamma_model,
      gamma_parameters,
      false,
      &inactivation_probability_plateau);

 return 2.0 * M_PI * integral_m2 + M_PI * gsl_pow_2(0.01*a0_m) * inactivation_probability_plateau;

}


int AT_KatzModel_inactivation_cross_section_m2(
    const long   n,
    const double E_MeV_u[],
    const long   particle_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double gamma_parameters[],
    const long   stop_power_source,
    double inactivation_cross_section_m2[]){

    const double D0_Gy = gamma_parameters[1];
    const double c = gamma_parameters[2];
    const double m = gamma_parameters[3];

    if (rdd_model == RDD_KatzExtTarget) {
        long i;
        const double a0_m = rdd_parameters[1];
        for (i = 0; i < n; i++) {
            const double max_electron_range_m = AT_max_electron_range_m(E_MeV_u[i], Water_Liquid, (int) er_model);
            const double KatzPoint_r_min_m = AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameters);
            const double Katz_point_coeff_Gy = AT_RDD_Katz_coeff_Gy_general(E_MeV_u[i], particle_no, Water_Liquid,
                                                                            er_model);
            const double r_max_m = GSL_MIN(a0_m, max_electron_range_m);

            double Katz_plateau_Gy = 0.0;
            double alpha = 0.0;
            if ((er_model == ER_Waligorski) || (er_model == ER_Edmund)) { // "new" Katz RDD
                alpha = AT_ER_PowerLaw_alpha(E_MeV_u[i]);
                Katz_plateau_Gy = AT_RDD_Katz_PowerLawER_Daverage_Gy(KatzPoint_r_min_m, r_max_m, max_electron_range_m,
                                                                     alpha, Katz_point_coeff_Gy);
            } else if (er_model == ER_ButtsKatz) { // "old" Katz RDD
                Katz_plateau_Gy = AT_RDD_Katz_LinearER_Daverage_Gy(KatzPoint_r_min_m, r_max_m, max_electron_range_m,
                                                                   Katz_point_coeff_Gy);
            }

            inactivation_cross_section_m2[i] = AT_KatzModel_KatzExtTarget_inactivation_cross_section_m2(
                    a0_m,
                    KatzPoint_r_min_m,
                    max_electron_range_m,
                    er_model,
                    alpha,
                    Katz_plateau_Gy,
                    Katz_point_coeff_Gy,
                    D0_Gy,
                    c,
                    m);
        }

        return EXIT_SUCCESS;
    }

    if (rdd_model == RDD_CucinottaExtTarget) {
        long i;
        const double density_g_cm3 = AT_density_g_cm3_from_material_no(Water_Liquid);
        const double density_kg_m3 = density_g_cm3 * 1000.0;

        const double a0_m = rdd_parameters[1]; // AT_RDD_a0_m( max_electron_range_m, rdd_model, rdd_parameters );

        for (i = 0; i < n; i++) {

            const double max_electron_range_m = AT_max_electron_range_m(E_MeV_u[i], Water_Liquid, (int) er_model);
            const double KatzPoint_r_min_m = AT_RDD_r_min_m(max_electron_range_m, rdd_model, rdd_parameters);
            const double Katz_point_coeff_Gy = AT_RDD_Katz_coeff_Gy_general(E_MeV_u[i], particle_no, Water_Liquid,
                                                                            er_model);
            const double r_max_m = GSL_MIN(a0_m, max_electron_range_m);
            double LET_MeV_cm2_g;
            AT_Mass_Stopping_Power_with_no(stop_power_source,
                                           1,
                                           &E_MeV_u[i],
                                           &particle_no,
                                           Water_Liquid,
                                           &LET_MeV_cm2_g);
            const double LET_J_m = LET_MeV_cm2_g * density_g_cm3 * 100.0 * MeV_to_J; // [MeV / cm] -> [J/m]
            const double beta = AT_beta_from_E_single(E_MeV_u[i]);
            const double C_norm = AT_RDD_Cucinotta_Cnorm(KatzPoint_r_min_m, max_electron_range_m, beta, density_kg_m3,
                                                         LET_J_m, Katz_point_coeff_Gy);
            double Cucinotta_plateau_Gy = AT_RDD_Cucinotta_Ddelta_average_Gy(KatzPoint_r_min_m, r_max_m,
                                                                             max_electron_range_m, beta,
                                                                             Katz_point_coeff_Gy);
            Cucinotta_plateau_Gy += C_norm *
                                    AT_RDD_Cucinotta_Dexc_average_Gy(KatzPoint_r_min_m, r_max_m, max_electron_range_m,
                                                                     beta, Katz_point_coeff_Gy);

            inactivation_cross_section_m2[i] = AT_KatzModel_CucinottaExtTarget_inactivation_cross_section_m2(
                    a0_m,
                    KatzPoint_r_min_m,
                    max_electron_range_m,
                    beta,
                    C_norm,
                    Cucinotta_plateau_Gy,
                    Katz_point_coeff_Gy,
                    D0_Gy,
                    c,
                    m);
        }

        return EXIT_SUCCESS;
    }

#ifndef NDEBUG
    char rdd_name[200];
    AT_RDD_name_from_number(rdd_model, rdd_name);
    fprintf(stderr, "RDD model %ld [%s] not supported\n", rdd_model, rdd_name);
#endif
    return EXIT_FAILURE;
}


double AT_KatzModel_KatzExtTarget_ButtsKatz_TrackWidth(
		const double z2kappabeta2,
		const double m) {

	const double    am[6] = { 1.5,   2.0,  2.5,   3.0,   4.0,   5.0};
	const double conx1[6] = { 4.314, 4.6,  4.82,  5.0,   5.29,  5.51};
	const double cony1[6] = { 6.55,  5.25, 4.7,   4.6,   4.4,   4.33};
	const double conx2[6] = { 4.4,   6.3,  8.3,  10.3,  12.8,  14.3};
	const double cony2[6] = { 6.65,  6.05, 5.95,  5.85,  5.75,  5.65};

	int i;
	double factor,ylo,yhi,yval,con,track;

	if ( (m<1.5) || (m>5.0) ){
		return -1;
	} else{

		i = 1;
		while (m >= am[i]) {
			i++;
		}
		i--;

		factor = (m - am[i])/(am[i+1]-am[i]);
		if (z2kappabeta2 < conx1[i+1]){
			return 1.0;
		} else if (z2kappabeta2 >= conx2[i+1]){
			//    for z2/kb2 high enough that both y values for low & high
			//    m can be calculated using a translated y=(1-exp(-z2/kb2))
			ylo = (cony2[i]*7.7687e-1)/(1-exp(-(conx2[i]*1.5e0)/z2kappabeta2));
			yhi = (cony2[i+1]*7.7687e-1)/(1-exp(-(conx2[i+1]*1.5e0)/z2kappabeta2));
		} else if (z2kappabeta2 > conx2[i]){
			//    for values where y for low m must be calculated using trans-
			//    lated y=(1-exp(-z2/kb2)) and high m y=straight line approximation
			ylo = (cony2[i]*7.7687e-1)/(1-exp(-(conx2[i]*1.5e0)/z2kappabeta2));
			yhi = (((cony2[i+1]-cony1[i+1])/(conx2[i+1]-conx1[i+1]))*(z2kappabeta2-conx1[i+1]))+cony1[i+1];
		} else {
			//    if z2/kb2 does not fit into any of above intervals, then it must
			//    be in region where both curves are approximated by straight lines
			//	 std::cout << "both straight lines" << std::endl;
			ylo = (((cony2[i]-cony1[i])/(conx2[i]-conx1[i]))*(z2kappabeta2-conx1[i]))+cony1[i];
			yhi = (((cony2[i+1]-cony1[i+1])/(conx2[i+1]-conx1[i+1]))*(z2kappabeta2-conx1[i+1]))+cony1[i+i];
		}

		//    having found y values for low and high m, perform interpolation
		//    on the logs of these values to find the value for the unknown m
		//    an interpolation is also done on the z2/kb2 constant which is the
		//    start of the track width region, to find this for the unknown m
		yval = exp(log(ylo)+(factor*(log(yhi)-log(ylo))));
		con = exp(log(cony1[i])+(factor*(log(cony1[i+1])-log(cony1[i]))));
		track = yval/con;
	}
	return track;
}


double AT_KatzModel_KatzExtTarget_Zhang_TrackWidth(
		const double z2kappabeta2,
		const double m) {

	const double XA = -log(1.0-pow(0.98,1.0/m));
	const double YA = 0.98;

	const double gi[6] = { -0.223198, 6.79973, -9.08647, 7.62573, -3.42814,  0.632902};
	const double hi[6] = { -1.03134, 5.06363, -9.39776, 9.3199, -4.75329, 0.97801};

	const double logm = log(m);

	const double XB = exp( gi[0] + gi[1]*logm + gi[2]*logm*logm + gi[3]*pow(logm,3) + gi[4]*pow(logm,4) + gi[5]*pow(logm,5));
	const double YB = exp( hi[0] + hi[1]*logm + hi[2]*logm*logm + hi[3]*pow(logm,3) + hi[4]*pow(logm,4) + hi[5]*pow(logm,5));

	double result = -1;
	if ( (m<1.5) || (m>3.5) ){
		return -1;
	} else{
		if( z2kappabeta2 >= XB ){
			result = YB * 0.8209 / ( 1.0 - exp(-XB * 1.72 / z2kappabeta2 ));
		} else {
			result = ((YA-YB)/(XA-XB))*(z2kappabeta2 - XA) + YA;
		}
	}
	return result;
}


double AT_KatzModel_inactivation_cross_section_approximation_m2(
		const double E_MeV_u,
		const long   particle_no,
		const long   rdd_model,
		const long   er_model,
		const double m_number_of_targets,
		const double sigma0_m2,
		const double kappa){

	double result = -1.0;

	double beta = AT_beta_from_E_single(E_MeV_u);
	double zeff = AT_effective_charge_from_beta_single(beta, AT_Z_from_particle_no_single(particle_no));
	double z2kappabeta2 = gsl_pow_2( zeff / beta ) / kappa;

	double Pi = pow( 1.0 - exp( - z2kappabeta2) , m_number_of_targets);

	if( rdd_model == RDD_KatzExtTarget && er_model == ER_ButtsKatz ){

		double factor = 1;
		if( Pi > 0.98 ){
			factor = AT_KatzModel_KatzExtTarget_ButtsKatz_TrackWidth( z2kappabeta2, m_number_of_targets );
		} else {
			factor = Pi;
		}
		result = factor * sigma0_m2;
	}

	if( rdd_model == RDD_KatzExtTarget && ((er_model == ER_Waligorski) || (er_model == ER_Edmund)) ){

		double factor = 1;
		if( Pi > 0.98 ){
			factor = AT_KatzModel_KatzExtTarget_Zhang_TrackWidth( z2kappabeta2, m_number_of_targets );
		} else {
			factor = Pi;
		}
		result = factor * sigma0_m2;
	}

	return result;
}



double AT_KatzModel_single_field_survival_from_sigma(
    const double dose_Gy,
	const double E_MeV_u,
    const long   particle_no,
    const double sigma_m2,
    const double D0_Gy,
    const double m,
    const double sigma0_m2,
    const long   stopping_power_source_no){

	double fluence_cm2 = AT_fluence_cm2_from_dose_Gy_single( E_MeV_u, particle_no, dose_Gy, Water_Liquid, stopping_power_source_no ); /* fluence * LET -> dose * density */

	assert(sigma0_m2 > 0);
	assert(sigma_m2 > 0);

	/* fraction of dose delivered in ion kill mode */
	double ion_kill_mode_fraction = sigma_m2 / sigma0_m2;
	if( ion_kill_mode_fraction > 1.0) ion_kill_mode_fraction = 1.0;

	double gamma_kill_dose = (1. - ion_kill_mode_fraction) * dose_Gy;


	double gamma_kill_mode_survival;
	double ion_kill_mode_survival;

	/* ion kill mode survival gives exponential SF curve */
	ion_kill_mode_survival = exp(- sigma_m2 * 1e4 * fluence_cm2);  /* depends on D0, m and a0; not on kappa !*/

	assert(D0_Gy > 0 );

	/* gamma kill mode survival gives shouldered SF curve */
	if( ion_kill_mode_fraction > 0.98 ){
		gamma_kill_mode_survival = 1.0;
		/* this if is quite artificial, for ion_kill_mode_fraction greater than 0.98
		 * gamma_kill_mode_survival is so close to 1 that it does not make sense to calculate it
		 */
	} else {
		gamma_kill_mode_survival = 1. - pow( 1. - exp(- gamma_kill_dose / D0_Gy), m);  /* depends on D0, m and kappa; not on a0 ! */
	}

	/* finally survival as a product of ion and gamma kill modes */
	return ion_kill_mode_survival * gamma_kill_mode_survival;
}



/* TODO implement old Katz with kappa and sigma instead of track-width here */
double AT_KatzModel_single_field_survival(
    const double dose_Gy,
	const double E_MeV_u,
    const long   particle_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double D0_Gy,
    const double m,
    const double sigma0_m2,
    const bool   use_approximation,
    const double kappa,
    const long   stopping_power_source_no){

	assert( sigma0_m2 > 0);

	/* single particle inactivation cross section calculation */
	double inactivation_cross_section_m2 = 0.0;
	double gamma_parameters[5] = {1., D0_Gy, 1., m, 0.};
	double survival = -1.0;

	if(  use_approximation == false ){
		int status = AT_KatzModel_inactivation_cross_section_m2(
				1,
				&E_MeV_u,
				particle_no,
				rdd_model,
				rdd_parameters,
				er_model,
				gamma_parameters,
				stopping_power_source_no,
				&inactivation_cross_section_m2);    /* here we use D0, m and a0 */

		if( status != EXIT_SUCCESS ){
#ifndef NDEBUG
			fprintf(stderr, "Problem with evaluating inactivation cross section\n");
#endif
			return -1.0;
		}
	} else {
		inactivation_cross_section_m2 = AT_KatzModel_inactivation_cross_section_approximation_m2(
                E_MeV_u,
                particle_no,
                rdd_model,
                er_model,
                m,
                sigma0_m2,
                kappa);    /* here we use sigma0, m and kappa */
	}


    survival = AT_KatzModel_single_field_survival_from_sigma(dose_Gy,
                                                             E_MeV_u,
                                                             particle_no,
                                                             inactivation_cross_section_m2,
                                                             D0_Gy,
                                                             m,
                                                             sigma0_m2,
                                                             stopping_power_source_no);

	return survival;
}


typedef struct {
    double E_MeV_u;
    long particle_no;
    long rdd_model;
    double *rdd_parameters;
    long er_model;
    double D0_Gy;
    double m;
    double sigma0_m2;
    bool use_approximation;
    double kappa;
    long stopping_power_source_no;
    double survival;
} _AT_KatzModel_survival_root_params;


double _AT_KatzModel_survival_root(double x, void *params) {
    assert(params != NULL);

    _AT_KatzModel_survival_root_params *current_params = (_AT_KatzModel_survival_root_params *) params;

    double result = AT_KatzModel_single_field_survival(x,
                                                       current_params->E_MeV_u,
                                                       current_params->particle_no,
                                                       current_params->rdd_model,
                                                       current_params->rdd_parameters,
                                                       current_params->er_model,
                                                       current_params->D0_Gy,
                                                       current_params->m,
                                                       current_params->sigma0_m2,
                                                       current_params->use_approximation,
                                                       current_params->kappa,
                                                       current_params->stopping_power_source_no);

    result -= current_params->survival;

    return result;
}


double AT_KatzModel_single_field_rbe(const double E_MeV_u,
                                     const long particle_no,
                                     const long rdd_model,
                                     const double rdd_parameters[],
                                     const long er_model,
                                     const double D0_Gy,
                                     const double m,
                                     const double sigma0_m2,
                                     const bool use_approximation,
                                     const double kappa,
                                     const long stopping_power_source_no,
                                     const double survival) {

    // set default RBE to negative value
    double rbe = -1;

    // set default survival to 0.1
    double current_survival = survival;
    if (current_survival < 0.0)
        current_survival = 0.1;

    // calculate dose for which reference radiation predicts given survival level
    // formula below can be derived from following equation (mtarget model with m,D0 parameters ):
    //   current_survival = 1 - (1 - exp( - D_ref_Gy / D0_Gy))^m
    double D_ref_Gy = -D0_Gy * log(1. - pow(1. - current_survival, 1. / m));

    // calculate dose for which irradiation with ions predicts given survival level
    // it is impossible to analytically inverse Katz model formulation, therefore we use GSL numerical solver
    _AT_KatzModel_survival_root_params current_params;
    current_params.E_MeV_u = E_MeV_u;
    current_params.particle_no = particle_no;
    current_params.rdd_model = rdd_model;
    current_params.rdd_parameters = (double *) rdd_parameters;
    current_params.er_model = er_model;
    current_params.D0_Gy = D0_Gy;
    current_params.m = m;
    current_params.sigma0_m2 = sigma0_m2;
    current_params.use_approximation = use_approximation;
    current_params.kappa = kappa;
    current_params.stopping_power_source_no = stopping_power_source_no;
    current_params.survival = current_survival;

    int status;
    int iter = 0, max_iter = 1000;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    gsl_function F;

    F.function = &_AT_KatzModel_survival_root; // function to solve for (D): f(D) = SF_ion(D) - SF_level
    F.params = (void *) (&current_params);

    double x_root = D_ref_Gy; // solver initial guess

    // solver lower and upper bounds
    double x_lower = 0.0;
    double x_upper = 10000 * D_ref_Gy;

    if (GSL_FN_EVAL(&F, x_lower) * GSL_FN_EVAL(&F, x_upper) > 0) {
        return -1.0;
    }

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, x_lower, x_upper);

    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        x_root = gsl_root_fsolver_root(s);
        x_lower = gsl_root_fsolver_x_lower(s);
        x_upper = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_delta(x_lower, x_upper, 0, 1e-6);
//        printf("iter %d, lower %g, upper %g, status %d\n", iter, x_lower, x_upper, status);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    double D_ion_Gy = x_root;

//    printf("D_ion = %g [Gy]\n", D_ion_Gy);

    rbe = D_ref_Gy / D_ion_Gy;

    return rbe;
}


int AT_KatzModel_mixed_field_survival(
	const long   number_of_items,
    double fluence_cm2[],
	const double E_MeV_u[],
    const long   particle_no[],
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double D0_Gy,
    const double m,
    const double sigma0_m2,
    const bool   use_approximation,
    const double kappa,
    const long   stopping_power_source_no,
    double * survival){

	assert( sigma0_m2 > 0);

	/* single particle inactivation cross section calculation */
	double inactivation_cross_section_m2 = 0.0;
	double gamma_parameters[5] = {1., D0_Gy, 1., m, 0.};

	long i = 0;
	double sum1 = 0.0;
	double sum2 = 0.0;

#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(+:sum1,sum2) private(inactivation_cross_section_m2)
#endif
	for( i = 0 ; i < number_of_items; i++){

#ifdef HAVE_OPENMP
		if( fluence_cm2[i] == 0 || E_MeV_u[i] == 0){
			inactivation_cross_section_m2 = 0.0;
		} else {
			inactivation_cross_section_m2 = AT_KatzModel_inactivation_cross_section_approximation_m2(E_MeV_u[i], particle_no[i], rdd_model, er_model, m, sigma0_m2, kappa );
		}
#else
		if( use_approximation ){
			inactivation_cross_section_m2 = AT_KatzModel_inactivation_cross_section_approximation_m2(E_MeV_u[i], particle_no[i], rdd_model, er_model, m, sigma0_m2, kappa );
		} else {
			int status = AT_KatzModel_inactivation_cross_section_m2(
					1,
					&(E_MeV_u[i]),
					particle_no[i],
					rdd_model,
					rdd_parameters,
					er_model,
					gamma_parameters,
					stopping_power_source_no,
					&inactivation_cross_section_m2);    /* here we use D0, m and a0 */

			if( status != EXIT_SUCCESS ){
#ifndef NDEBUG
				fprintf(stderr, "Problem with evaluating inactivation cross section\n");
#endif
				return status;
			}
		}
#endif

		double Pi = 1.0;
		if( inactivation_cross_section_m2 < sigma0_m2 ){
			Pi = inactivation_cross_section_m2 / sigma0_m2;
		};

		double Di = 0.0;
		if( fluence_cm2[i] > 0 && E_MeV_u[i] > 0){
			Di = AT_dose_Gy_from_fluence_cm2_single( E_MeV_u[i], particle_no[i], fluence_cm2[i], Water_Liquid, stopping_power_source_no);
		}

		sum1 += inactivation_cross_section_m2 * fluence_cm2[i] * 1e4;
		sum2 += (1.0 - Pi) * Di;

	}

	double IonKill = exp( - sum1 );
	double GammaKill = 1.0 - pow( 1 - exp(- sum2 / D0_Gy), m);

	*survival = IonKill * GammaKill;

	return EXIT_SUCCESS;
}



int AT_KatzModel_single_field_survival_optimized_for_fluence_vector(
	const long   number_of_items,
    const double dose_Gy[],
	const double E_MeV_u,
    const long   particle_no,
    const long   rdd_model,
    const double rdd_parameters[],
    const long   er_model,
    const double D0_Gy,
    const double m,
    const double sigma0_m2,
    const bool   use_approximation,
    const double kappa,
    const long stopping_power_source_no,
    double * survival){

	assert( sigma0_m2 > 0);

	/* single particle inactivation cross section calculation */
	double sigma_m2 = 0.0;
	double gamma_parameters[5] = {1., D0_Gy, 1., m, 0.};

	if(  use_approximation == false ){
		int status = AT_KatzModel_inactivation_cross_section_m2(
				1,
				&E_MeV_u,
				particle_no,
				rdd_model,
				rdd_parameters,
				er_model,
				gamma_parameters,
				stopping_power_source_no,
				&sigma_m2);    /* here we use D0, m and a0 */

#ifndef NDEBUG
		printf("Inactivation cross section = %g\n", sigma_m2);
#endif

		if( status != EXIT_SUCCESS ){
#ifndef NDEBUG
			fprintf(stderr, "Problem with evaluating inactivation cross section\n");
#endif
			return status;
		}
	} else {
        sigma_m2 = AT_KatzModel_inactivation_cross_section_approximation_m2(
                E_MeV_u,
                particle_no,
                rdd_model,
                er_model,
                m,
                sigma0_m2,
                kappa);    /* here we use sigma0, m and kappa */
	}

	long i;
	for( i = 0 ; i < number_of_items ; i++){
		survival[i] = AT_KatzModel_single_field_survival_from_sigma(dose_Gy[i],
                                                                    E_MeV_u,
                                                                    particle_no,
                                                                    sigma_m2,
                                                                    D0_Gy,
                                                                    m,
                                                                    sigma0_m2,
                                                                    stopping_power_source_no);
	}

	return EXIT_SUCCESS;
}


double AT_sI_int( double  r_m,
    void*   params){
  assert( params != NULL );
  double  P = AT_P_RDD(r_m, params);
  return (r_m * P);
}


double AT_P_RDD( double  r_m,
    void*   params)
{
  double  D_Gy;
  long   n_tmp              = 1;
  assert( params != NULL );
  AT_P_RDD_parameters* par  = (AT_P_RDD_parameters*)params;
  double  fr_m               = r_m;

  // Call RDD to get dose
  AT_D_RDD_Gy  (  n_tmp,
      &fr_m,
      *(par->E_MeV_u),
      *(par->particle_no),
      *(par->material_no),
      *(par->rdd_model),
      par->rdd_parameters,
      *(par->er_model),
      PSTAR,
      &D_Gy);

  long gamma_model = GR_GeneralTarget;
  double P;
  AT_gamma_response(  n_tmp,
      &D_Gy,
      gamma_model,
      par->gamma_parameters,
      false,
      // return
      &P);
  return P;
}
