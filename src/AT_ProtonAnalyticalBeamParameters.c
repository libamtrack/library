/**
 * @brief Proton analytical models of dose, LET and RBE
 */

/*
 *    AT_ProtonAnalyticalBeamParameters.c
 *    ==============
 *
 *    Created on: 11.02.2019
 *    Creator: grzanka
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

#include "AT_ProtonAnalyticalBeamParameters.h"

typedef struct {
    double E_MeV_u;
    double sigma_E_MeV_u;
    long material_no;
    double eps;
} _AT_dose_Bortfeld_Gy_negative_params;


double _AT_dose_Bortfeld_Gy_negative(double x, void *params) {
    assert(params != NULL);

    _AT_dose_Bortfeld_Gy_negative_params *current_params = (_AT_dose_Bortfeld_Gy_negative_params *) params;

    double result = -AT_dose_Bortfeld_Gy_single(x,
                                                1e9,
                                                current_params->E_MeV_u,
                                                current_params->sigma_E_MeV_u,
                                                current_params->material_no,
                                                current_params->eps);


    return result;
}


double AT_max_location_Bortfeld_cm(const double E_MeV_u,
                                   const double sigma_E_MeV_u,
                                   const long material_no,
                                   const double eps) {

    double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])

    double max_location_guess_cm = alpha * pow(E_MeV_u, p);  // range in [cm]

    _AT_dose_Bortfeld_Gy_negative_params current_params;
    current_params.E_MeV_u = E_MeV_u;
    current_params.sigma_E_MeV_u = sigma_E_MeV_u;
    current_params.material_no = material_no;
    current_params.eps = eps;

    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;

    double x_minimum = max_location_guess_cm;
    double x_lower = 0.0;
    double x_upper = 2.0 * max_location_guess_cm;

    gsl_function F;

    F.function = &_AT_dose_Bortfeld_Gy_negative;
    F.params = (void *) (&current_params);

    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc(T);
    gsl_min_fminimizer_set(s, &F, x_minimum, x_lower, x_upper);

//    printf("using %s method\n",
//           gsl_min_fminimizer_name(s));
//
//    printf("%5s [%9s, %9s] %9s %10s %9s\n",
//           "iter", "lower", "upper", "min",
//           "err", "err(est)");
//
//    printf("%5d [%.7f, %.7f] %.7f %.7f\n",
//           iter, x_lower, x_upper,
//           x_minimum, x_upper - x_lower);

    do {
        iter++;
        status = gsl_min_fminimizer_iterate(s);

        x_minimum = gsl_min_fminimizer_x_minimum(s);
        x_lower = gsl_min_fminimizer_x_lower(s);
        x_upper = gsl_min_fminimizer_x_upper(s);

        status = gsl_min_test_interval(x_lower, x_upper, 1e-6, 0.0);

//        if (status == GSL_SUCCESS)
//            printf("Converged:\n");
//
//        printf("%5d [%.7f, %.7f] "
//               "%.7f %.7f\n",
//               iter, x_lower, x_upper,
//               x_minimum, x_upper - x_lower);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_min_fminimizer_free(s);

    return x_minimum;

}

typedef struct {
    double E_MeV_u;
    double sigma_E_MeV_u;
    long material_no;
    double dose_cut_Gy;
    double eps;
} _AT_dose_Bortfeld_Gy_root_params;

double _AT_dose_Bortfeld_Gy_root(double x, void *params) {
    assert(params != NULL);

    _AT_dose_Bortfeld_Gy_root_params *current_params = (_AT_dose_Bortfeld_Gy_root_params *) params;

    double result = AT_dose_Bortfeld_Gy_single(-x,
                                               1e9,
                                               current_params->E_MeV_u,
                                               current_params->sigma_E_MeV_u,
                                               current_params->material_no,
                                               current_params->eps);

    result -= current_params->dose_cut_Gy;

//    printf("f(%g) = %g\n", x, result);

    return result;
}


double AT_range_Bortfeld_cm(const double E_MeV_u,
                            const double sigma_E_MeV_u,
                            const long material_no,
                            const double eps,
                            const double dose_drop,
                            const short search_direction) {

    double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])

    double range_guess_cm = alpha * pow(E_MeV_u, p);  // range in [cm]

    double current_dose_drop = dose_drop;
    if (current_dose_drop < 0.0)
        current_dose_drop = 0.8;

    _AT_dose_Bortfeld_Gy_root_params current_params;
    current_params.E_MeV_u = E_MeV_u;
    current_params.sigma_E_MeV_u = sigma_E_MeV_u;
    current_params.material_no = material_no;
    current_params.eps = eps;
    current_params.dose_cut_Gy = 0.0;

    double max_location_cm = AT_max_location_Bortfeld_cm(E_MeV_u,
                                                         sigma_E_MeV_u,
                                                         material_no,
                                                         eps);

    current_params.dose_cut_Gy =
            current_dose_drop * _AT_dose_Bortfeld_Gy_root(-max_location_cm, (void *) (&current_params));

    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    gsl_function F;

    F.function = &_AT_dose_Bortfeld_Gy_root;
    F.params = (void *) (&current_params);


    double x_root = -range_guess_cm;
    double x_lower = -2.0 * range_guess_cm;
    double x_upper = -max_location_cm;

    if (search_direction < 0) {
        x_lower = -max_location_cm;
        x_upper = 0.0;
        x_root = 0.5 * (x_lower + x_upper);

        if (GSL_FN_EVAL(&F, x_lower) * GSL_FN_EVAL(&F, x_upper) > 0) {
            return 0.0;
        }
    }


    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, x_lower, x_upper);

//    printf ("using %s method\n",
//            gsl_root_fsolver_name (s));
//
//    printf ("%5s [%9s, %9s] %9s %10s %9s\n",
//            "iter", "lower", "upper", "root",
//            "err", "err(est)");
//
//    printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//            iter, x_lower, x_upper,
//            x_root, x_upper - x_lower);

    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        x_root = gsl_root_fsolver_root(s);
        x_lower = gsl_root_fsolver_x_lower(s);
        x_upper = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lower, x_upper,
                                        0, 1e-6);

//        if (status == GSL_SUCCESS)
//            printf ("Converged:\n");
//
//        printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//                iter, x_lower, x_upper,
//                x_root, x_upper - x_lower);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return -x_root;

}


double AT_fwhm_Bortfeld_cm(const double E_MeV_u,
                           const double sigma_E_MeV_u,
                           const long material_no,
                           const double eps) {

    double distal_50 = AT_range_Bortfeld_cm(E_MeV_u,
                                            sigma_E_MeV_u,
                                            material_no,
                                            eps,
                                            0.5,
                                            1);

    double proximal_50 = AT_range_Bortfeld_cm(E_MeV_u,
                                              sigma_E_MeV_u,
                                              material_no,
                                              eps,
                                              0.5,
                                              -1);

    double result;
    if (proximal_50 == 0) {
        result = 1.0 / 0.0;
    } else {
        result = distal_50 - proximal_50;
    }

    return result;
}


double AT_max_plateau_Bortfeld(const double E_MeV_u,
                               const double sigma_E_MeV_u,
                               const long material_no,
                               const double eps) {

    double entrance_dose_Gy = AT_dose_Bortfeld_Gy_single(0.0,
                                                         1e8,
                                                         E_MeV_u,
                                                         sigma_E_MeV_u,
                                                         material_no,
                                                         eps);

    double max_position_cm = AT_max_location_Bortfeld_cm(E_MeV_u,
                                                         sigma_E_MeV_u,
                                                         material_no,
                                                         eps);

    double max_dose_Gy = AT_dose_Bortfeld_Gy_single(max_position_cm,
                                                    1e8,
                                                    E_MeV_u,
                                                    sigma_E_MeV_u,
                                                    material_no,
                                                    eps);

    return max_dose_Gy / entrance_dose_Gy;

}

/**
 * TODO
 */
typedef struct {
    double range_cm;
    double sigma_E_MeV_u;
    long material_no;
    double dose_drop;
    double eps;
} _AT_range_Bortfeld_Gy_root_params;


double _AT_range_Bortfeld_Gy(double x, void *params) {
    assert(params != NULL);

    _AT_range_Bortfeld_Gy_root_params *current_params = (_AT_range_Bortfeld_Gy_root_params *) params;


    double result = AT_range_Bortfeld_cm(x,
                                         current_params->sigma_E_MeV_u,
                                         current_params->material_no,
                                         current_params->eps,
                                         current_params->dose_drop,
                                         1);

    result -= current_params->range_cm;

    return result;
}


double AT_energy_Bortfeld_MeV_u(const double range_cm,
                                const double sigma_E_MeV_u,
                                const long material_no,
                                const double eps,
                                const double dose_drop) {

    if (range_cm == 0) {
        return 0.0;
    }

    double p = AT_p_MeV_from_material_no(material_no); // exponent of range-energy relation
    double alpha = AT_alpha_g_cm2_MeV_from_material_no(material_no); // proportionality factor (0.0022 cm/MeV^p in [1])

    double E_guess_MeV_u = pow(range_cm / alpha, 1.0 / p);

    double current_dose_drop = dose_drop;
    if (current_dose_drop < 0.0)
        current_dose_drop = 0.8;

    _AT_range_Bortfeld_Gy_root_params current_params;
    current_params.range_cm = range_cm;
    current_params.sigma_E_MeV_u = sigma_E_MeV_u;
    current_params.material_no = material_no;
    current_params.eps = eps;
    current_params.dose_drop = dose_drop;

    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double x_root = E_guess_MeV_u;
    double x_lower = 0.5 * E_guess_MeV_u;
    double x_upper = 2.0 * E_guess_MeV_u;
    gsl_function F;

    F.function = &_AT_range_Bortfeld_Gy;
    F.params = (void *) (&current_params);

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, x_lower, x_upper);

//    printf("using %s method\n",
//           gsl_root_fsolver_name(s));
//
//    printf("%5s [%9s, %9s] %9s %10s %9s\n",
//           "iter", "lower", "upper", "root",
//           "err", "err(est)");
//
//    printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//           iter, x_lower, x_upper,
//           x_root, x_upper - x_lower);

    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        x_root = gsl_root_fsolver_root(s);
        x_lower = gsl_root_fsolver_x_lower(s);
        x_upper = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lower, x_upper,
                                        0, 1e-6);

//        if (status == GSL_SUCCESS)
//            printf("Converged:\n");
//
//        printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//               iter, x_lower, x_upper,
//               x_root, x_upper - x_lower);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return x_root;

}


/**
 * TODO
 */
typedef struct {
    double range_cm;
    double fwhm_cm;
    double max_plateau;
    long material_no;
    double eps;
    double dose_drop;
} _AT_func_params;


// real -> [0, 0.2]
double constraint(double x) {
    if (x > 10) {
        return 0.2 * (M_1_PI * atan(10.0) + 0.5);
    } else {
        return 0.2 * (M_1_PI * atan(x) + 0.5);
    }
}

// [0, 0.2] -> real
double deconstraint(double y) {
    return tan(M_PI * (y / 0.2 - 0.5));
}


int
func_f(const gsl_vector *x, void *params, gsl_vector *f) {
    double E_MeV_u = gsl_vector_get(x, 0);
    double sigma_E_MeV_u = gsl_vector_get(x, 1);
    double eps = constraint(gsl_vector_get(x, 2));

    _AT_func_params *current_params = (_AT_func_params *) params;

    double current_range_cm = AT_range_Bortfeld_cm(E_MeV_u,
                                                   sigma_E_MeV_u,
                                                   current_params->material_no,
                                                   eps,
                                                   current_params->dose_drop,
                                                   1);

    double current_fwhm_cm = AT_fwhm_Bortfeld_cm(E_MeV_u,
                                                 sigma_E_MeV_u,
                                                 current_params->material_no,
                                                 eps);

    double current_max_plateau = AT_max_plateau_Bortfeld(E_MeV_u,
                                                         sigma_E_MeV_u,
                                                         current_params->material_no,
                                                         eps);

    gsl_vector_set(f, 0, 10.0 * (current_range_cm - (current_params->range_cm)) / (current_params->range_cm));
    gsl_vector_set(f, 1, 1.0 * (current_fwhm_cm - (current_params->fwhm_cm)) / (current_params->fwhm_cm));
    gsl_vector_set(f, 2, 1.0 * (current_max_plateau - (current_params->max_plateau)) / (current_params->max_plateau));

    return GSL_SUCCESS;
}

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w) {
    gsl_vector *x = gsl_multifit_nlinear_position(w);

    /* print out current location */
    printf("%f %f %f\t",
           gsl_vector_get(x, 0),
           gsl_vector_get(x, 1),
           constraint(gsl_vector_get(x, 2)));

    _AT_func_params *current_params = (_AT_func_params *) (params);


    double current_range_cm = AT_range_Bortfeld_cm(gsl_vector_get(x, 0),
                                                   gsl_vector_get(x, 1),
                                                   current_params->material_no,
                                                   constraint(gsl_vector_get(x, 2)),
                                                   current_params->dose_drop,
                                                   1);

    double current_fwhm_cm = AT_fwhm_Bortfeld_cm(gsl_vector_get(x, 0),
                                                 gsl_vector_get(x, 1),
                                                 current_params->material_no,
                                                 constraint(gsl_vector_get(x, 2)));

    double current_max_plateau = AT_max_plateau_Bortfeld(gsl_vector_get(x, 0),
                                                         gsl_vector_get(x, 1),
                                                         current_params->material_no,
                                                         constraint(gsl_vector_get(x, 2)));

    /* print out current location */
    printf("range = %g , fwhm = %g , max/plat = %g\n", current_range_cm, current_fwhm_cm, current_max_plateau);

}

void AT_fit_Bortfeld(const double range_cm,
                     const double fwhm_cm,
                     const double max_to_plateau,
                     const long material_no,
                     const double dose_drop) {

    printf("TARGER range = %g , fwhm = %g , max/plat = %g\n", range_cm, fwhm_cm, max_to_plateau);

    /* problem dimensions */
    const size_t n = 3;
    const size_t p = 3;

    const size_t max_iter = 200;
    const double xtol = 1.0e-8;
    const double gtol = 1.0e-8;
    const double ftol = 1.0e-8;

    /* starting point */
    gsl_vector *x0 = gsl_vector_alloc(p);
    gsl_vector_set(x0, 0, 100.0);
    gsl_vector_set(x0, 1, 1.5);
    gsl_vector_set(x0, 2, deconstraint(0.02));


    /* define function to be minimized */
    _AT_func_params current_params;
    current_params.range_cm = range_cm;
    current_params.fwhm_cm = fwhm_cm;
    current_params.max_plateau = max_to_plateau;
    current_params.material_no = material_no;
    current_params.eps = -1;
    current_params.dose_drop = dose_drop;

    gsl_multifit_nlinear_fdf fdf;
    fdf.f = func_f;
    fdf.df = NULL;
    fdf.fvv = NULL;
    fdf.n = n;
    fdf.p = p;
    fdf.params = (void *) (&current_params); // TODO

    /* output */
    int info; // reasons of convergence

    gsl_multifit_nlinear_parameters params = gsl_multifit_nlinear_default_parameters();
    params.trs = gsl_multifit_nlinear_trs_lmaccel;

    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *work =
            gsl_multifit_nlinear_alloc(T, &params, n, p);

    gsl_vector *f = gsl_multifit_nlinear_residual(work);
//    gsl_vector * x = gsl_multifit_nlinear_position(work);
    double chisq0, chisq, rcond;

    /* initialize solver */
    gsl_multifit_nlinear_init(x0, &fdf, work);

    /* store initial cost */
    gsl_blas_ddot(f, f, &chisq0);

    /* iterate until convergence */
    gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                                callback, (void *) (&current_params), &info, work);

    /* compute covariance of best fit parameters */
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc(p, p);
    J = gsl_multifit_nlinear_jac(work);
    gsl_multifit_nlinear_covar(J, 0.0, covar);

//    /* store final cost */
    gsl_blas_ddot(f, f, &chisq);
//
//    /* store cond(J(x)) */
//    gsl_multifit_nlinear_rcond(&rcond, work);
//
//    /* print summary */

    double dof = n - p;
    double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

    fprintf(stderr, "E       = %.5f +/- %.5f\n", gsl_vector_get(work->x, 0), c * sqrt(gsl_matrix_get(covar, 0, 0)));
    fprintf(stderr, "delta E = %.5f +/- %.5f\n", gsl_vector_get(work->x, 1), c * sqrt(gsl_matrix_get(covar, 1, 1)));
    fprintf(stderr, "eps     = %.5f +/- %.5f\n", constraint(gsl_vector_get(work->x, 2)),
            c * sqrt(gsl_matrix_get(covar, 2, 2)));

//
    fprintf(stderr, "NITER         = %zu\n", gsl_multifit_nlinear_niter(work));
    fprintf(stderr, "NFEV          = %zu\n", fdf.nevalf);
    fprintf(stderr, "NJEV          = %zu\n", fdf.nevaldf);
    fprintf(stderr, "NAEV          = %zu\n", fdf.nevalfvv);
    fprintf(stderr, "initial cost  = %.12e\n", chisq0);
    fprintf(stderr, "final cost    = %.12e\n", chisq);
//    fprintf(stderr, "final x       = (%.12e, %.12e)\n", gsl_vector_get(x, 0), gsl_vector_get(x, 1));
//    fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);
//
//    printf("\n\n");

    gsl_multifit_nlinear_free(work);

}