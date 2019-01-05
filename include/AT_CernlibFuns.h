//Mock of funcions from thirdparty/cernlib
#ifndef LIBRARY_AT_CERNLIBFUNS_H
#define LIBRARY_AT_CERNLIBFUNS_H

#include <stdio.h>

/**
 * TODO
 * @param[in]  lambda_landau
 * @return     Landau 
 */
double CL_denlan(double lambda_landau);

/**
 * TODO
 * @param[in]  rnd
 * @return     Landau random number
 */
double CL_ranlan(double rnd);

/**
 * TODO
 * @param[in]  kappa
 * @param[out]  beta2
 */
void CL_vavset(double kappa, double beta2);

/**
 * TODO
 * @param[in]  lambda_vavilov
 * @return     Vavilov 
 */
double CL_vavden(double lambda_vavilov);

/**
 * TODO
 * @param[in]  kappa
 * @param[in]  beta2
 * @param[in]  rnd
 * @return     Vavilov 
 */
double CL_vavran(double kappa, double beta2, double rnd);

#endif //LIBRARY_AT_CERNLIBFUNS_H
