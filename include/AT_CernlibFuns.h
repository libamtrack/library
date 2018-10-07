//
// Created by osboxes on 04.09.18.
//

#ifndef LIBRARY_AT_CERNLIBFUNS_H
#define LIBRARY_AT_CERNLIBFUNS_H

#include <stdio.h>

double CL_denlan(double lambda_landau);

double CL_ranlan(double rnd);

void CL_vavset(double kappa, double beta2);

double CL_vavden(double lambda_vavilov);

double CL_vavran(double kappa, double beta2, double rnd);

#endif //LIBRARY_AT_CERNLIBFUNS_H
