#include "AT_CernlibFuns.h"

double CL_denlan(double lambda_landau){
	extern float denlan_();
	float lambda_landau_f = (float)lambda_landau;
	return((double)denlan_(&lambda_landau_f));
}

double CL_ranlan(double rnd){
	extern float ranlan_();
	float rnd_f = (float)rnd;
	return((double)ranlan_(&rnd_f));
}


void CL_vavset(double kappa, double beta2){
	float kappa_f = (float)kappa;
	float beta2_f = (float)beta2;
        extern void vavset_();
	int mode = 1;
	vavset_(&kappa_f, &beta2_f, &mode);
	return;
}

double CL_vavden(double lambda_vavilov){
	extern float vavden_();
	float lambda_vavilov_f = (float)lambda_vavilov;
	return((double)vavden_(&lambda_vavilov_f));
}

double CL_vavran(double kappa, double beta2, double rnd){
	extern float vavran_();
	float kappa_f = (float)kappa;
	float beta2_f = (float)beta2;
	float rnd_f = (float)rnd;
	return((double)vavran_(&kappa_f, &beta2_f, &rnd_f));
}

