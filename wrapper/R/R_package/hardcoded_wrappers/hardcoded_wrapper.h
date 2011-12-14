#ifndef AT_HARDCODED_WRAPPER_H_
#define AT_HARDCODED_WRAPPER_H_

#include "AT_SPC.h"
#include "AT_DataParticle.h"
#include "AT_GammaResponse.h"
#include "AT_DataMaterial.h"

void AT_gamma_response_R( const int*  n,
    const float*  d_Gy,
    const int*    gamma_model,
    const float*  gamma_parameter,
    const int*	  lethal_events_mode,
    float*        S);

void AT_particle_name_from_particle_no_R(const int* particle_no,
    char** particle_name);

void AT_particle_no_from_particle_name_R(const char** particle_name,
    int* particle_no);

void AT_material_name_from_material_no_R( const int* material_no,
    char** material_name);

void AT_material_no_from_material_name_R( const char** material_name,
	int* material_no);


#endif /* AT_HARDCODED_WRAPPER_H_ */
