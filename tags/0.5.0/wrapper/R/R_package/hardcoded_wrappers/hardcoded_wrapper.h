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

void AT_SPC_get_size_from_filename_R(const char** filename,
    int* spc_size);

void AT_SPC_read_data_from_filename_R( const char** filename,
		const int* n,
		int* depth_step,
		double* depth_g_cm2,
		double* E_MeV_u,
		double* DE_MeV_u,
		int* particle_no,
		double* fluence_cm2,
        int* n_bins_read);

#endif /* AT_HARDCODED_WRAPPER_H_ */
