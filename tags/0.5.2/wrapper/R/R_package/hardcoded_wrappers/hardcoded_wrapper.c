#include "hardcoded_wrapper.h"

void AT_gamma_response_R( const int*  n,
    const float*  d_Gy,
    const int*    gamma_model,
    const float*  gamma_parameter,
    const int*    lethal_events_mode,
	float*        S){

  /* int -> long conversion */
  const long n_long = (long)(*n);
  const long gamma_model_long = (long)(*gamma_model);

  /* int -> bool conversion */
  const bool lethal_events_mode_bool = (bool)(*lethal_events_mode);

  /* float -> double conversion */
  double * d_Gy_double = (double*)calloc(n_long,sizeof(double));
  long i;
  for(i = 0 ; i < n_long ; i++){
    d_Gy_double[i] = (double)d_Gy[i];
  }
  long number_of_gamma_parameters = 0;
  if ( gamma_model_long == GR_GeneralTarget ){
    while  (gamma_parameter[number_of_gamma_parameters] != 0){
      number_of_gamma_parameters  += 4;
    }
    number_of_gamma_parameters++; /* to include also trailing zero */
  } else {
    number_of_gamma_parameters = AT_Gamma_number_of_parameters( gamma_model_long );
  }
  double * gamma_parameter_double  = (double*)calloc(number_of_gamma_parameters,sizeof(double));
  for(i = 0 ; i < number_of_gamma_parameters ; i++){
    gamma_parameter_double[i] = (double)gamma_parameter[i];
  }

  /* place for results */
  double * S_double = (double*)calloc(n_long,sizeof(double));

  AT_gamma_response(n_long,
      d_Gy_double,
      gamma_model_long,
      gamma_parameter_double,
      lethal_events_mode_bool,
      S_double);

  /* double -> float conversion (results) */
  for(i = 0 ; i < n_long ; i++){
    S[i] = (float)S_double[i];
  }

  free(d_Gy_double);
  free(S_double);
}

void AT_particle_name_from_particle_no_R(const int* particle_no,
    char** particle_name){

  const long n = 1;
  const long particle_no_long = (long)(*particle_no);
  char particle_name_str[1][6];

  strcpy(particle_name_str[0], *particle_name);

  AT_particle_name_from_particle_no( n,
      &particle_no_long,
      particle_name_str);

  strcpy(*particle_name, particle_name_str[0]);
}


void AT_particle_no_from_particle_name_R(const char** particle_name,
    int* particle_no){
  const long n = 1;
  char * particle_name_str = (char*)calloc(PARTICLE_NAME_NCHAR, sizeof(char));
  strcpy(particle_name_str, *particle_name);

  long particle_no_long = 0;

  AT_particle_no_from_particle_name( n,
      &particle_name_str,
      &particle_no_long);

  free( particle_name_str);

  *particle_no = (int)particle_no_long;
}

void AT_material_name_from_material_no_R( const int* material_no,
    char** material_name){

	const long material_no_long = (long)(*material_no);
	char * material_name_str = (char*)calloc(MATERIAL_NAME_LENGTH, sizeof(char));

	AT_material_name_from_number( material_no_long,
	      material_name_str);

	strcpy(*material_name, material_name_str);

	free(material_name_str);
}

void AT_material_no_from_material_name_R( const char** material_name,
	int* material_no){

	  char * material_name_str = (char*)calloc(MATERIAL_NAME_LENGTH, sizeof(char));
	  strcpy(material_name_str, *material_name);

	  long material_no_long = 0;

	  material_no_long	= AT_material_number_from_name( material_name_str);

	  free( material_name_str);

	  *material_no = (int)material_no_long;
}


