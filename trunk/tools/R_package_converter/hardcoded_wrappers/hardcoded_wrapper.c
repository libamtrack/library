#include "hardcoded_wrapper.h"

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

