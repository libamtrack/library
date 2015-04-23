#include "AT_StoppingPowerData.h"

AT_stopping_power_functions_struct AT_stopping_power_functions =
		{
//				{ FromFile, Bethe, PSTAR, ICRU},
				{ "FromFile", "Bethe", "PSTAR", "ICRU"},    /**< source_name */
				{ &AT_FromFile_wrapper, &AT_Bethe_wrapper, &AT_PSTAR_wrapper, &AT_ICRU_wrapper}
		};


int AT_stopping_power_source_model_name_from_number( const long source_no, char* source_name){
//
//	assert( source_no > 0);
//	assert( source_no <= STOPPING_POWER_SOURCE_N);
//	assert(AT_stopping_power_sources.stopping_power_source_no[source_no-1] == source_no);
//	if( source_no <= 0)
//		return -1;
//	if( source_no > STOPPING_POWER_SOURCE_N)
//		return -1;
//    strcpy(source_name, AT_stopping_power_sources.stopping_power_source_name[source_no-1]);
    return AT_Success;
}


long AT_stopping_power_source_model_number_from_name( const char* source_name ){
//
//	long  match;
//	const long n_tmp = 1;
//
//	assert( source_name != NULL);
//
//	find_elements_char(  &source_name,
//			n_tmp,
//			AT_stopping_power_sources.stopping_power_source_name,
//			&match);
//
//	return (match+1);
	return(0);
}
