#ifndef AT_GAMMARESPONSE_H_
#define AT_GAMMARESPONSE_H_

/**
 * @brief Gamma response models
 */

/*
 *    AT_GammaResponse.h
 *    ==================
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "gsl/gsl_pow_int.h"
#include "gsl/gsl_sf_gamma.h"

#include "AT_Constants.h"
#include "AT_NumericalRoutines.h"


/**
 * Gamma Response code numbers
 */
enum AT_GammaResponseModels{
  GR_Test                  = 1,      /**< no parameters */
      GR_GeneralTarget     = 2,      /**< TODO */
      GR_Radioluminescence = 3,      /**< 0 - Smax, 1 - D0, 2 - dyn */
      GR_ExpSaturation     = 4,      /**< 0 - Smax, 1 - D0 */
      GR_LinQuad           = 5,      /**< 0 - alpha, 1 - beta, 2 - D0 */
      GR_Geiss			   = 6,      /**< TODO */
	  GR_AM_DSBEnhancement = 7       /**< DSB enhancement factor (A. Mairani) */
};


#define GR_DATA_N    7

#define GR_MAX_NUMBER_OF_PARAMETERS  9
//TODO Find better solution for array size of parameters

/**
 * @struct AT_GR_data_struct
 */
typedef struct {
  const long     n;
  const long     GR_no[GR_DATA_N];
  const long     n_parameters[GR_DATA_N];
  const char*    parameter_name[GR_DATA_N][GR_MAX_NUMBER_OF_PARAMETERS];
  const double   parameter_default[GR_DATA_N][GR_MAX_NUMBER_OF_PARAMETERS];
  const char*    GR_name[GR_DATA_N];
} AT_GR_data_struct;


/**
 * TODO
 */
static const AT_GR_data_struct AT_GR_Data = {
    GR_DATA_N,
    {  GR_Test,
       GR_GeneralTarget,
       GR_Radioluminescence,
       GR_ExpSaturation,
       GR_LinQuad,
       GR_Geiss,
       GR_AM_DSBEnhancement
    },
    {  2,
	   5,
	   3,
	   2,
	   3,
	   5,
	   0},
    {  {"",      "",      "",      "",   ""},
       {"S_max", "D0_Gy", "c",     "m",  ""},
       {"S_max", "D0_Gy", "dyn",   "",   ""},
       {"S_max", "D0_Gy", "",      "",   ""},
       {"alpha", "beta",  "D0_Gy", "",   ""},
       {"const", "k1",    "a1",    "k2", "a2"},
	   {"",      "",      "",      "",   ""}
    },
    {  {0.,  0.,  0.,   0.,  0.},
       {1.,  10., 1.,   1.,  0.},
       {1.,  10., 5.,   0.,  0.},
       {1.,  10., 0.,   0.,  0.},
       {0.2, 0.02, 10.,  0.,  0.},
       {6e3, 0.8, 3e-4, 0.2, 1e-6},
       {0.,  0.,  0.,   0.,  0.}
	  },
    {  "simple test gamma response",
       "generalized multi-target/multi-hit gamma response",
       "radioluminescence gamma response",
       "exp.-sat. gamma response (obsolete, use gen. target/hit instead)",
       "linear-quadratic gamma response",
       "Geiss model (1997)",
       "DSB enhancement factor (A. Mairani)"
    }
};


/**
 * Get index of gamma response model in AT_GR_Data for given Gamma_no
 * (currently for example gamma response model with number 2 has index 1)
 *
 * @param[in] Gamma_no  gamma response model number
 * @return          gamma response model index in AT_GR_Data table
 */
long AT_Gamma_index_from_material_number( const long Gamma_no );


/**
 * Returns name of the gamma response model from model number
 *
 * @param[in]  Gamma_no    gamma response model number
 * @param[out] Gamma_name  string containing gamma response model name
 */
void AT_Gamma_name_from_number(  const long Gamma_no,
    char* Gamma_name);


/**
 * Returns number of parameters of the gamma response model from model number
 *
 * @param[in]   Gamma_no   gamma response model number
 * return                  number of GR parameters
 */
long AT_Gamma_number_of_parameters( const long Gamma_no );


/**
 * Returns a system (detector or cells) response for given doses
 * according to the chosen gamma response model
 *
 * @param[in]  number_of_doses  number of doses given in vector d_Gy
 * @param[in]  d_Gy             doses in Gy (array of size number_of_doses)
 * @param[in]  gamma_model      gamma response model index
 * @param[in]  gamma_parameter  vector holding necessary parameters for the chose gamma response model (array of size 9)
 * @param[in]  lethal_event_mode  if true computation is done in lethal event mode
 * @param[out] response         gamma responses (array of size number_of_doses)
 */
void AT_gamma_response( const long  number_of_doses,
    const double   d_Gy[],
    const long     gamma_model,
    const double   gamma_parameter[],
    const bool	   lethal_event_mode,
    double         S[]);


/**
 * Returns the detector / cell gamma response for dose distribution
 *
 * @param[in]  number_of_bins            number of bins in the dose histogram
 * @param[in]  dose_Gy_bin_position      midpoint doses for histogram in Gy (array of size number_of_bins)
 * @param[in]  dose_Gy_bin_width         bin widths for histogram in Gy (array of size number_of_bins)
 * @param[in]  dose_bin_frequency        dose frequencies for histogram (array of size number_of_bins)
 * @param[in]  gamma_model               gamma response model index
 * @param[in]  gamma_parameter           vector holding necessary parameters for the chose gamma response model (array of size GR_MAX_NUMBER_OF_PARAMETERS)
 * @param[in]  lethal_events_mode        if true computation is done in lethal event mode
 * @return     response
 */
double AT_get_gamma_response_for_average_dose(  const long  number_of_bins,
    const double   dose_Gy_bin_position[],
    const double   dose_Gy_bin_width[],
    const double   dose_bin_frequency[],
    const long     gamma_model,
    const double   gamma_parameter[],
    const bool     lethal_events_mode);


/**
 * Returns the detector / cell gamma response for dose distribution
 *
 * @param[in]  number_of_bins            number of bins in the dose histogram
 * @param[in]  dose_Gy_bin_position      midpoint doses for histogram in Gy (array of size number_of_bins)
 * @param[in]  dose_bin_frequency        dose frequencies for histogram (array of size number_of_bins)
 * @param[in]  gamma_model               gamma response model index
 * @param[in]  gamma_parameter           vector holding necessary parameters for the chose gamma response model (array of size GR_MAX_NUMBER_OF_PARAMETERS)
 * @param[in]  lethal_events_mode        if true computation is done in lethal event mode
 * @param[out] response_bin_frequency    resulting response frequencies (array of size number_of_bins)
 */
void AT_get_response_distribution_from_dose_distribution(  const long  number_of_bins,
    const double   dose_Gy_bin_position[],
    const double   dose_bin_frequency[],
    const long     gamma_model,
    const double   gamma_parameter[],
    const bool     lethal_events_mode,
    double         response_bin_frequency[] );


/**
 * Returns the ion response from an ion response distribution
 *
 * @param[in]  number_of_bins                number of bins in the dose histogram
 * @param[in]  dose_Gy_bin_width             bin widths for histogram in Gy (array of size number_of_bins)
 * @param[in]  dose_bin_frequency            dose frequencies for histogram (array of size number_of_bins)
 * @param[in]  ion_response_bin_frequency    ion response frequencies (array of size number_of_bins)
 * @return     response
 */
double AT_get_ion_response_from_response_distribution(  const long  number_of_bins,
	const double   dose_Gy_bin_width[],
    const double   dose_bin_frequency[],
    const double   ion_response_bin_frequency[]);


/**
 * Returns ion response for dose distribution
 *
 * @param[in]  number_of_bins            number of bins in the dose histogram
 * @param[in]  dose_Gy_bin_position      midpoint doses for histogram in Gy (array of size number_of_bins)
 * @param[in]  dose_bin_frequency        dose frequencies for histogram (array of size number_of_bins)
 * @param[in]  gamma_model               gamma response model index
 * @param[in]  gamma_parameter           vector holding necessary parameters for the chose gamma response model (array of size GR_MAX_NUMBER_OF_PARAMETERS)
 * @param[in]  lethal_events_mode        if true computation is done in lethal event mode
 * @return     resulting ion response
 */
double AT_get_ion_response_from_dose_distribution(  const long  number_of_bins,
	    const double   dose_Gy_bin_position[],
	    const double   dose_Gy_bin_width[],
	    const double   dose_bin_frequency[],
	    const long     gamma_model,
	    const double   gamma_parameter[],
	    const bool     lethal_events_mode);


/**
 * Returns relative efficiency for dose distribution
 *
 * @param[in]  number_of_bins            number of bins in the dose histogram
 * @param[in]  dose_Gy_bin_position      midpoint doses for histogram in Gy (array of size number_of_bins)
 * @param[in]  dose_bin_frequency        dose frequencies for histogram (array of size number_of_bins)
 * @param[in]  gamma_model               gamma response model index
 * @param[in]  gamma_parameter           vector holding necessary parameters for the chose gamma response model (array of size GR_MAX_NUMBER_OF_PARAMETERS)
 * @param[in]  lethal_events_mode        if true computation is done in lethal event mode
 * @return     relative efficiency
 */
double AT_get_ion_efficiency_from_dose_distribution(  const long  number_of_bins,
	    const double   dose_Gy_bin_position[],
	    const double   dose_Gy_bin_width[],
	    const double   dose_bin_frequency[],
	    const long     gamma_model,
	    const double   gamma_parameter[],
	    const bool     lethal_events_mode);


/**
 * Returns relative efficiency from ion response distribution
 *
 * @param[in]  number_of_bins            number of bins in the dose histogram
 * @param[in]  dose_Gy_bin_position      midpoint doses for histogram in Gy (array of size number_of_bins)
 * @param[in]  dose_bin_frequency        dose frequencies for histogram (array of size number_of_bins)
 * @param[in]  ion_response_bin_frequency    ion response frequencies (array of size number_of_bins)
 * @param[in]  gamma_model               gamma response model index
 * @param[in]  gamma_parameter           vector holding necessary parameters for the chose gamma response model (array of size GR_MAX_NUMBER_OF_PARAMETERS)
 * @param[in]  lethal_events_mode        if true computation is done in lethal event mode
 */
double AT_get_ion_efficiency_from_response_distribution(  const long  number_of_bins,
	    const double   dose_Gy_bin_position[],
	    const double   dose_Gy_bin_width[],
	    const double   dose_bin_frequency[],
	    const double   ion_response_bin_frequency[],
	    const long     gamma_model,
	    const double   gamma_parameter[],
	    const bool     lethal_events_mode);


/**
 * Returns the detector / cell gamma response for a local dose distribution
 * according to the chosen gamma response model, used by response model
 * routines in AmTrack.c
 *
 * @param[in]  number_of_bins      number of bins in given local dose distribution
 * @param[in]  d_Gy                local dose bin position in Gy (array of size number_of_bins)
 * @param[in]  dd_Gy               local dose bin width in Gy (array of size number_of_bins)
 * @param[in]  f                   local dose frequency (array of size number_of_bins)
 * @param[in]  f0                  frequency of zero local dose
 * @param[in]  gamma_model         gamma response model index
 * @param[in]  gamma_parameter     vector holding necessary parameters for the chose gamma response model (array of size GR_MAX_NUMBER_OF_PARAMETERS)
 * @param[in]  lethal_events_mode  if true, allows to do calculations for cell survival
 *    @see  AmTrack.c/AT_IGK
 * @param[out] S                   gamma responses for given bins (array of size number_of_bins)
 * @param[out] S_HCP               HCP response for given local dose distribution (expectation value of S distribution)
 * @param[out] S_gamma             gamma response for given local dose distribution (gamma response of expectation value of d distribution)
 * @param[out] efficiency          RE = S_HCP/S_gamma for given local dose distribution
 */
void AT_get_gamma_response(  const long  number_of_bins,
    const double   d_Gy[],
    const double   dd_Gy[],
    const double   f[],
    const double   f0, // TODO parameter not used here, might be removed
    const long     gamma_model,
    const double   gamma_parameter[],
    const bool     lethal_events_mode,
    double         S[],
    double*        S_HCP,
    double*        S_gamma,
    double*        efficiency);


#endif // AT_GAMMARESPONSE_H_
