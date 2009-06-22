#include <stdio.h>
#include <stdlib.h>

#define _WINDOWS // [_LINUX or _WINDOWS] : in Linux we have isnan function while in Windows we have _isnan
#define _S // [_S or _R] in S we can pass long type to the function via as.single, but in R we pass int type
//#define _DEBUG // debugging printouts
#define _SOLVER	// use SOLVER instead of analytical inversion in r_RDD_m

#include "SGP_Constants.h"
#include "SGP_Data.h"
#include "SGP_Functions.h"
#include "SGP_RDD.h"
#include "SGP_SuccessiveConvolutions.h"
#include "SGP_GammaResponse.h"
#include "SGP_FileOperations.h"
#include "SGP_Transport.h"

/*
BOOL APIENTRY DllMain( HANDLE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
    return TRUE;
}
*/

