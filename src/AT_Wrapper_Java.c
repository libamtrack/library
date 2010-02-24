#ifdef _JAVA

#include "AT_Wrapper_Java.h"

/**
 * Wrapper for Java Interface Test
 * TODO: Move to AT_Wrapper_Java.c
 */
JNIEXPORT jfloatArray JNICALL Java_AmTrack_AmTrackTestGUI_ATconvertBeamParameters
  (JNIEnv *pEnv, jclass clazz, jfloat fluenceCm2, jfloat sigmaCm, jfloat n, jfloat fwhmMm)
{
  long  n_data          =       1;
  float fluence_cm2     =       (float)fluenceCm2;
  float sigma_cm        =       (float)sigmaCm;
  float N               =       (float)n;
  float FWHM_mm         =       (float)fwhmMm;

  AT_convert_beam_parameters(     &n_data,
                                  &fluence_cm2,
                                  &sigma_cm,
                                  &N,
                                  &FWHM_mm);

  jfloatArray results;
  results              =       (*pEnv)->NewFloatArray(pEnv, 4);
  jfloat* presults     =       (*pEnv)->GetFloatArrayElements(pEnv, results, 0);
  presults[0]          =       (jfloat)fluence_cm2;
  presults[1]          =       (jfloat)sigma_cm;
  presults[2]          =       (jfloat)N;
  presults[3]          =       (jfloat)FWHM_mm;
  (*pEnv)->SetFloatArrayRegion(pEnv, results, 0, 4, presults);
  return results;
}

#endif
