 /* example.i */
 %module example
 %include "typemaps.i"

 %{
 /* Put header files here or function declarations like below */
 #include "example.h"
 #include "AmTrack.h" 
 #include "AT_ElectronRange.h"
 %}
 
 
 extern int fact(int n);
 extern int my_mod(int x, int y);
 extern char *get_time();
 extern void AT_PrintName(  void);
 extern int AT_GetNumber(  void);

%include "arrays_java.i"
%apply float[] {float *};
    
extern void AT_max_electron_range_m( const long  n,
    const float*  E_MeV_u,
    const long  material_no,
    const long   er_model,
    float*  max_electron_range_m);

