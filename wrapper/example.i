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
 
 /* If you need to interface new function, add it here */
 

/************************ R ******************************/
/* SWIG does not support (yet) passing of arrays between */
/* R and C. Only simple data types are supported         */

/********************** JAVA *****************************/
/* Following trick will enable passing arrays as pointers */
/* Arrays of long (long* pointers) are not supported on 64-bit JVM */
/* Arrays of floats, doubles are fully supported */
/* Single value passed by value is supported for any type */ 

%include "arrays_java.i"
%apply double[] {double *};
    
extern void AT_max_electron_ranges_m( const long number_of_particles,
    const double E_MeV_u[],
    const long   material_no,
    const long   er_model,
    double  max_electron_range_m[]);

