AC_INIT([src/AT_Algorithms_CPP.c])

# Following lines kindly supplied by Dirk Eddelbuettel
## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
   # Use gsl-config for header and linker arguments
   GSL_CFLAGS=`${GSL_CONFIG} --cflags`
   GSL_LIBS=`${GSL_CONFIG} --libs`
else
   AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

AC_MSG_CHECKING([if GSL version >= 1.8])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <string.h>
#include <gsl/gsl_version.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit (minor < 8);
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 1.8])
else
  AC_MSG_RESULT([yes])
  AC_SUBST([GSL_CFLAGS])
  AC_SUBST([GSL_LIBS])
fi

# Check if libmath (CERNLIB) and appropriate header are available
AC_CHECK_HEADER([cfortran.h], [], AC_DEFINE(HAVE_NO_CERNLIB, 1, [No cfortran.h found.]))
AC_CHECK_HEADER([gen.h], [], AC_DEFINE(HAVE_NO_CERNLIB, 1, [No cfortran.h found.]))
AC_CHECK_LIB(mathlib, vavden_, [], AC_DEFINE(HAVE_NO_CERNLIB, 1, [No cfortran.h found.]))
	
if test "x${HAVE_NO_CERNLIB}" != ""; then
	AC_SUBST([CERNLIB_CFLAGS], ["-DHAVE_CERNLIB"])
	AC_SUBST([CERNLIB_LIBS], ["-lmathlib"])
else
	AC_MSG_RESULT([No mathlib CERNLIB or headers found, some energy loss distribution functions will not be available.])
	AC_SUBST([CERNLIB_CFLAGS])
	AC_SUBST([CERNLIB_LIBS])
fi

AC_OUTPUT(src/Makevars)

