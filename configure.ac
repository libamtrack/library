AC_PREREQ([2.68])
m4_define(gsl_required_version, 1.8)
AC_INIT([libamtrack], [0.6.3], [s.greilich@dkfz.de], [libamtrack])
AC_PACKAGE_URL([libamtrack.dkfz.org])
# Code names: 
# minor - 4: Ocelot, 5: Wombat, 6: Armadillo
# micro - 0: Black, 1: Blue, 2: Red, 3: Green, 4: Yellow, 5: Violet
AC_DEFINE(CODE_NAME, ["Green Armadillo"], [Cryptonym to faciliate versioning for users])
# Status: "Development" or "Release" (will trigger some different actions, e.g. with R package)
AC_DEFINE(CODE_STATUS, ["Release"], [Status of the code])

LT_INIT
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/AT_Algorithms_CPP.c])
AC_CONFIG_SRCDIR([test/C/AT_test.c])
AC_CONFIG_SRCDIR([thirdparty/cernlib/AT_CernlibFuns.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_PROG_F77([gfortran g77])


# check for OpenMP flags
AC_ARG_ENABLE(openmp,
          [AS_HELP_STRING([--enable-openmp],[use OpenMP directives for parallelism])],
          enable_openmp=$enableval, enable_openmp=no)
          
if test "$enable_openmp" = "yes"; then
      AC_DEFINE(HAVE_OPENMP,1,[Define to 1 to enable OpenMP.])
      OPENMP_CFLAGS="-fopenmp"
      OPENMP_LDFLAGS="-lgomp"
fi

AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_LDFLAGS)

# Checks for libraries, add user-set path to gsl if necessary (e.g. on Windows)
AC_CHECK_LIB([m],[cos])

if test "${GSL_LIB}" != ""; then
    	GSL_CFLAGS = "-I${GSL_LIB}/include"
	GSL_LIBS = "-L${GSL_LIB}/lib -lgsl -lgslcblas"
else
	AC_CHECK_LIB([gslcblas],[cblas_dgemm])
	AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
fi

# Checks for essential header files
AC_CHECK_HEADER([malloc.h])
AC_CHECK_HEADER([stdint.h])
AC_CHECK_HEADER([stdlib.h])
AC_CHECK_HEADER([string.h])


# Checks for header files.
AC_CHECK_HEADERS([malloc.h stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow sqrt strstr strtol])

# Check for system and set flag for not doing R-package test in mingw
AC_CANONICAL_HOST
RTEST_FLAG=""
AS_CASE([$host_os],
        [mingw*],
        [
           RTEST_FLAG="--no-test"
		   AC_MSG_NOTICE([Running on MinGW, will suppress check of R package])
        ])
AC_SUBST([RTEST_FLAG])

AC_CONFIG_FILES(Makefile src/Makefile test/Makefile test/C/Makefile example/Makefile example/basic_plots/Makefile example/demo/Makefile thirdparty/cernlib/Makefile)
AC_OUTPUT
