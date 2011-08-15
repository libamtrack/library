AC_PREREQ([2.67])
AC_INIT([libamtrack], [0.5.2], [s.greilich@dkfz.de])
LT_INIT
#AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/AT_Algorithms_CPP.c])
AC_CONFIG_SRCDIR([test/C/AT_test.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_MAKE_SET

# Checks for libraries, add uset-set path to gsl if necessary (e.g. on Windows)
if test "${LIB_GSL}" != ""; then
	CFLAGS="${CFLAGS} -I${LIB_GSL}/include"
	LDFLAGS="${LDFLAGS} -L${LIB_GSL}/lib -lgsl"
fi
#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_SEARCH_LIBS(cos, m)
AC_SEARCH_LIBS([cblas_dgemm], [gslcblas], [], [
  AC_MSG_ERROR([unable to find the cblas_dgemm() function in gslcblas])
])
AC_SEARCH_LIBS([gsl_blas_dgemm], [gsl], [], [
  AC_MSG_ERROR([unable to find the gsl_blas_dgemm() function in gsl])
])

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

AC_CONFIG_FILES(Makefile src/Makefile test/Makefile test/C/Makefile example/Makefile example/basic_plots/Makefile example/demo/Makefile)
AC_OUTPUT
