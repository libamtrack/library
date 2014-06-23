AC_PREREQ([2.68])
m4_define(gsl_required_version, 1.8)
AC_INIT([libamtrack], [0.5.5], [s.greilich@dkfz.de], [libamtrack])
AC_PACKAGE_URL([libamtrack.dkfz.org])
# Code names: 
# minor - 4: Ocelot, 5: Wombat
# micro - 0: Black, 1: Blue, 2: Red, 3: Green, 4: Yellow, 5: Violet
AC_DEFINE(CODE_NAME, ["Violet Wombat"], [Cryptonym to faciliate versioning for users])
# Status: "Development" or "Release" (will trigger some different actions, e.g. with R package)
AC_DEFINE(CODE_STATUS, ["Development"], [Status of the code])

LT_INIT
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/AT_Algorithms_CPP.c])
AC_CONFIG_SRCDIR([test/C/AT_test.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL


## HOW TO GET revision in here? But not trying to get it during a configure
## process on a tarball distribution??
##
## Solution 1 (use svn:keywords):
# Set revision automatically by svn:keywords
# AC_DEFINE(SVN_REVISION, $Rev$, [SVN Revision])
## Not optimal, as resulting string is not very handy and autoconf.ac has to be commited in
## order to update
##
## Solution 2 (use svnversion):
AC_PATH_PROG(svnversioncommand, svnversion, [AC_MSG_ERROR([unable to find svnversion])])
if test "X$svnversioncommand" = "X" || test `$svnversioncommand -n '.'` = "exported"; then
AC_DEFINE_UNQUOTED(SVN_VERSION, ["`cat saved_svn_version.txt`"], [repository svn version])
else
AC_DEFINE_UNQUOTED(SVN_VERSION, ["`svnversion -n`"], [repository svn version])
`svnversion -n > saved_svn_version.txt`
fi 
## Also not optimal as it will be used during the configure-process on a user machine, too.
## Additional critics (from Brian Tanner's Research Homepage):
## The SVN_VERSION will only get updated when you run configure.  
## If you are like me and "make distclean" before committing your latest revisions, 
## then this is fine, because next time you want to build your code you'll have 
## to ./configure and the values will be updated.  I guess the downside is that if you 
## test some changes, make distclean, then commit, it won't show that the version has 
## changed.  You'll have to ./configure and then make distclean one more time before committing.  
## I can live with that.  Someone more savvy might point out that there probably a hook you 
## could write somewhere to update saved_svn_version.txt at some point in the process.

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
if test "${GSL_LIB}" != ""; then
    	GSL_CFLAGS = "-I${GSL_LIB}/include"
	GSL_LIBS = "-L${GSL_LIB}/lib -lgsl -lgslcblas"
else
    AX_PATH_GSL(1.8, 
    		AC_DEFINE(HAVE_GSL_1_8, 1, [gsl>=1.8 found]), 
    		AC_MSG_ERROR([No gsl>=1.8 found!]))
fi

# Look for libm
AC_SEARCH_LIBS([cos], [m], [], [AC_MSG_ERROR([unable to find the cos() function])])


# Checks for essential header files
AC_CHECK_HEADER([malloc.h])
AC_CHECK_HEADER([stdint.h])
AC_CHECK_HEADER([stdlib.h])
AC_CHECK_HEADER([string.h])


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
