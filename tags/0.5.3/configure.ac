AC_PREREQ([2.67])
AC_INIT([libamtrack], [0.5.3], [s.greilich@dkfz.de], [libamtrack])
AC_PACKAGE_URL([libamtrack.dkfz.org])
# Code names: 
# minor - 4: Ocelot, 5: Wombat
# micro - 0: Black, 1: Blue, 2: Red, 3: Green
AC_DEFINE(CODE_NAME, ["Green Wombat"], [Cryptonym to faciliate versioning for users])
# Status: "Development" or "Release" (will trigger some different actions, e.g. with R package)
AC_DEFINE(CODE_STATUS, ["Release"], [Status of the code])



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
# AC_DEFINE(SVN_REVISION, $Rev: 1187 $, [SVN Revision])
## Not optimal, as resulting string is not very handy and autoconf.ac has to be commited in
## order to update
##
## Solution 2 (use svnversion):
AC_PATH_PROG(svnversioncommand, svnversion)
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
if test "${LIB_GSL}" != ""; then
    GSL_CFLAGS="-I${LIB_GSL}/include"
	GSL_LIBS="-L${LIB_GSL}/lib -lgsl"
else
        AX_PATH_GSL(1.8, AC_DEFINE(HAVE_GSL_1_8, 1, [gsl>=1.8 found]), AC_MSG_ERROR([No gsl>=1.8 found!]))
fi

AC_SEARCH_LIBS(cos, m)

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
