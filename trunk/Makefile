#****************************************************************************
#
# Makefile for libamtrack, sgreilich, 2009/09/23
#
#****************************************************************************

############################################################
# USER SETTINGS
############################################################
## Please specify your OS here
#  OS        = Linux (or)
#  OS        = Windows
OS        = Linux
## Please set the correct path to your GSL installation here:
#  most likely:
#  GSLPATH   = /usr (Linux)
#  GSLPATH   = C:\Program Files\GnuWin32 (Windows)
GSLPATH   = /usr
## Under Windows please specify the path to your MinGW installation
#  leave empty under Linux: MINGWPATH = 
#  for Windows most likely:
#  MINGWPATH = C:\Programme\MinGW (Windows)
MINGWPATH = 
############################################################

CFLAGS    = -Wall -c -O3 -fPIC
LFLAGS    = -lm -lgsl -lgslcblas

ifeq ($(OS),Linux)
NAMELIB   = libamtrack.so
NAMEEXE   = AmTrack
DSEP      = /
SRCDIR    = ./src
INCLDIR   = ./include
RMCMD     = rm
GCCDIR    = 
else
NAMELIB   = libamtrack.dll
NAMEEXE   = AmTrack.exe
DSEP      = \\
SRCDIR    = .\src
INCLDIR   = .\include
RMCMD     = del
GCCDIR    = $(MINGWPATH)\bin$(DSEP)
endif

LIBCOBJS  = $(SRCDIR)$(DSEP)AmTrack.c $(SRCDIR)$(DSEP)AT_Constants.c $(SRCDIR)$(DSEP)AT_GammaResponse.c $(SRCDIR)$(DSEP)AT_DataLET.c $(SRCDIR)$(DSEP)AT_DataMaterial.c $(SRCDIR)$(DSEP)AT_DataParticle.c $(SRCDIR)$(DSEP)AT_ElectronRange.c $(SRCDIR)$(DSEP)AT_GammaResponse.c $(SRCDIR)$(DSEP)AT_NumericalRoutines.c $(SRCDIR)$(DSEP)AT_PhysicsRoutines.c $(SRCDIR)$(DSEP)AT_RDD.c $(SRCDIR)$(DSEP)AT_RDD_Simple.c $(SRCDIR)$(DSEP)AT_RDD_ShellAveraged.c $(SRCDIR)$(DSEP)AT_RDD_ExtendedTarget.c $(SRCDIR)$(DSEP)AT_KatzModel.c $(SRCDIR)$(DSEP)AT_SuccessiveConvolutions.c $(SRCDIR)$(DSEP)AT_Wrapper_R.c

LIBHOBJS  = $(INCLDIR)$(DSEP)AmTrack.h $(INCLDIR)$(DSEP)AT_Constants.h $(INCLDIR)$(DSEP)AT_GammaResponse.h $(INCLDIR)$(DSEP)AT_DataLET.h $(INCLDIR)$(DSEP)AT_DataMaterial.h $(INCLDIR)$(DSEP)AT_DataParticle.h $(INCLDIR)$(DSEP)AT_ElectronRange.h $(INCLDIR)$(DSEP)AT_GammaResponse.h $(INCLDIR)$(DSEP)AT_NumericalRoutines.h $(INCLDIR)$(DSEP)AT_PhysicsRoutines.h $(INCLDIR)$(DSEP)AT_RDD.h $(INCLDIR)$(DSEP)AT_RDD_Simple.h $(INCLDIR)$(DSEP)AT_RDD_ShellAveraged.h $(INCLDIR)$(DSEP)AT_RDD_ExtendedTarget.h $(INCLDIR)$(DSEP)AT_KatzModel.h $(INCLDIR)$(DSEP)AT_SuccessiveConvolutions.h $(INCLDIR)$(DSEP)AT_Wrapper_R.h

LIBOBJS  = AmTrack.o AT_Constants.o AT_DataLET.o AT_DataMaterial.o AT_DataParticle.o AT_ElectronRange.o AT_GammaResponse.o AT_NumericalRoutines.o AT_PhysicsRoutines.o AT_RDD.o AT_RDD_Simple.o AT_RDD_ShellAveraged.o AT_RDD_ExtendedTarget.o AT_SuccessiveConvolutions.o AT_KatzModel.o AT_Wrapper_R.o

all:$(LIBOBJS)
		echo $(GSLPATH)$(DSEP)
		$(GCCDIR)gcc -L$(GSLPATH)$(DSEP)lib -shared $(LIBOBJS) -o $(NAMELIB) $(LFLAGS) 
		$(RMCMD) *.o

UI:AT_UI.o $(LIBOBJS)
		$(GCCDIR)gcc -L$(GSLPATH)$(DSEP)lib $(LIBOBJS) AT_UI.o -o $(NAMEEXE) $(LFLAGS) 
		$(RMCMD) *.o

test:AT_test.o $(LIBOBJS)
		$(GCCDIR)gcc -L$(GSLPATH)$(DSEP)lib $(LIBOBJS) AT_test.o -o AT_test.exe $(LFLAGS) 
		$(RMCMD) *.o

clean:
		- $(RMCMD) *.o
		- $(RMCMD) $(NAMELIB)
		- $(RMCMD) $(NAMEEXE)
		- $(RMCMD) AT_test.exe

$(LIBOBJS):$(LIBCOBJS) $(LIBHOBJS)
		$(GCCDIR)gcc -I$(INCLDIR) -I$(GSLPATH)$(DSEP)include $(CFLAGS) $(LIBCOBJS)

AT_UI.o:test$(DSEP)C$(DSEP)AT_UI.c 
		$(GCCDIR)gcc -I$(INCLDIR) -I$(GSLPATH)$(DSEP)include $(CFLAGS) test$(DSEP)C$(DSEP)AT_UI.c

AT_test.o:test$(DSEP)C$(DSEP)AT_test.c 
		$(GCCDIR)gcc -I$(INCLDIR) -I$(GSLPATH)$(DSEP)include $(CFLAGS) test$(DSEP)C$(DSEP)AT_test.c

