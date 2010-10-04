#****************************************************************************
# Makefile for libamtrack
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
#  MINGWPATH = C:\Program Files\MinGW (Windows)
MINGWPATH = 
############################################################

CFLAGS    = -Wall -c -O3 -fPIC
LFLAGS    = -lm -lgsl -lgslcblas

ifeq ($(OS),Linux)
NAMELIB   = lib/libamtrack.so
DSEP      = /
SRCDIR    = ./src
INCLDIR   = ./include
RMCMD     = rm
GCC       = gcc 
GSLINCLUDE = "$(GSLPATH)/include"
GSLLIB     = "$(GSLPATH)/lib"
else
NAMELIB    = lib\libamtrack.dll
DSEP       = \\
SRCDIR     = .\src
INCLDIR    = .\include
RMCMD      = del
GCC        = "$(MINGWPATH)\bin\gcc.exe"
GSLINCLUDE = "$(GSLPATH)\include"
GSLLIB     = "$(GSLPATH)\lib"
endif

LIBCOBJS  = $(SRCDIR)$(DSEP)AmTrack.c $(SRCDIR)$(DSEP)AT_Error.c $(SRCDIR)$(DSEP)AT_Constants.c $(SRCDIR)$(DSEP)AT_GammaResponse.c $(SRCDIR)$(DSEP)AT_DataLET.c $(SRCDIR)$(DSEP)AT_DataMaterial.c $(SRCDIR)$(DSEP)AT_DataParticle.c $(SRCDIR)$(DSEP)AT_ElectronRange.c $(SRCDIR)$(DSEP)AT_GammaResponse.c $(SRCDIR)$(DSEP)AT_NumericalRoutines.c $(SRCDIR)$(DSEP)AT_PhysicsRoutines.c $(SRCDIR)$(DSEP)AT_RDD.c $(SRCDIR)$(DSEP)AT_RDD_Simple.c $(SRCDIR)$(DSEP)AT_RDD_ShellAveraged.c $(SRCDIR)$(DSEP)AT_RDD_ExtendedTarget.c $(SRCDIR)$(DSEP)AT_KatzModel.c $(SRCDIR)$(DSEP)AT_SuccessiveConvolutions.c $(SRCDIR)$(DSEP)AT_Wrapper_R.c $(SRCDIR)$(DSEP)AT_Histograms.c $(SRCDIR)$(DSEP)AT_Algorithms_GSM.c

LIBHOBJS  = $(INCLDIR)$(DSEP)AmTrack.h $(INCLDIR)$(DSEP)AT_Error.h $(INCLDIR)$(DSEP)AT_Constants.h $(INCLDIR)$(DSEP)AT_GammaResponse.h $(INCLDIR)$(DSEP)AT_DataLET.h $(INCLDIR)$(DSEP)AT_DataMaterial.h $(INCLDIR)$(DSEP)AT_DataParticle.h $(INCLDIR)$(DSEP)AT_ElectronRange.h $(INCLDIR)$(DSEP)AT_GammaResponse.h $(INCLDIR)$(DSEP)AT_NumericalRoutines.h $(INCLDIR)$(DSEP)AT_PhysicsRoutines.h $(INCLDIR)$(DSEP)AT_RDD.h $(INCLDIR)$(DSEP)AT_RDD_Simple.h $(INCLDIR)$(DSEP)AT_RDD_ShellAveraged.h $(INCLDIR)$(DSEP)AT_RDD_ExtendedTarget.h $(INCLDIR)$(DSEP)AT_KatzModel.h $(INCLDIR)$(DSEP)AT_SuccessiveConvolutions.h $(INCLDIR)$(DSEP)AT_Wrapper_R.h $(INCLDIR)$(DSEP)AT_Histograms.h $(INCLDIR)$(DSEP)AT_Algorithms_GSM.h

LIBOBJS  = AmTrack.o AT_Error.o AT_Constants.o AT_DataLET.o AT_DataMaterial.o AT_DataParticle.o AT_ElectronRange.o AT_GammaResponse.o AT_NumericalRoutines.o AT_PhysicsRoutines.o AT_RDD.o AT_RDD_Simple.o AT_RDD_ShellAveraged.o AT_RDD_ExtendedTarget.o AT_SuccessiveConvolutions.o AT_KatzModel.o AT_Wrapper_R.o AT_Histograms.o AT_Algorithms_GSM.o

all:$(LIBOBJS)
		echo $(GSLPATH)$(DSEP)
		$(GCC) -L$(GSLLIB) -shared $(LIBOBJS) -o $(NAMELIB) $(LFLAGS) 
		$(RMCMD) *.o

clean:
		- $(RMCMD) *.o
		- $(RMCMD) $(NAMELIB)

$(LIBOBJS):$(LIBCOBJS) $(LIBHOBJS)
		$(GCC) -I$(INCLDIR) -I$(GSLINCLUDE) $(CFLAGS) $(LIBCOBJS)
