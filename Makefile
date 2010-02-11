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


LIBCOBJS  = $(SRCDIR)$(DSEP)AmTrack.c $(SRCDIR)$(DSEP)AT_Constants.c $(SRCDIR)$(DSEP)AT_GammaResponse.c $(SRCDIR)$(DSEP)AT_DataLET.c $(SRCDIR)$(DSEP)AT_DataMaterial.c $(SRCDIR)$(DSEP)AT_DataParticle.c $(SRCDIR)$(DSEP)AT_ElectronRange.c $(SRCDIR)$(DSEP)AT_GammaResponse.c $(SRCDIR)$(DSEP)AT_NumericalRoutines.c $(SRCDIR)$(DSEP)AT_PhysicsRoutines.c $(SRCDIR)$(DSEP)AT_RDD.c $(SRCDIR)$(DSEP)AT_SuccessiveConvolutions.c $(SRCDIR)$(DSEP)AT_UI.c

LIBHOBJS  = $(INCLDIR)$(DSEP)AmTrack.h $(INCLDIR)$(DSEP)AT_Constants.h $(INCLDIR)$(DSEP)AT_GammaResponse.h $(INCLDIR)$(DSEP)AT_DataLET.h $(INCLDIR)$(DSEP)AT_DataMaterial.h $(INCLDIR)$(DSEP)AT_DataParticle.h $(INCLDIR)$(DSEP)AT_ElectronRange.h $(INCLDIR)$(DSEP)AT_GammaResponse.h $(INCLDIR)$(DSEP)AT_NumericalRoutines.h $(INCLDIR)$(DSEP)AT_PhysicsRoutines.h $(INCLDIR)$(DSEP)AT_RDD.h $(INCLDIR)$(DSEP)AT_SuccessiveConvolutions.h 

LIBOBJS  = AmTrack.o AT_Constants.o AT_DataLET.o AT_DataMaterial.o AT_DataParticle.o AT_ElectronRange.o AT_GammaResponse.o AT_NumericalRoutines.o AT_PhysicsRoutines.o AT_RDD.o AT_SuccessiveConvolutions.o


all:$(LIBOBJS)
		echo $(GSLPATH)$(DSEP)
		$(GCCDIR)gcc -L$(GSLPATH)$(DSEP)lib -shared $(LIBOBJS) -o $(NAMELIB) $(LFLAGS) 
		$(RMCMD) *.o

UI:AT_UI.o $(LIBOBJS)
		$(GCCDIR)gcc -L$(GSLPATH)$(DSEP)lib $(LIBOBJS) AT_UI.o -o $(NAMEEXE) $(LFLAGS) 
		$(RMCMD) *.o

clean:
		$(RMCMD) *.o
		$(RMCMD) $(NAMELIB)
		$(RMCMD) $(NAMEEXE)

$(LIBOBJS):$(LIBCOBJS) $(LIBHOBJS)
		$(GCCDIR)gcc -I$(INCLDIR) -I$(GSLPATH)$(DSEP)include $(CFLAGS) $(LIBCOBJS)


AT_UI.o:$(SRCDIR)$(DSEP)AT_UI.c 
		$(GCCDIR)gcc -I$(INCLDIR) -I$(GSLPATH)$(DSEP)include $(CFLAGS) $(SRCDIR)$(DSEP)AT_UI.c

# BELOW: sniplets from Leszeks makefile 2009/12/27
# TODO: merge makefiles
#
#LIBOBJFILES  = AT_Constants.o AT_FileOperations.o AT_Functions.o AT_GammaResponse.o AT_ParabolicCylinderFunction.o AT_RDD.o AT_SuccessiveConvolutions.o AT_Transport.o AT_Utils.o AmTrack.o 
#
#OBJFILES := $(patsubst %.c,%.o,$(wildcard *.c))
#
#all: libAT.so
#		echo "Linking"
#
#clean:
#		rm -f *.o
#
#libAT.so: $(LIBOBJFILES)
#		gcc $(LDFLAGS) -I../include/ -o libAT.so $(LIBOBJFILES)
#
#$(OBJFILES): *.c
#		gcc $(CFLAGS) -I../include/ -c -o $@ $<


