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
OS        = Windows
## Please set the correct path to your GSL installation here:
#  most likely:
#  GSLPATH   = /usr (Linux)
#  GSLPATH   = C:\Program Files\GnuWin32 (Windows)
GSLPATH   = C:\Programme\GnuWin32
## Under Windows please specify the path to your MinGW installation
#  leave empty under Linux: MINGWPATH = 
#  for Windows most likely:
#  MINGWPATH = C:\Programme\MinGW (Windows)
MINGWPATH = C:\Programme\MinGW
############################################################

CFLAGS    = -Wall -c -O3
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


LIBCOBJS  = $(SRCDIR)$(DSEP)AmTrack.c $(SRCDIR)$(DSEP)AT_Constants.c $(SRCDIR)$(DSEP)AT_FileOperations.c $(SRCDIR)$(DSEP)AT_Functions.c $(SRCDIR)$(DSEP)AT_GammaResponse.c $(SRCDIR)$(DSEP)AT_ParabolicCylinderFunction.c $(SRCDIR)$(DSEP)AT_RDD.c $(SRCDIR)$(DSEP)AT_SuccessiveConvolutions.c $(SRCDIR)$(DSEP)AT_Transport.c $(SRCDIR)$(DSEP)AT_Utils.c

LIBHOBJS  = $(INCLDIR)$(DSEP)AmTrack.h $(INCLDIR)$(DSEP)AT_Constants.h $(INCLDIR)$(DSEP)AT_FileOperations.h $(INCLDIR)$(DSEP)AT_Functions.h $(INCLDIR)$(DSEP)AT_GammaResponse.h $(INCLDIR)$(DSEP)AT_ParabolicCylinderFunction.h $(INCLDIR)$(DSEP)AT_RDD.h $(INCLDIR)$(DSEP)AT_SuccessiveConvolutions.h $(INCLDIR)$(DSEP)AT_Transport.h $(INCLDIR)$(DSEP)AT_Utils.h

LIBOBJS  = AmTrack.o AT_Constants.o AT_FileOperations.o AT_Functions.o AT_GammaResponse.o AT_ParabolicCylinderFunction.o AT_RDD.o AT_SuccessiveConvolutions.o AT_Transport.o AT_Utils.o


all:$(LIBOBJS)
		echo $(GSLPATH)$(DSEP)
		$(GCCDIR)gcc -L$(GSLPATH)$(DSEP)lib -shared $(LIBOBJS) -o $(NAMELIB) $(LFLAGS) 
		$(RMCMD) *.o


UI:AT_UI.o
		$(GCCDIR)gcc -L$(GSLPATH)$(DSEP)lib AT_UI.o -o $(NAMEEXE) $(LFLAGS) 
		$(RMCMD) *.o

clean:
		$(RMCMD) *.o
		$(RMCMD) $(NAMELIB)
		$(RMCMD) $(NAMEEXE)

$(LIBOBJS):$(LIBCOBJS) $(LIBHOBJS)
		$(GCCDIR)gcc -I$(INCLDIR) -I$(GSLPATH)$(DSEP)include $(CFLAGS) $(LIBCOBJS)


AT_UI.o:$(SRCDIR)$(DSEP)AT_UI.c 
		$(GCCDIR)gcc -I$(INCLDIR) -I$(GSLPATH)$(DSEP)include $(CFLAGS) $(SRCDIR)$(DSEP)AT_UI.c

