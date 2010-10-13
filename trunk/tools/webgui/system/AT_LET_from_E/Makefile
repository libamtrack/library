#****************************************************************************
# Makefile for AT_LET_from_E
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

CFLAGS    = -Wall -c -O3
LFLAGS    = -lm -lgsl -lgslcblas

ifeq ($(OS),Linux)
STATICLIB   = libamtrack.a
LIBDIR    = ../../../../lib/
DSEP      = /
INCLDIR   = ../../../../include/
RMCMD     = rm
GCC       = gcc 
GSLINCLUDE = "$(GSLPATH)/include"
GSLLIB     = "$(GSLPATH)/lib"
else
STATICLIB    = libamtrack.lib
LIBDIR    = ..\..\..\..\lib\
DSEP       = \\
INCLDIR    = ..\..\..\..\include\
RMCMD      = del
GCC        = "$(MINGWPATH)\bin\gcc.exe"
GSLINCLUDE = "$(GSLPATH)\include"
GSLLIB     = "$(GSLPATH)\lib"
endif


AT_LET_from_E.o:AT_LET_from_E.c 
		$(GCC) -I$(GSLINCLUDE) -I$(INCLDIR) $(CFLAGS) AT_LET_from_E.c


all:AT_LET_from_E.o
		$(GCC) -L$(GSLLIB) -L$(LIBDIR) AT_LET_from_E.o $(LFLAGS) -lamtrack -o AT_LET_from_E.exe

static:AT_LET_from_E.c
		$(GCC) -I$(INCLDIR) -L$(GSLLIB) AT_LET_from_E.c $(LIBDIR)$(STATICLIB) $(LFLAGS) -o AT_LET_from_E.exe

clean:
		- $(RMCMD) *.o
		- $(RMCMD) AT_LET_from_E.exe
