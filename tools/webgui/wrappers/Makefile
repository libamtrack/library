#****************************************************************************
# Makefile for webgui wrappers
#****************************************************************************

############################################################
# USER SETTINGS
############################################################
## OS settings detection
TRG       = $(shell uname)
ifeq ($(TRG),Linux)
OS = Linux
else
OS = Windows
endif
## Please set the correct path to your GSL installation here:
#  most likely:
ifeq ($(OS),Linux)
GSLPATH   = /usr
else
# GSLPATH   = C:\Program Files\GnuWin32
GSLPATH   = C:\Programme\GnuWin32
endif
## Please specify the path to your MinGW installation (Windows only)
ifeq ($(OS),Linux)
MINGWPATH = 
else
#  MINGWPATH = C:\Program Files\MinGW
MINGWPATH = "C:\Programme\MinGW"
endif
############################################################

CFLAGS    = -Wall -c -O3
LFLAGS    = -lm -lgsl -lgslcblas

ifeq ($(OS),Linux)
STATICLIB   = libamtrack.a
LIBDIR    = ../../../lib/
DSEP      = /
INCLDIR   = ../../../include/
RMCMD     = rm
GCC       = gcc 
GSLINCLUDE = "$(GSLPATH)/include"
GSLLIB     = "$(GSLPATH)/lib"
else
STATICLIB    = libamtrack.lib
LIBDIR    = ..\..\..\lib\
DSEP       = \\
INCLDIR    = ..\..\..\include\
RMCMD      = del
GCC        = "$(MINGWPATH)\bin\gcc.exe"
GSLINCLUDE = "$(GSLPATH)\include"
GSLLIB     = "$(GSLPATH)\lib"
endif

SOURCES=AT_CPPSC.c AT_CSDA_range.c AT_D_RDD_GY.c AT_E_from_beta.c AT_E_from_gamma.c AT_fluence_from_dose.c AT_LET_from_E.c AT_Katz.c
OBJECTS=$(SOURCES:.c=.o)
EXECS=$(SOURCES:.c=.exe)

all: $(EXECS) $(SOURCES)


%.exe: %.o
		$(GCC) -L$(GSLLIB) -L$(LIBDIR) $< $(LFLAGS) -lamtrack -o $@


%.o: %.c
		$(GCC) -I$(GSLINCLUDE) -I$(INCLDIR) $(CFLAGS) $< -o $@


clean:
		- $(RMCMD) $(EXECS)
