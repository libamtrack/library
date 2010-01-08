#****************************************************************************
#
# Makefile for libamtrack, lgrzanka, 2009/12/27
#
#****************************************************************************

CFLAGS    = -Wall -c -O3 -fPIC
LDFLAGS    = -lm -lgsl -lgslcblas


#LIBCOBJS  = $(SRCDIR)$(DSEP)AmTrack.c $(SRCDIR)$(DSEP)AT_Constants.c $(SRCDIR)$(DSEP)AT_FileOperations.c $(SRCDIR)$(DSEP)AT_Functions.c $(SRCDIR)$(DSEP)AT_GammaResponse.c $(SRCDIR)$(DSEP)AT_ParabolicCylinderFunction.c $(SRCDIR)$(DSEP)AT_RDD.c $(SRCDIR)$(DSEP)AT_SuccessiveConvolutions.c $(SRCDIR)$(DSEP)AT_Transport.c $(SRCDIR)$(DSEP)AT_Utils.c

#LIBCOBJS = $(SRCDIR)$(DSEP)*.c

#LIBHOBJS  = $(INCLDIR)$(DSEP)AmTrack.h $(INCLDIR)$(DSEP)AT_Constants.h $(INCLDIR)$(DSEP)AT_FileOperations.h $(INCLDIR)$(DSEP)AT_Functions.h $(INCLDIR)$(DSEP)AT_GammaResponse.h $(INCLDIR)$(DSEP)AT_ParabolicCylinderFunction.h $(INCLDIR)$(DSEP)AT_RDD.h $(INCLDIR)$(DSEP)AT_SuccessiveConvolutions.h $(INCLDIR)$(DSEP)AT_Transport.h $(INCLDIR)$(DSEP)
#AT_Utils.h

#LIBHOBJS  = $(INCLDIR)$(DSEP)*.h


LIBOBJFILES  = AT_Constants.o AT_FileOperations.o AT_Functions.o AT_GammaResponse.o AT_ParabolicCylinderFunction.o AT_RDD.o AT_SuccessiveConvolutions.o AT_Transport.o AT_Utils.o AmTrack.o 


OBJFILES := $(patsubst %.c,%.o,$(wildcard *.c))

all: libAT.so
		echo "Linking"


clean:
		rm -f *.o


libAT.so: $(LIBOBJFILES)
		gcc $(LDFLAGS) -I../include/ -o libAT.so $(LIBOBJFILES)


$(OBJFILES): *.c
		gcc $(CFLAGS) -I../include/ -c -o $@ $<


