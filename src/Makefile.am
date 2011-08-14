lib_LTLIBRARIES = libamtrack.la
libamtrack_la_SOURCES = AT.c AT_Algorithms_CPP.c AT_Algorithms_IGK.c AT_Error.c AT_Constants.c \
          AT_GammaResponse.c AT_DataStoppingPower.c AT_DataRange.c AT_DataLET.c \
          AT_DataMaterial.c AT_DataParticle.c AT_ElectronRange.c \
          AT_NumericalRoutines.c AT_PhysicsRoutines.c AT_RDD.c AT_RDD_Simple.c \
          AT_RDD_ShellAveraged.c AT_RDD_ExtendedTarget.c AT_KatzModel.c AT_SuccessiveConvolutions.c \
          AT_Histograms.c AT_Algorithms_GSM.c AT_SPC.c
include_HEADERS= ../include/AT.h ../include/AT_Algorithms_CPP.h ../include/AT_Algorithms_IGK.h ../include/AT_Error.h ../include/AT_Constants.h \
 ../include/AT_GammaResponse.h ../include/AT_DataStoppingPower.h ../include/AT_DataRange.h ../include/AT_DataLET.h \
 ../include/AT_DataMaterial.h ../include/AT_DataParticle.h ../include/AT_ElectronRange.h \
 ../include/AT_NumericalRoutines.h ../include/AT_PhysicsRoutines.h ../include/AT_RDD.h ../include/AT_RDD_Simple.h \
 ../include/AT_RDD_ShellAveraged.h ../include/AT_RDD_ExtendedTarget.h ../include/AT_KatzModel.h ../include/AT_SuccessiveConvolutions.h \
 ../include/AT_Histograms.h ../include/AT_Algorithms_GSM.h ../include/AT_SPC.h
libamtrack_la_CFLAGS =-I../include -Wall -std=gnu99 -fPIC
libamtrack_la_LDFLAGS = -version-info 0:0:0

bin_PROGRAMS = amtrack_test
amtrack_test_SOURCES = AT_Test.c
amtrack_test_CFLAGS =-I../include -Wall -std=gnu99 -fPIC
amtrack_test_LDADD = libamtrack.la