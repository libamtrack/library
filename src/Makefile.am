lib_LTLIBRARIES = libamtrack.la
libamtrack_la_SOURCES = AT.c AT_Algorithms_CPP.c AT_Algorithms_IGK.c AT_Error.c AT_Constants.c \
          AT_GammaResponse.c AT_DataRange.c \
          AT_DataMaterial.c AT_DataParticle.c AT_ElectronRange.c \
          AT_NumericalRoutines.c AT_PhysicsRoutines.c AT_RDD.c AT_RDD_Simple.c AT_RDD_Tabulated.c\
          AT_RDD_ShellAveraged.c AT_RDD_ExtendedTarget.c AT_KatzModel.c AT_SuccessiveConvolutions.c \
          AT_Histograms.c AT_Algorithms_GSM.c AT_SPC.c AT_MultipleCoulombScattering.c AT_EnergyLoss.c \
          AT_StoppingPowerDataBethe.c AT_StoppingPowerDataFromFile.c \
          AT_StoppingPowerDataPSTAR.c AT_StoppingPowerDataICRU.c \
          AT_StoppingPowerData.c AT_StoppingPower.c
include_HEADERS= ../include/AT.h ../include/AT_Algorithms_CPP.h ../include/AT_Algorithms_IGK.h ../include/AT_Error.h ../include/AT_Constants.h \
 ../include/AT_GammaResponse.h ../include/AT_DataRange.h \
 ../include/AT_DataMaterial.h ../include/AT_DataParticle.h ../include/AT_ElectronRange.h \
 ../include/AT_NumericalRoutines.h ../include/AT_PhysicsRoutines.h ../include/AT_RDD.h ../include/AT_RDD_Simple.h  ../include/AT_RDD_Tabulated.h\
 ../include/AT_RDD_ShellAveraged.h ../include/AT_RDD_ExtendedTarget.h ../include/AT_KatzModel.h ../include/AT_SuccessiveConvolutions.h \
 ../include/AT_Histograms.h ../include/AT_Algorithms_GSM.h ../include/AT_SPC.h ../include/AT_EnergyLoss.h \
 ../include/AT_MultipleCoulombScattering.h \
 ../include/AT_StoppingPowerDataFromFile.h \
 ../include/AT_StoppingPowerDataBethe.h \
 ../include/AT_StoppingPowerDataPSTAR.h \
 ../include/AT_StoppingPowerDataICRU.h \
 ../include/AT_StoppingPowerData.h \
 ../include/AT_StoppingPower.h

libamtrack_la_CFLAGS =-I../include -I../thirdparty/cernlib -Wall -pedantic -std=c99 -fPIC $(GSL_CFLAGS) $(OPENMP_CFLAGS)
libamtrack_la_LDFLAGS = -version-info 0:0:0 $(GSL_LIBS) ../thirdparty/cernlib/libcernlibfuns.la -no-undefined
