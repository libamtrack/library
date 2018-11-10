cmake_minimum_required(VERSION 3.7.2)
project(library C)

set(CMAKE_C_STANDARD 99)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(/home/osboxes/gsl-latest/gsl-2.5/)
add_library(at STATIC
        src/AT.c
        src/AT_Algorithms_CPP.c
        src/AT_Algorithms_GSM.c
        src/AT_Algorithms_IGK.c
        src/AT_CernlibFuns.c
        src/AT_Constants.c
        src/AT_DataMaterial.c
        src/AT_DataParticle.c
        src/AT_DataRange.c
        src/AT_ElectronRange.c
        src/AT_EnergyLoss.c
        src/AT_Error.c
        src/AT_GammaResponse.c
        src/AT_Histograms.c
        src/AT_KatzModel.c
        src/AT_MultipleCoulombScattering.c
        src/AT_NumericalRoutines.c
        src/AT_PhysicsRoutines.c
        src/AT_RDD.c
        src/AT_RDD_ExtendedTarget.c
        src/AT_RDD_ShellAveraged.c
        src/AT_RDD_Simple.c
        src/AT_RDD_Tabulated.c
        src/AT_SPC.c
        src/AT_StoppingPower.c
        src/AT_StoppingPowerData.c
        src/AT_StoppingPowerDataBethe.c
        src/AT_StoppingPowerDataFromFile.c
        src/AT_StoppingPowerDataICRU.c
        src/AT_StoppingPowerDataPSTAR.c
        src/AT_SuccessiveConvolutions.c)
install(TARGETS at DESTINATION lib)
include_directories(${CMAKE_SOURCE_DIR}/inc)
