cmake_minimum_required(VERSION 3.7.2)
project(library C)
set(CMAKE_C_STANDARD 99)

# FIND_PACKAGE(GSL REQUIRED) - doesn't work with `emcmake cmake`, while working with pure `cmake`

file(GLOB SOURCE_FILES "../../src/*.c")

include_directories( ${GSL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/../../include )

add_library(at STATIC ${SOURCE_FILES})
install(TARGETS at DESTINATION ../../lib)
