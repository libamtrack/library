cmake_minimum_required(VERSION 3.7.2)
project(library C)
set(CMAKE_C_STANDARD 99)

set(GSL_INCLUDE_DIR "$ENV{GSL_INCLUDE_DIR}")
FIND_PACKAGE(GSL REQUIRED)

file(GLOB SOURCE_FILES "../../src/*.c")

include_directories( ${GSL_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/../../include )

add_library(at STATIC ${SOURCE_FILES})
install(TARGETS at DESTINATION ../../lib)

target_link_libraries(at ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})