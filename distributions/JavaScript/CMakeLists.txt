cmake_minimum_required(VERSION 3.7.2)
project(library C)
set(CMAKE_C_STANDARD 99)

FIND_PACKAGE(GSL REQUIRED)

file(GLOB SOURCE_FILES "../../src/*.c")

include_directories( ${CMAKE_SOURCE_DIR}/../../include ${GSL_INCLUDE_DIRS})

add_library(at STATIC ${SOURCE_FILES})
install(TARGETS at DESTINATION ../../lib)
target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
