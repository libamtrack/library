cmake_minimum_required(VERSION 2.8)

project(libamtrack-test)

set(CMAKE_C_FLAGS "-std=c99")

file(GLOB SOURCE_FILES "*.c")

include_directories( "../../include")

add_executable (${PROJECT_NAME} ${SOURCE_FILES} )

target_link_libraries( ${PROJECT_NAME} LINK_PUBLIC amtrack ${GSL_LIBRARIES})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
