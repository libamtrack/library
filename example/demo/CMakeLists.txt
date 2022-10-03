cmake_minimum_required(VERSION 2.8)

project(amtrack_demo)

file(GLOB SOURCE_FILES "*.c")

include_directories( "../../include")

add_executable (${PROJECT_NAME} ${SOURCE_FILES} )

target_link_libraries( ${PROJECT_NAME} LINK_PUBLIC amtrack ${GSL_LIBRARIES})

# Visual Studio does not need or want you to explicitly request linking the math library.
# You must avoid adding it as a link library when building for Windows.
# see https://stackoverflow.com/questions/54935559/linking-math-library-in-cmake-file-on-windows-and-linux
IF (NOT WIN32)
  target_link_libraries (${PROJECT_NAME} m)
ENDIF()


install(TARGETS ${PROJECT_NAME} DESTINATION bin)
