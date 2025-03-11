cmake_minimum_required(VERSION 3.16)

project(amtrack LANGUAGES C)  # Explicitly set language to C

# Add option for building development files
option(BUILD_DEV "Install development files (headers, CMake package config, etc.)" OFF)

# Dependencies
find_package(GSL REQUIRED)

# Collect source and header files
file(GLOB SOURCE_FILES "src/*.c")
file(GLOB HEADER_FILES "include/*.h")

# Define the shared library
add_library(amtrack SHARED ${SOURCE_FILES})

# Include directories
target_include_directories(amtrack PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${GSL_INCLUDE_DIRS}
)

# Link against GSL and GSL CBLAS
target_link_libraries(amtrack PUBLIC GSL::gsl GSL::gslcblas)

# Require at least C11
target_compile_features(amtrack PUBLIC c_std_11)

# Install shared library
install(TARGETS amtrack
    EXPORT amtrackTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

# Install development files if requested
if(BUILD_DEV)
    message(STATUS "Development files will be installed.")

    # Install headers
    install(DIRECTORY include/ DESTINATION include)

    # Generate and install CMake package config file inline (instead of an external .in file)
    file(WRITE "${CMAKE_BINARY_DIR}/libamtrackConfig.cmake" "
@PACKAGE_INIT@
include(\"\${CMAKE_CURRENT_LIST_DIR}/amtrackTargets.cmake\")
")

    install(FILES "${CMAKE_BINARY_DIR}/libamtrackConfig.cmake"
        DESTINATION lib/cmake/libamtrack
    )

    # Install CMake export file
    install(EXPORT amtrackTargets
        FILE amtrackTargets.cmake
        NAMESPACE amtrack::
        DESTINATION lib/cmake/libamtrack
    )
endif()

# Build sub-projects
add_subdirectory("example/demo")
add_subdirectory("example/basic_plots")
add_subdirectory("test/C")

# CPack Installer (optional)
include(CPack)
