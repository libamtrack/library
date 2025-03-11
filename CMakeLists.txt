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

# Ensure GSL is found
if (NOT GSL_FOUND)
    message(FATAL_ERROR "GSL not found! Install it using vcpkg or system package manager.")
endif()

# Include directories
target_include_directories(amtrack PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${GSL_INCLUDE_DIRS}
)

# Link against GSL and GSL CBLAS
target_link_libraries(amtrack PUBLIC GSL::gsl GSL::gslcblas)

# On Windows, explicitly link `getopt`
if (WIN32)
    find_library(GETOPT_LIB NAMES getopt PATHS ${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-windows/lib)
    find_path(GETOPT_INCLUDE_DIR getopt.h PATHS ${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-windows/include)

    if (NOT GETOPT_LIB)
        message(FATAL_ERROR "getopt.lib not found! Make sure it's installed via vcpkg.")
    endif()

    target_include_directories(amtrack PUBLIC ${GETOPT_INCLUDE_DIR})
    target_link_libraries(amtrack PUBLIC ${GETOPT_LIB})
endif()

# Require at least C11
target_compile_features(amtrack PUBLIC c_std_11)

# Ensure `amtrack` library is built before executables
add_dependencies(amtrack amtrack)

# Set correct output directories for the library
set_target_properties(amtrack PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

if (WIN32)
    # Export all symbols on Windows
    set_target_properties(amtrack PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()

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
    install(DIRECTORY include/ DESTINATION include)

    # Generate and install CMake package config file inline
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

# Build sub-projects (executables)
add_subdirectory("example/demo")
add_subdirectory("example/basic_plots")
add_subdirectory("test/C")

# Ensure executables link to amtrack and depend on it
foreach(target amtrack_demo amtrack_plots amtrack_test)
    if (TARGET ${target})
        add_dependencies(${target} amtrack)
        target_link_libraries(${target} PRIVATE amtrack)
        target_link_directories(${target} PRIVATE ${CMAKE_BINARY_DIR}/lib)
    endif()
endforeach()

# CPack Installer (optional)
include(CPack)
