cmake_minimum_required(VERSION 3.16)

# -----------------------------------------------------------------------------
# Retrieve Git version information (if available)
# -----------------------------------------------------------------------------
find_program(GIT_EXECUTABLE git)
if (GIT_EXECUTABLE)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
else()
    set(GIT_VERSION "0.0.0-unknown")
endif()

# Remove "v" prefix from tag if present (e.g., v1.2.3 → 1.2.3)
string(REGEX REPLACE "^v" "" GIT_VERSION "${GIT_VERSION}")

# Extract Major and Minor version numbers; default patch to 0
string(REGEX MATCH "([0-9]+)\\.([0-9]+)" VERSION_MATCH ${GIT_VERSION})
set(PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})
set(PROJECT_VERSION_MINOR ${CMAKE_MATCH_2})
set(PROJECT_VERSION_PATCH "0")

# Detect extra commit info (number of commits since tag and hash)
string(REGEX MATCH "[0-9]+\\.[0-9]+-([0-9]+)-g([a-f0-9]+)" GIT_COMMIT_MATCH ${GIT_VERSION})
if (CMAKE_MATCH_1)
    set(PROJECT_VERSION_EXTRA "+${CMAKE_MATCH_1}-g${CMAKE_MATCH_2}")
    set(PROJECT_GIT_HASH ${CMAKE_MATCH_2})
else()
    set(PROJECT_VERSION_EXTRA "")
    set(PROJECT_GIT_HASH "release")
endif()

# -----------------------------------------------------------------------------
# Define project with computed version
# -----------------------------------------------------------------------------
project(libamtrack VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
message(STATUS "Building libamtrack version ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}${PROJECT_VERSION_EXTRA}")

# -----------------------------------------------------------------------------
# Generate version header file and version info file
# -----------------------------------------------------------------------------
set(VERSION_HEADER_PATH "${CMAKE_BINARY_DIR}/include/AT_Version.h")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/AT_Version.h.in
    ${VERSION_HEADER_PATH}
    @ONLY
)
include_directories("${CMAKE_BINARY_DIR}/include")

file(WRITE "${CMAKE_BINARY_DIR}/VERSION" "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}-${PROJECT_GIT_HASH}")

# -----------------------------------------------------------------------------
# Create CMake package configuration files for versioning
# -----------------------------------------------------------------------------
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/libamtrackConfigVersion.cmake"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}
    COMPATIBILITY AnyNewerVersion
)


# -----------------------------------------------------------------------------
# Define options for development and examples
# -----------------------------------------------------------------------------
option(BUILD_DEV "Install development files (headers, CMake package config, etc.)" OFF)
option(BUILD_EXAMPLES "Build example projects" ON)

# # -----------------------------------------------------------------------------
# # Locate Dependencies
# # -----------------------------------------------------------------------------
# find_package(GSL REQUIRED)
# if (NOT GSL_FOUND)
#     message(FATAL_ERROR "GSL not found! Install it using vcpkg or your system package manager.")
# endif()

###############################################################################
# Configure GSL dependency
###############################################################################
include(FetchContent)
FetchContent_Declare(
  gsl
  GIT_REPOSITORY https://github.com/ampl/gsl.git
  GIT_TAG v2.7.0
)

set(GSL_BUILD_TESTS OFF CACHE BOOL "Disable GSL tests" FORCE)
set(GSL_BUILD_DOCS OFF CACHE BOOL "Disable GSL documentation" FORCE)
set(NO_AMPL_BINDINGS 1 CACHE BOOL "Disable AMPL bindings" FORCE)
set(GSL_DISABLE_TESTS 1 CACHE BOOL "Disable GSL tests" FORCE)
set(DOCUMENTATION OFF CACHE BOOL "Disable GSL documentation" FORCE)

FetchContent_MakeAvailable(gsl)

add_library(GSL::gsl ALIAS gsl)
add_library(GSL::gslcblas ALIAS gsl)

# -----------------------------------------------------------------------------
# Collect source and header files
# -----------------------------------------------------------------------------
file(GLOB SOURCE_FILES "src/*.c")
file(GLOB HEADER_FILES "include/*.h")

# -----------------------------------------------------------------------------
# Create the shared library target
# -----------------------------------------------------------------------------
add_library(amtrack SHARED ${SOURCE_FILES})
target_include_directories(amtrack PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${GSL_INCLUDE_DIRS}
)
target_link_libraries(amtrack PUBLIC GSL::gsl GSL::gslcblas)
target_compile_features(amtrack PUBLIC c_std_11)

# -----------------------------------------------------------------------------
# Windows-specific settings for examples (getopt dependency)
# -----------------------------------------------------------------------------
if (WIN32 AND BUILD_EXAMPLES)
    find_library(GETOPT_LIB NAMES getopt PATHS ${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-windows/lib)
    find_path(GETOPT_INCLUDE_DIR getopt.h PATHS ${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-windows/include)
    if (NOT GETOPT_LIB)
        message(FATAL_ERROR "getopt.lib not found! Make sure it's installed via vcpkg.")
    endif()
    target_include_directories(amtrack PUBLIC ${GETOPT_INCLUDE_DIR})
    target_link_libraries(amtrack PUBLIC ${GETOPT_LIB})
endif()

# -----------------------------------------------------------------------------
# Specify output directories for targets
# -----------------------------------------------------------------------------
set_target_properties(amtrack PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Export all symbols on Windows to simplify DLL creation
if (WIN32)
    set_target_properties(amtrack PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()

# -----------------------------------------------------------------------------
# Installation of the libamtrack target
# -----------------------------------------------------------------------------
option(LIBAMTRACK_INSTALL "Install libamtrack targets" ON)
if (LIBAMTRACK_INSTALL)
    install(FILES "${VERSION_HEADER_PATH}" DESTINATION include)
    install(FILES "${CMAKE_BINARY_DIR}/VERSION" DESTINATION share/libamtrack)
    install(FILES "${CMAKE_BINARY_DIR}/libamtrackConfigVersion.cmake"
        DESTINATION lib/cmake/libamtrack
    )
    install(TARGETS amtrack
        EXPORT amtrackTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
    )
endif()

# -----------------------------------------------------------------------------
# Installation of development files (headers and CMake config)
# -----------------------------------------------------------------------------
if (BUILD_DEV)
    message(STATUS "Development files will be installed.")
    install(DIRECTORY include/ DESTINATION include)

    file(WRITE "${CMAKE_BINARY_DIR}/libamtrackConfig.cmake" "
@PACKAGE_INIT@
include(\"\${CMAKE_CURRENT_LIST_DIR}/amtrackTargets.cmake\")
")
    install(FILES "${CMAKE_BINARY_DIR}/libamtrackConfig.cmake"
        DESTINATION lib/cmake/libamtrack
    )
    install(EXPORT amtrackTargets
        FILE amtrackTargets.cmake
        NAMESPACE amtrack::
        DESTINATION lib/cmake/libamtrack
    )
endif()

# -----------------------------------------------------------------------------
# Add examples and tests if enabled
# -----------------------------------------------------------------------------
if (BUILD_EXAMPLES)
    add_subdirectory("example/demo")
    add_subdirectory("example/basic_plots")
    add_subdirectory("test/C")

    foreach(target amtrack_demo amtrack_plots amtrack_test)
        if (TARGET ${target})
            add_dependencies(${target} amtrack)
            target_link_libraries(${target} PRIVATE amtrack)
            target_link_directories(${target} PRIVATE ${CMAKE_BINARY_DIR}/lib)
        endif()
    endforeach()
endif()

# -----------------------------------------------------------------------------
# Include CPack for packaging support
# -----------------------------------------------------------------------------
include(CPack)