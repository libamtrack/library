cmake_minimum_required(VERSION 3.16)  # Align with the main CMakeLists.txt

project(amtrack_plots LANGUAGES C)

# Collect all C source files in the current directory
file(GLOB SOURCE_FILES "*.c")

# Define the executable
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# Enable C99 standard explicitly
target_compile_options(${PROJECT_NAME} PRIVATE -std=c99)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
    $<INSTALL_INTERFACE:include>  # Ensures headers are found correctly in installed builds
)

# Link against GSL, amtrack, and math library (`-lm`)
target_link_libraries(${PROJECT_NAME} PRIVATE amtrack GSL::gsl GSL::gslcblas)

# Ensure `libm` (`-lm`) is linked only when needed
if(NOT WIN32)
    target_link_libraries(${PROJECT_NAME} PRIVATE m)
endif()

# Install the binary to `bin/`
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
