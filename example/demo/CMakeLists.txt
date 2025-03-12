cmake_minimum_required(VERSION 3.16)  # Ensure consistency with main CMakeLists.txt

project(amtrack_demo LANGUAGES C)

# Collect all C source files in the current directory
file(GLOB SOURCE_FILES "*.c")

# Define the executable
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# Include directories (use `PRIVATE` since they are only needed at compile time)
target_include_directories(${PROJECT_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
    $<INSTALL_INTERFACE:include>  # Ensures installed projects find headers correctly
)

# Link against amtrack and GSL
target_link_libraries(${PROJECT_NAME} PRIVATE amtrack GSL::gsl GSL::gslcblas)

# Ensure the math library (`libm`) is only linked when needed
if(NOT WIN32)
    target_link_libraries(${PROJECT_NAME} PRIVATE m)
endif()

# Install the demo binary to `bin`
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
