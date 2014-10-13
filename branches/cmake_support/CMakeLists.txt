cmake_minimum_required(VERSION 2.8)

project(amtrack)

# TODO: add support for Debug and Release targets

# TODO: add generation of R wrapper, useful hints:
# http://stackoverflow.com/questions/4135936/ide-setup-for-package-development-with-c-code-integrated


# TODO: add reading of svn version number
#include(FindSubversion)
#if( EXISTS ${CMAKE_SOURCE_DIR}/.svn/ )
#    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} libamtrack)
#endif()

# dependencies - i.e. GSL
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
FIND_PACKAGE(GSL REQUIRED)


if (GSL_FOUND)
    # source and header files
    file(GLOB SOURCE_FILES "src/*.c")
    file(GLOB HEADER_FILES "include/*.h")

    # specify include directories
    include_directories( "include" ${GSL_INCLUDE_DIRS})

    # specify target type
    add_library (${PROJECT_NAME} SHARED ${SOURCE_FILES} )

    install(TARGETS ${PROJECT_NAME} DESTINATION lib)

    target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})

    # Build sub-projects
    add_subdirectory( "example/demo")
    add_subdirectory( "example/basic_plots")
    add_subdirectory( "test/C")

    # TODO: build a CPack driven installer package
    include (InstallRequiredSystemLibraries)
    #set (CPACK_RESOURCE_FILE_LICENSE
    #     "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
    #set (CPACK_PACKAGE_VERSION_MAJOR "${Tutorial_VERSION_MAJOR}")
    #set (CPACK_PACKAGE_VERSION_MINOR "${Tutorial_VERSION_MINOR}")
    include (CPack)

else (GSL_FOUND)
    message( FATAL_ERROR "FindGSL: Could not find GSL headers or library" )
endif (GSL_FOUND)


