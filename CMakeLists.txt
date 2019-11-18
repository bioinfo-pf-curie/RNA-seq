# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.10.2)
project (RNA-seq LANGUAGES NONE)


set(CMAKE_BUILD_TYPE Release)
set(ap_install_prefix_subdirectory "" CACHE STRING "if not null, the pipeline will be installed in CMAKE_INSTALL_PREFIX/ap_install_prefix_subdirectory")
set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/${ap_install_prefix_subdirectory})

# force out-of-source build
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
        message(FATAL_ERROR "In-source builds not allowed.
        Please make a new directory (called a build directory) and run CMake from there.
        You may need to remove CMakeCache.txt." )
endif()


option(ap_use_containers "Allow the usage of singularity images (profile Containers)" ON)
option(ap_use_toolspath "Allow the usage of tools installed in existing directories (profile toolsPath) " OFF)


install(DIRECTORY . DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS PATTERN "conf/templates" EXCLUDE PATTERN ".git*" EXCLUDE PATTERN "CMakeLists.txt" EXCLUDE PATTERN "*.example" EXCLUDE PATTERN)

add_subdirectory (conf/templates)

