# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.10.2)
project (pipeline LANGUAGES NONE)


set(CMAKE_BUILD_TYPE Release)

# force out-of-source build
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
        message(FATAL_ERROR "In-source builds not allowed.
        Please make a new directory (called a build directory) and run CMake from there.
        You may need to remove CMakeCache.txt." )
endif()


option(AP_USE_CONTAINERS "Allow the usage of singularity images (profile Containers)" ON)
option(AP_USE_TOOLSPATH "Allow the usage of tools installed in existing directories (profile toolsPath) " OFF)


install(DIRECTORY . DESTINATION ${CMAKE_INSTALL_PREFIX} PATTERN "conf/templates" EXCLUDE PATTERN ".git*" EXCLUDE PATTERN "CMakeLists.txt" EXCLUDE)

add_subdirectory (conf/templates)

