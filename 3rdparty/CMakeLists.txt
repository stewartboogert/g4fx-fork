# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(geant4fxExternal)

include(ExternalProject)

#----------------------------------------------------------------------------
# Add external 3rd party libraries
#

message(STATUS ${CMAKE_INSTALL_PREFIX}/3rdparty)

ExternalProject_Add(catch2
        URL ${PROJECT_SOURCE_DIR}/../3rdparty/Catch2-3.3.2.tar.gz
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/catch2
        CMAKE_ARGS  -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        BUILD_COMMAND make install -j4)

ExternalProject_Add(tbb
        URL ${PROJECT_SOURCE_DIR}/../3rdparty/oneTBB-2021.9.0.tar.gz
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tbb
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        BUILD_COMMAND make install -j4)

ExternalProject_Add(embree
        URL ${PROJECT_SOURCE_DIR}/../3rdparty/embree-4.0.1.tar.gz
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/embree
        CMAKE_ARGS -DTBB_DIR=${CMAKE_INSTALL_PREFIX}/3rdparty/tbb -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        BUILD_COMMAND make install -j4
        )

add_dependencies(embree tbb)
