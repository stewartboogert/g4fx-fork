

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(geant4fx)

include(ExternalProject)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
if(WITH_GEANT4_VIS)
  find_package(Geant4 REQUIRED gdml vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

#----------------------------------------------------------------------------
# Add external 3rd party libraries
#

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
find_package(CGAL 5 REQUIRED)
find_package(VTK REQUIRED)
find_package(Catch2 REQUIRED)
find_package(OPENCASCADE  COMPONENTS FoundationClasses ModelingAlgorithms
                                     ModelingData DataExchange Visualization
             REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})


file(GLOB geant4AppSources ${PROJECT_SOURCE_DIR}/src/geant4App/*.cc)
file(GLOB geant4AppHeaders ${PROJECT_SOURCE_DIR}/include/geant4App/*.hh)
add_library(geant4App SHARED ${geant4AppSources} ${geant4AppHeaders})
target_compile_options(geant4App PRIVATE -Werror )
target_link_libraries(geant4App ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES})


file(GLOB occCADSources ${PROJECT_SOURCE_DIR}/src/occCAD/*.cc)
file(GLOB occCADHeaders ${PROJECT_SOURCE_DIR}/include/occCAD/*.hh)
add_library(occCAD SHARED ${occCADSources} ${occCADHeaders})
target_compile_options(occCAD PRIVATE -Werror )
target_include_directories(occCAD PUBLIC ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(occCAD ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} ${OpenCASCADE_DataExchange_LIBRARIES})

#add_library(core ${sources} ${headers})
#target_compile_options(core PRIVATE -Werror )
#target_link_libraries(core ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES})
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(load_gdml ${PROJECT_SOURCE_DIR}/src/test/load_gdml.cc )
target_compile_options(load_gdml PRIVATE -Werror  )
target_link_libraries(load_gdml ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} geant4App)

add_executable(load_step ${PROJECT_SOURCE_DIR}/src/test/load_step.cc )
target_compile_options(load_step PRIVATE -Werror  )
target_link_libraries(load_step ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} occCAD)


#add_executable(tess_test ${PROJECT_SOURCE_DIR}/src/test/tess_test.cc )
#target_compile_options(tess_test PRIVATE -Werror  )
#target_link_libraries(tess_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

#add_executable(vol_mesh_test ${PROJECT_SOURCE_DIR}/src/test/vol_mesh_test.cc )
#target_compile_options(vol_mesh_test PRIVATE -Werror  )
#target_link_libraries(vol_mesh_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core )

#add_executable(raymarch_test ${PROJECT_SOURCE_DIR}/src/test/raymarch_test.cc )
#target_compile_options(raymarch_test PRIVATE -Werror  )
#target_link_libraries(raymarch_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

#add_executable(dicom_test ${PROJECT_SOURCE_DIR}/src/test/dicom_test.cc )
#target_compile_options(dicom_test PRIVATE -Werror  )
#target_link_libraries(dicom_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

#add_executable(tracking_test ${PROJECT_SOURCE_DIR}/src/test/tracking_test.cc )
#target_compile_options(dicom_test PRIVATE -Werror  )
#target_link_libraries(tracking_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

#add_executable(bsp_test ${PROJECT_SOURCE_DIR}/src/test/bsp_test.cc ${sources} ${headers})
#target_compile_options(bsp_test PRIVATE -Werror  )
#target_link_libraries(bsp_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)


add_executable(EmbreeTests ${PROJECT_SOURCE_DIR}/src/test/testEmbree.cpp)
target_link_libraries(EmbreeTests PRIVATE Catch2::Catch2WithMain)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build G01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

# copy over test files to build path (good for testing)
file(COPY ${PROJECT_SOURCE_DIR}/testfiles DESTINATION ${PROJECT_BINARY_DIR}/)

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(G01 DEPENDS load_gdml)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS load_gdml EmbreeTests DESTINATION bin)

# copy over test files to install path
INSTALL(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/testfiles/ DESTINATION ${INSTALL_PATH}/)
