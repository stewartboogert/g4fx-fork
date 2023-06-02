

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
ExternalProject_Add(tbb URL 
  ${PROJECT_SOURCE_DIR}/3rdparty/oneTBB-2021.9.0.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tbb
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/3rdparty
  BUILD_COMMAND make install -j2)

ExternalProject_Add(embree URL 
  ${PROJECT_SOURCE_DIR}/3rdparty/embree-4.0.1.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/embree
  CMAKE_ARGS -DTBB_DIR=${CMAKE_INSTALL_PREFIX}/3rdparty/tbb -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/3rdparty
  BUILD_COMMAND make install -j2
  )

ADD_DEPENDENCIES(embree tbb)
#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
find_package(CGAL 5 REQUIRED)
find_package(VTK REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


add_library(core ${sources} ${headers})
target_compile_options(core PRIVATE -Werror )
target_link_libraries(core ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES})
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(load_gdml load_gdml.cc )
target_compile_options(load_gdml PRIVATE -Werror  )
target_link_libraries(load_gdml ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

add_executable(tess_test tess_test.cc )
target_compile_options(tess_test PRIVATE -Werror  )
target_link_libraries(tess_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

add_executable(vol_mesh_test vol_mesh_test.cc )
target_compile_options(vol_mesh_test PRIVATE -Werror  )
target_link_libraries(vol_mesh_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

add_executable(raymarch_test raymarch_test.cc )
target_compile_options(raymarch_test PRIVATE -Werror  )
target_link_libraries(raymarch_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

add_executable(dicom_test dicom_test.cc )
target_compile_options(dicom_test PRIVATE -Werror  )
target_link_libraries(dicom_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

add_executable(tracking_test tracking_test.cc )
target_compile_options(dicom_test PRIVATE -Werror  )
target_link_libraries(tracking_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

add_executable(bsp_test bsp_test.cc ${sources} ${headers})
target_compile_options(bsp_test PRIVATE -Werror  )
target_link_libraries(bsp_test ${Geant4_LIBRARIES} CGAL::CGAL ${VTK_LIBRARIES} core)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build G01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(GDML_FILES
    T001_Union.gdml T002_Subtraction.gdml T003_Intersection.gdml T004_UnionN.gdml T005_solids.gdml T006_torus.geo T006_torus.mesh
    T006_torus.off T006_torus.stl T006_torus.vtk T007_teapot.obj T007_teapot.stl vis.mac)

foreach(_script ${GDML_FILES})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(G01 DEPENDS load_gdml)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS load_gdml DESTINATION bin)

