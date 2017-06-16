#----------------------------------------------------------------------------
# Setup the project
#

#Enforce cmake version
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
#cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

#Set project name
project(mtestgeant4)

#Choose C++ standard
set(CMAKE_CXX_STANDARD 11) #Use with cmake version >= 3.1
#add_compile_options(-std=c++11) #Use if with cmake version < 3.1


#----------------------------------------------------------------------------
# Setup GEANT4
#

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Setup ROOT
#

#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT PATHS $ENV{ROOTSYS}/etc/cmake )# REQUIRED COMPONENTS RIO)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 and ROOT libraries
#

add_executable(example example.cc ${sources} ${headers})
target_link_libraries(example ${Geant4_LIBRARIES})
target_link_libraries(example ${ROOT_LIBRARIES})


#----------------------------------------------------------------------------
# Copy all scripts to the build directory. This is so that we can run the 
# executable directly because it relies on these scripts being in the current 
# working directory.
#

set(EXAMPLE_SCRIPTS
  example.mac
  vis.mac
  )

foreach(_script ${EXAMPLE_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#

add_custom_target(mtestgeant4 DEPENDS example)


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#

install(TARGETS example DESTINATION bin)


