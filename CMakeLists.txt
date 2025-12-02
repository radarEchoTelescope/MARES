#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(MS_single_event)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
#----------------------------------------------------------------------------

# Setup include directory for this project
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/lib/IceRayTracing/namespace/woROOT/)
include_directories(${PROJECT_SOURCE_DIR}/lib/libconfig/bin/include/)

find_package(PkgConfig)
find_package(GSL REQUIRED)
pkg_check_modules(LIBCONFIG_PLUSPLUS libconfig++)
if(NOT LIBCONFIG_PLUSPLUS_FOUND)
    set(LIBCONFIG_PLUSPLUS ${PROJECT_SOURCE_DIR}/lib/libconfig/bin/lib/libconfig++.so)
    message( "Using local libconfig at: " ${LIBCONFIG_PLUSPLUS})
endif()



#----------------------------------------------------------------------------
# Build MARES library
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs

#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

set( LIB_SOURCES
    ${PROJECT_SOURCE_DIR}/lib/IceRayTracing/namespace/woROOT/IceRayTracing.cc
    ${PROJECT_SOURCE_DIR}/src/antenna.cc
    ${PROJECT_SOURCE_DIR}/src/cascade.cc
    ${PROJECT_SOURCE_DIR}/src/cascade1D.cc
    ${PROJECT_SOURCE_DIR}/src/scatter.cc
    ${PROJECT_SOURCE_DIR}/src/scatter_propagation.cc
    ${PROJECT_SOURCE_DIR}/src/settings_math.cc
    ${PROJECT_SOURCE_DIR}/src/settings.cc
    ${PROJECT_SOURCE_DIR}/src/settings_defaults.cc

)

set( LIB_HEADERS
    ${PROJECT_SOURCE_DIR}/lib/IceRayTracing/namespace/woROOT/IceRayTracing.hh
    ${PROJECT_SOURCE_DIR}/include/antenna.hh
    ${PROJECT_SOURCE_DIR}/include/cascade.hh
    ${PROJECT_SOURCE_DIR}/include/cascade1D.hh
    ${PROJECT_SOURCE_DIR}/include/scatter.hh
    ${PROJECT_SOURCE_DIR}/include/settings.hh
    ${PROJECT_SOURCE_DIR}/include/settings_units_constants.hh
    ${PROJECT_SOURCE_DIR}/lib/libconfig/bin/include/libconfig.h++

)

add_library(libMARES ${LIB_SOURCES} ${LIB_HEADERS})

#----------------------------------------------------------------------------
# Add the executable, and link it to the MARES library
#
# message( "This is the message " ${PROJECT_SOURCE_DIR})

set( EXEC_SOURCES
    ${PROJECT_SOURCE_DIR}/MARES_single_point.cpp
    ${PROJECT_SOURCE_DIR}/MARES_single_event.cpp
    ${PROJECT_SOURCE_DIR}/MARES_beam_event.cpp
    ${PROJECT_SOURCE_DIR}/MARES_single_event_extern.cpp
    ${PROJECT_SOURCE_DIR}/MARES.cpp
 )


foreach( testsourcefile ${EXEC_SOURCES} )
    # Grab the filename to produce the executable
    get_filename_component(testname ${testsourcefile} NAME_WE)
    # message( "This is the executable's name " ${testname})
    add_executable( ${testname} ${testsourcefile} )
    # Make sure MARES is linked to each app
    target_link_libraries( ${testname} 
        libMARES
       # ${LIBCONFIG_PLUSPLUS_LIBRARIES}
        ${LIBCONFIG_PLUSPLUS} 
        GSL::gsl
    )
endforeach( testsourcefile ${EXEC_SOURCES} )
