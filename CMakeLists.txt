#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MS_single_event)

#----------------------------------------------------------------------------

# Setup include directory for this project
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
# message( "This is the message " ${PROJECT_SOURCE_DIR})
add_executable(MS_single_event MS_single_event.cpp ${sources} ${headers})
# target_link_libraries(MS_single_event gsl)
target_link_libraries(MS_single_event)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS ice_shelf DESTINATION bin)