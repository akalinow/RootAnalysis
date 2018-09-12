cmake_minimum_required (VERSION 2.6)
project (RootAnalysis)
set( CMAKE_VERBOSE_MAKEFILE on )

# Include Boost as an imported target
find_package(Boost REQUIRED COMPONENTS regex)
add_library(boost INTERFACE IMPORTED)
set_property(TARGET boost PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})

##https://root.cern.ch/how/integrate-root-my-project-cmake
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO GenVector Physics)
include(${ROOT_USE_FILE})

if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS 6.0.0)
  SET( CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS} -Wall -Werror -O3 -std=c++11 -fopenmp")
  MESSAGE("Using -std=c++11 for compiler older than 6.0.0")
else()
  SET( CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS} -Wall -Werror -O3 -std=c++14 -fopenmp")
  MESSAGE("Using -std=c++14 for compiler newer or equal to 6.0.0")
endif()

add_subdirectory(src)
add_subdirectory(HTTAnalysis)
add_subdirectory(Pythia8Interface)
#add_subdirectory(SVFitAnalysis)
