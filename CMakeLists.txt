cmake_minimum_required (VERSION 3.2)
project (RootAnalysis)
set( CMAKE_VERBOSE_MAKEFILE on )

# Include Boost as an imported target
#On lxplus one needs a newer boost version thant available by default
find_package(Boost REQUIRED COMPONENTS regex)
add_library(boost INTERFACE IMPORTED)
set_property(TARGET boost PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
MESSAGE("Boost directory: " ${Boost_LIBRARY_DIRS})

##https://root.cern.ch/how/integrate-root-my-project-cmake
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO GenVector Physics RooFitCore RooFit)
include(${ROOT_USE_FILE})

if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS 6.0.0)
  SET( CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS} -Wall -Werror -O3 -std=c++11 -fopenmp")
  MESSAGE("Using -std=c++11 for compiler older than 6.0.0")
elseif("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS 7.3.0)
  SET( CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS} -Wall -Werror -O3 -std=c++14 -fopenmp")
  MESSAGE("Using -std=c++14 for compiler older than 7.3.0")
else()
  SET( CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS} -Wall -Werror -O3 -std=c++17 -fopenmp")
  MESSAGE("Using -std=c++17 for compiler newer or equal 7.3.0")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX $ENV{PWD}  CACHE PATH "..." FORCE)
endif()

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

MESSAGE("Installation prefix set to " ${CMAKE_INSTALL_PREFIX} )

add_subdirectory(Core/src)
add_subdirectory(OMTF)

