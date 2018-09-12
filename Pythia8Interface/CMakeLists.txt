project(Pythia8Interface LANGUAGES CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../HTTAnalysis)

add_library(Pythia8 SHARED IMPORTED)
set_property(TARGET Pythia8 PROPERTY IMPORTED_LOCATION $ENV{PYTHIA8_DIR}/lib/libpythia8.so)

add_library(Pythia8Interface SHARED Pythia8Interface)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(Pythia8Interface PUBLIC
                           ${CMAKE_SOURCE_DIR}/include
                           $ENV{PYTHIA8_DIR}/include
			   ${ROOTINCDIR})
  
# Depend on a library that we defined in the top-level file
target_link_libraries(Pythia8Interface RootAnalysisCore HTTAnalysis boost ${ROOT_LIBRARIES} ${PYTHIA8_LIBRARY})

# Define an executable
add_executable(generateEvents generateEvents.cc)

target_include_directories(generateEvents PUBLIC
    ${CMAKE_SOURCE_DIR}/include   
    ${ROOTINCDIR})

target_link_libraries(generateEvents Pythia8Interface RootAnalysisCore HTTAnalysis ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} ${PYTHIA8_LIBRARY} Proof RooFitCore RooFit EG EGPythia8 Pythia8)


