project(HTTAnalysis LANGUAGES CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(G__HTTEvent HTTEvent.h LINKDEF LinkDef.h)

add_library(HTTAnalysis SHARED ChannelSpecifics.cc  EventProxyHTT.cc 
  HTTAnalyzer.cc PropertyEnum.cc
  HTTAnalyzerCPmethods.cc HTTAnalyzerUtilityMethods.cc HTTHistograms.cc HTTSynchNTuple.cc
  MuTauSpecifics.cc TauTauSpecifics.cc Tools.cc BTagCalibrationStandalone.cpp HTTEvent.cxx
  MLAnalyzer.cc
  G__HTTEvent.cxx)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(HTTAnalysis PUBLIC
    ${CMAKE_SOURCE_DIR}/include   
    ${ROOTINCDIR})
  
# Depend on a library that we defined in the top-level file
target_link_libraries(HTTAnalysis RootAnalysisCore boost ${ROOT_LIBRARIES})

# Define an executable
add_executable(httAnalysis HTTAnalysis.cc)

target_include_directories(httAnalysis PUBLIC
    ${CMAKE_SOURCE_DIR}/include   
    ${ROOTINCDIR})

target_link_libraries(httAnalysis HTTAnalysis RootAnalysisCore ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} Proof RooFitCore RooFit)


