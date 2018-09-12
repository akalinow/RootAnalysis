project(SVFitAnalysis LANGUAGES CXX)

##Add SVFit library
set (ClassicSVfit_PATH "/scratch_local/akalinow/CMS/SVFit/development/")
set (ClassicSVfitLib_PATH "/scratch_local/akalinow/CMS/SVFit/development/")

include_directories(${ClassicSVfit_PATH})
link_directories(${ClassicSVfitLib_PATH})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(G1__HTTEvent HTTEvent.h LINKDEF LinkDef.h)

add_library(SVFitAnalysis SHARED ChannelSpecifics.cc  EventProxyHTT.cc svfitAnalyzer.cc
  svfitAnalyzerUtilityMethods.cc svfitHistograms.cc
  MuTauSpecifics.cc TauTauSpecifics.cc  MuMuSpecifics.cc
  Tools.cc  HTTEvent.cxx MLAnalyzer.cc  PropertyEnum.cc
  G1__HTTEvent.cxx)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(SVFitAnalysis PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${ROOTINCDIR})

# Depend on a library that we defined in the top-level file
target_link_libraries(SVFitAnalysis RootAnalysisCore boost ${ROOT_LIBRARIES} ClassicSVfit)

# Define an executable
add_executable(svfitAnalysis svfitAnalysis.cc)

target_include_directories(svfitAnalysis PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${ROOTINCDIR})

target_link_libraries(svfitAnalysis SVFitAnalysis RootAnalysisCore ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} Proof RooFitCore RooFit)
