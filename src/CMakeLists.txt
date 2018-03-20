project(RootAnalysisCore LANGUAGES CXX)

add_library(RootAnalysisCore SHARED EventProxyBase.cc AnalysisHistograms.cc Analyzer.cc
  ObjectMessenger.cc strbitset.cc TFileDirectory.cc TFileService.cc
  TH1AddDirectorySentry.cc TreeAnalyzer.cc SummaryAnalyzer.cc
  commonUtils.cc)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(RootAnalysisCore PUBLIC
    ${CMAKE_SOURCE_DIR}/include   
    ${ROOTINCDIR})

# Depend on a library that we defined in the top-level file
target_link_libraries(RootAnalysisCore boost ${ROOT_LIBRARIES})

