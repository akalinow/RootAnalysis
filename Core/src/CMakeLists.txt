project(Core LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include)
include_directories(${ROOTINCDIR})

add_library(FrameworkCore SHARED EventProxyBase.cc AnalysisHistograms.cc Analyzer.cc
  ObjectMessenger.cc strbitset.cc TFileDirectory.cc TFileService.cc
  TH1AddDirectorySentry.cc TreeAnalyzer.cc SummaryAnalyzer.cc
  commonUtils.cc)

install(TARGETS FrameworkCore
        LIBRARY DESTINATION lib)
