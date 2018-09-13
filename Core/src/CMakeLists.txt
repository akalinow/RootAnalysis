project(Core LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include)

add_library(Core SHARED EventProxyBase.cc AnalysisHistograms.cc Analyzer.cc
  ObjectMessenger.cc strbitset.cc TFileDirectory.cc TFileService.cc
  TH1AddDirectorySentry.cc TreeAnalyzer.cc SummaryAnalyzer.cc
  commonUtils.cc)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(Core PUBLIC ${ROOTINCDIR})

# Depend on a library that we defined in the top-level file
target_link_libraries(Core boost ${ROOT_LIBRARIES})

install(TARGETS Core
        LIBRARY DESTINATION lib)
