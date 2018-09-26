project(HTauTau_Analysis LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include ${CMAKE_SOURCE_DIR}/Core/include ${CMAKE_SOURCE_DIR}/HTauTau/DataFormats/include)
include_directories(${ROOTINCDIR})

add_library(Analysis SHARED ChannelSpecifics.cc HTTAnalyzer.cc
  HTTAnalyzerCPmethods.cc HTTHistograms.cc HTTSynchNTuple.cc
  MuTauSpecifics.cc TauTauSpecifics.cc Tools.cc BTagCalibrationStandalone.cpp MLAnalyzer.cc)

# Define an executable
add_executable(httAnalysis HTTAnalysis.cc)

target_link_libraries(httAnalysis Analysis ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} Proof RooFitCore RooFit FrameworkCore DataFormats)

install(TARGETS Analysis httAnalysis 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
install(FILES $<1:${PROJECT_SOURCE_DIR}/../config/> DESTINATION config)
install(FILES $<1:${PROJECT_SOURCE_DIR}/../python/> DESTINATION python)
install(DIRECTORY DESTINATION fig_png)
