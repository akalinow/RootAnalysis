project(HSCP_Analysis LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include ${CMAKE_SOURCE_DIR}/Core/include ${CMAKE_SOURCE_DIR}/HSCP/DataFormats/include)
include_directories(${ROOTINCDIR})

add_library(HSCPAnalysis SHARED HSCPAnalyzer.cc HSCPHistograms.cc)

# Define an executable
add_executable(hscpAnalysis HSCPAnalysis.cc)

target_link_libraries(hscpAnalysis HSCPAnalysis ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} Proof RooFitCore RooFit FrameworkCore HSCPDataFormats)

install(TARGETS HSCPAnalysis hscpAnalysis 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
install(FILES $<1:${PROJECT_SOURCE_DIR}/../config/> DESTINATION config)
install(DIRECTORY DESTINATION fig_png)
