project(OMTF_Analysis LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include ${CMAKE_SOURCE_DIR}/Core/include ${CMAKE_SOURCE_DIR}/OMTF/DataFormats/include)
include_directories(${ROOTINCDIR})

add_library(OMTFAnalysis SHARED GMTAnalyzer.cc GMTHistograms.cc)

# Define an executable
add_executable(omtfAnalysis GMTAnalysis.cc)

target_link_libraries(omtfAnalysis OMTFAnalysis ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} FrameworkCore OMTFDataFormats)

install(TARGETS OMTFAnalysis omtfAnalysis 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
install(FILES $<1:${PROJECT_SOURCE_DIR}/../config/> DESTINATION config)
install(DIRECTORY DESTINATION fig_png)
install(DIRECTORY DESTINATION fig_eps)
