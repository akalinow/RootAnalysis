project(OMTF_Analysis LANGUAGES CXX)

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/../include 
    ${CMAKE_SOURCE_DIR}/Core/include 
    ${CMAKE_SOURCE_DIR}/OMTF/DataFormats/include
    ${ROOTINCDIR}
)

# Add library
add_library(OMTFAnalysis SHARED OMTFAnalyzer.cc OMTFHistograms.cc)

# Add executable
add_executable(omtfAnalysis OMTFAnalysis.cc)

# Link libraries
target_link_libraries(
    omtfAnalysis 
    OMTFAnalysis 
    ${Boost_REGEX_LIBRARY} 
    ${ROOT_LIBRARIES} 
    FrameworkCore 
    OMTFDataFormats
)

# Install targets
install(
    TARGETS OMTFAnalysis omtfAnalysis 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
)

# Install files and directories
install(FILES $<1:${PROJECT_SOURCE_DIR}/../config/> DESTINATION config)
install(DIRECTORY DESTINATION fig_png)
install(DIRECTORY DESTINATION fig_eps)