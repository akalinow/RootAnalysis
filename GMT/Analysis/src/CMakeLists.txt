project(GMT_Analysis LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include ${CMAKE_SOURCE_DIR}/Core/include ${CMAKE_SOURCE_DIR}/GMT/DataFormats/include)
include_directories(${ROOTINCDIR})

add_library(GMTAnalysis SHARED GMTAnalyzer.cc GMTHistograms.cc)

# Define an executable
add_executable(gmtAnalysis GMTAnalysis.cc)

target_link_libraries(gmtAnalysis GMTAnalysis ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} FrameworkCore GMTDataFormats)

install(TARGETS GMTAnalysis gmtAnalysis 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
install(FILES $<1:${PROJECT_SOURCE_DIR}/../config/> DESTINATION config)
install(DIRECTORY DESTINATION fig_png)
install(DIRECTORY DESTINATION fig_eps)
