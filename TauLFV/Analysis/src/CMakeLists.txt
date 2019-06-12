include_directories(${PROJECT_SOURCE_DIR}/Core/include )
include_directories(${PROJECT_SOURCE_DIR}/TauLFV/DataFormats/include )
include_directories(${PROJECT_SOURCE_DIR}/TauLFV/Analysis/include )

include_directories(${ROOTINCDIR})

add_library(Analysis SHARED TauLFVAnalyzer.cc TauLFVHistograms.cc Tools.cc)

# Define an executable
add_executable(lfvAnalysis TauLFVAnalysis.cc)

target_link_libraries(lfvAnalysis Analysis ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} Proof RooFitCore RooFit FrameworkCore DataFormats)

install(TARGETS Analysis lfvAnalysis 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
install(FILES $<1:${PROJECT_SOURCE_DIR}/TauLFV/Analysis/config/> DESTINATION config)
install(DIRECTORY DESTINATION fig_png)
