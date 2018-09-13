project(HTauTau_SVfitAnalysis LANGUAGES CXX)

##Add SVfit library
set (ClassicSVfit_PATH ${CMAKE_SOURCE_DIR}/../)
include_directories(${ClassicSVfit_PATH})

include_directories(${PROJECT_SOURCE_DIR}/../include ${CMAKE_SOURCE_DIR}/Core/include)
include_directories(${CMAKE_SOURCE_DIR}/HTauTau/DataFormats/include)
include_directories(${CMAKE_SOURCE_DIR}/HTauTau/Analysis/include)
include_directories(${ROOTINCDIR})

add_library(SVfitTF SHARED ${ClassicSVfit_PATH}/TauAnalysis/SVfitTF/src/crystalBall.cc)

if (NOT EXISTS ${ClassicSVfit_PATH}/TauAnalysis/SVfitTF/src/crystalBall.cc)
  MESSAGE(WARNING "SVfitTF not present. SVfitAnalysis will not be compiled")
  return()
endif()

if (NOT EXISTS ${ClassicSVfit_PATH}/TauAnalysis/ClassicSVfit/src/ClassicSVfit.cc)
  MESSAGE(WARNING "SVfitTF not present. SVfitAnalysis will not be compiled")
  return()
endif() 

add_library(ClassicSVfit SHARED ${ClassicSVfit_PATH}/TauAnalysis/ClassicSVfit/src/ClassicSVfit.cc
  ${ClassicSVfit_PATH}/TauAnalysis/ClassicSVfit/src/ClassicSVfitIntegrand.cc
  ${ClassicSVfit_PATH}/TauAnalysis/ClassicSVfit/src/MeasuredTauLepton.cc
  ${ClassicSVfit_PATH}/TauAnalysis/ClassicSVfit/src/svFitAuxFunctions.cc
  ${ClassicSVfit_PATH}/TauAnalysis/ClassicSVfit/src/svFitHistogramAdapter.cc
  ${ClassicSVfit_PATH}/TauAnalysis/ClassicSVfit/src/SVfitIntegratorMarkovChain.cc
  ${ClassicSVfit_PATH}/TauAnalysis/ClassicSVfit/src/FastMTT.cc)

add_library(SVfitAnalysis SHARED SVfitAnalyzer.cc)

# Depend on a library that we defined in the top-level file
target_link_libraries(SVfitAnalysis Core DataFormats Analysis ClassicSVfit SVfitTF ${ROOT_LIBRARIES})

# Define an executable
add_executable(svfitAnalysis SVfitAnalysis.cc)

target_link_libraries(svfitAnalysis SVfitAnalysis Core DataFormats Analysis ClassicSVfit SVfitTF
  ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} Proof RooFitCore RooFit)

install(TARGETS SVfitAnalysis ClassicSVfit SVfitTF svfitAnalysis
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib)
install(FILES $<1:${PROJECT_SOURCE_DIR}/../config/> DESTINATION config)
