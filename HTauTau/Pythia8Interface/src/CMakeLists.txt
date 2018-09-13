project(HTauTau_Pythia8Interface LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include ${CMAKE_SOURCE_DIR}/Core/include ${CMAKE_SOURCE_DIR}/HTauTau/DataFormats/include)
include_directories($ENV{PYTHIA8_DIR}/include)
include_directories(${ROOTINCDIR})

if (NOT EXISTS $ENV{PYTHIA8_DIR}/include/Pythia8/Pythia.h)
  MESSAGE(WARNING "Pythia8 not present. Pythia8Interface will not be compiled.")
  return()
endif()

add_library(Pythia8Interface SHARED Pythia8Interface.cc)

add_library(Pythia8 SHARED IMPORTED)
set_property(TARGET Pythia8 PROPERTY IMPORTED_LOCATION $ENV{PYTHIA8_DIR}/lib/libpythia8.so)
  
# Depend on a library that we defined in the top-level file
target_link_libraries(Pythia8Interface Core DataFormats ${ROOT_LIBRARIES})

# Define an executable
add_executable(generateEvents generateEvents.cc)

target_link_libraries(generateEvents Pythia8Interface Core DataFormats ${Boost_REGEX_LIBRARY} ${ROOT_LIBRARIES} ${PYTHIA8_LIBRARY} Proof RooFitCore RooFit EG EGPythia8 Pythia8)

set(CMAKE_INSTALL_RPATH $ENV{PYTHIA8_DIR}/lib)

install(TARGETS generateEvents Pythia8Interface
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
install(FILES $<1:${PROJECT_SOURCE_DIR}/../config/> DESTINATION config)
