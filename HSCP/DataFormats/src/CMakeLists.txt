project(HSCP_DataFormats LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include )

add_library(HSCPDataFormats SHARED HSCPEvent.cc EventProxyHSCP.cc)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(HSCPDataFormats PUBLIC ${ROOTINCDIR}  ${CMAKE_SOURCE_DIR}/Core/include)
  
# Depend on a library that we defined in the top-level file
target_link_libraries(HSCPDataFormats)

install(TARGETS HSCPDataFormats LIBRARY DESTINATION lib)

