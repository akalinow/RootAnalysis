project(OMTF_DataFormats LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/../include )

# Lines below are without adding new objects

# ROOT_GENERATE_DICTIONARY(G__OMTFObjects EventObj.h GenObj.h GenObjColl.h L1Obj.h L1ObjColl.h  LINKDEF ${PROJECT_SOURCE_DIR}/../include/LinkDef.h)
# # 
# add_library(OMTFDataFormats SHARED EventObj.cc GenObj.cc L1Obj.cc L1ObjColl.cc G__OMTFObjects.cxx EventProxyOMTF.cc) 
# # 
# With new objects

ROOT_GENERATE_DICTIONARY(G__OMTFObjects EventObj.h GenObj.h GenObjColl.h L1Obj.h L1ObjColl.h L1PhaseIIObj.h L1PhaseIIObjColl.h LINKDEF ${PROJECT_SOURCE_DIR}/../include/LinkDef.h)
# 
add_library(OMTFDataFormats SHARED EventObj.cc GenObj.cc L1Obj.cc L1PhaseIIObj.cc L1ObjColl.cc L1PhaseIIObj.cc L1PhaseIIObjColl.cc  G__OMTFObjects.cxx EventProxyOMTF.cc)
# 
# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(OMTFDataFormats PUBLIC ${ROOTINCDIR}  ${CMAKE_SOURCE_DIR}/Core/include)
  
# Depend on a library that we defined in the top-level file
target_link_libraries(OMTFDataFormats)

install(TARGETS OMTFDataFormats LIBRARY DESTINATION lib)

if(EXISTS libOMTFObjects_rdict.pcm)
  install(FILES libOMTFObjects_rdict.pcm libOMTFObjects.rootmap DESTINATION lib)
endif()

 