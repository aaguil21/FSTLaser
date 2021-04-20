ROOT_GENERATE_DICTIONARY(G__FstQAStudy
  ${CMAKE_CURRENT_SOURCE_DIR}/FstQAStudy.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCons.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstEvent.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstRawHit.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCluster.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstTrack.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/FstQAStudyLinkDef.h)

add_library(FstQAStudy SHARED FstQAStudy.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstEvent.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstRawHit.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCluster.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstTrack.cxx G__FstQAStudy.cxx)
target_link_libraries(FstQAStudy ${ROOT_LIBRARIES})

install(TARGETS FstQAStudy DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstQAStudy_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstQAStudy.rootmap" DESTINATION lib)

