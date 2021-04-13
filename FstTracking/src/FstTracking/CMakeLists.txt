ROOT_GENERATE_DICTIONARY(G__FstTracking
  ${CMAKE_CURRENT_SOURCE_DIR}/FstTracking.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCons.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstEvent.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstRawHit.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCluster.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstTrack.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/FstTrackingLinkDef.h)

add_library(FstTracking SHARED FstTracking.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstEvent.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstRawHit.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCluster.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstTrack.cxx G__FstTracking.cxx)
target_link_libraries(FstTracking ${ROOT_LIBRARIES})

install(TARGETS FstTracking DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstTracking_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstTracking.rootmap" DESTINATION lib)

