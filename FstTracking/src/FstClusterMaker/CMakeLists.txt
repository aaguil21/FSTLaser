ROOT_GENERATE_DICTIONARY(G__FstClusterMaker
  ${CMAKE_CURRENT_SOURCE_DIR}/FstClusterMaker.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCons.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstEvent.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstRawHit.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCluster.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstTrack.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/FstClusterMakerLinkDef.h)

add_library(FstClusterMaker SHARED FstClusterMaker.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstEvent.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstRawHit.cxx ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCluster.cxx  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstTrack.cxx G__FstClusterMaker.cxx)
target_link_libraries(FstClusterMaker ${ROOT_LIBRARIES})

install(TARGETS FstClusterMaker DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstClusterMaker_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstClusterMaker.rootmap" DESTINATION lib)

