# RawHit
ROOT_GENERATE_DICTIONARY(G__FstRawHit
  ${CMAKE_CURRENT_SOURCE_DIR}/FstRawHit.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FstCons.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/FstRawHitLinkDef.h)

add_library(FstRawHit SHARED FstRawHit.cxx G__FstRawHit.cxx)
target_link_libraries(FstRawHit ${ROOT_LIBRARIES})

install(TARGETS FstRawHit DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstRawHit_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstRawHit.rootmap" DESTINATION lib)

# Cluster
ROOT_GENERATE_DICTIONARY(G__FstCluster
  ${CMAKE_CURRENT_SOURCE_DIR}/FstCluster.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FstCons.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FstRawHit.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/FstClusterLinkDef.h)

add_library(FstCluster SHARED FstCluster.cxx ${CMAKE_CURRENT_SOURCE_DIR}/FstRawHit.cxx G__FstCluster.cxx)
target_link_libraries(FstCluster ${ROOT_LIBRARIES})

install(TARGETS FstCluster DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstCluster_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstCluster.rootmap" DESTINATION lib)

# Track
ROOT_GENERATE_DICTIONARY(G__FstTrack
  ${CMAKE_CURRENT_SOURCE_DIR}/FstTrack.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FstCons.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/FstTrackLinkDef.h)

add_library(FstTrack SHARED FstTrack.cxx G__FstTrack.cxx)
target_link_libraries(FstTrack ${ROOT_LIBRARIES})

install(TARGETS FstTrack DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstTrack_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstTrack.rootmap" DESTINATION lib)

# Event
ROOT_GENERATE_DICTIONARY(G__FstEvent
  ${CMAKE_CURRENT_SOURCE_DIR}/FstEvent.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FstCons.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FstRawHit.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FstCluster.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FstTrack.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/FstEventLinkDef.h)

add_library(FstEvent SHARED FstEvent.cxx ${CMAKE_CURRENT_SOURCE_DIR}/FstRawHit.cxx ${CMAKE_CURRENT_SOURCE_DIR}/FstCluster.cxx ${CMAKE_CURRENT_SOURCE_DIR}/FstTrack.cxx G__FstEvent.cxx)
target_link_libraries(FstEvent ${ROOT_LIBRARIES})

install(TARGETS FstEvent DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstEvent_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstEvent.rootmap" DESTINATION lib)

