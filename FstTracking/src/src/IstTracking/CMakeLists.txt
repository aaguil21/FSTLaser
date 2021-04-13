ROOT_GENERATE_DICTIONARY(G__IstTracking
  ${CMAKE_CURRENT_SOURCE_DIR}/IstTracking.h
  ${CMAKE_CURRENT_SOURCE_DIR}/IstTrackingCons.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/IstTrackingLinkDef.h)

add_library(IstTracking SHARED IstTracking.cxx G__IstTracking.cxx)
target_link_libraries(IstTracking ${ROOT_LIBRARIES})

install(TARGETS IstTracking DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libIstTracking_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libIstTracking.rootmap" DESTINATION lib)

