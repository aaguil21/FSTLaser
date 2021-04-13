ROOT_GENERATE_DICTIONARY(G__FstNoiseStudy
  ${CMAKE_CURRENT_SOURCE_DIR}/FstNoiseStudy.h
  ${CMAKE_SOURCE_DIR}/src/FstUtil/FstCons.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/FstNoiseStudyLinkDef.h)

add_library(FstNoiseStudy SHARED FstNoiseStudy.cxx G__FstNoiseStudy.cxx)
target_link_libraries(FstNoiseStudy ${ROOT_LIBRARIES})

install(TARGETS FstNoiseStudy DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstNoiseStudy_rdict.pcm" DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libFstNoiseStudy.rootmap" DESTINATION lib)

