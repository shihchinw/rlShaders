file(GLOB CPP_FILES *.cpp ext/*.cc)
file(GLOB H_FILES *.h ext/*.h)

file(GLOB EXT_FILES ext/*.h ext/*.cc)
source_group("ext" FILES ${EXT_FILES})

set(SOURCE_FILES ${CPP_FILES} ${H_FILES})

include_directories(${ARNOLD_INCLUDE_DIRS})

add_library(rlShaders MODULE ${SOURCE_FILES})

set_target_properties(rlShaders
                      PROPERTIES
                      COMPILE_FLAGS ${ARNOLD_COMPILE_FLAGS}
                      LINK_FLAGS ${ARNOLD_LINK_FLAGS}
                      PREFIX "")

target_link_libraries(rlShaders ${ARNOLD_LIBRARIES})
install(TARGETS rlShaders DESTINATION shaders)
install(FILES rlShaders.mtd DESTINATION shaders)
