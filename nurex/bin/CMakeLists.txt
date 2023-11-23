set(NUREX_APPS nurex_calculator)

foreach(entry ${NUREX_APPS})
    add_executable(${entry} ${entry}.cpp)
    target_include_directories(${entry} PRIVATE
                            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libs>) 
    target_link_libraries(${entry} nurex)
endforeach(entry in ${NUREX_APPS})

install (TARGETS ${NUREX_APPS} RUNTIME DESTINATION bin)