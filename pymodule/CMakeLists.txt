find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(_pynurex pynurex.cpp)
target_link_libraries(_pynurex PRIVATE nurex)
file(COPY pynurex DESTINATION .)
add_custom_command(TARGET _pynurex POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE_NAME:_pynurex> ${CMAKE_CURRENT_BINARY_DIR}/pynurex/$<TARGET_FILE_NAME:_pynurex>
)
# Install python bindings
install(
TARGETS _pynurex
LIBRARY DESTINATION "pynurex"
RUNTIME DESTINATION "pynurex"
ARCHIVE DESTINATION "pynurex")
