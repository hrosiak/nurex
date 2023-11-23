find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(_pynurex pynurex.cpp)
target_link_libraries(_pynurex PRIVATE nurex)
# Install python bindings
install(
TARGETS _pynurex
LIBRARY DESTINATION "pynurex"
RUNTIME DESTINATION "pynurex"
ARCHIVE DESTINATION "pynurex")
