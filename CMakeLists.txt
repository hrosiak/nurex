cmake_minimum_required(VERSION 3.14...3.25)
project(nurex VERSION "1.0")

############ options #############
option(BUILD_SHARED_LIBS "build as shared library" ON)
option(THREADS "Use multi-threading" ON)
option(PYTHON_MODULE "compile the Nurex python module" OFF)
option(TESTS "build tests" OFF)
option(EXAMPLES "build examples" OFF)
option(GSL_INTEGRATION "use GSL library for integrations" OFF)
option(GSL_INTERPOLATION "use GSL library for interpolation" OFF)
option(APPS "build nurex applications" OFF)
option(NO_FILESYSTEM "do not build file system functions" OFF)
option(EVAPORATION_BARRIER_PARAMETRIZED "default barrier is set to parametrize" OFF)
######## build type ############
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    #set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wfatal-errors -fno-omit-frame-pointer -fsanitize=address")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wfatal-errors -fsanitize=address")
endif()
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
MESSAGE(STATUS "install prefix: " ${CMAKE_INSTALL_PREFIX})
MESSAGE( STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )

############# Requirements ##################
if(GSL_INTEGRATION OR GSL_INTERPOLATION)
    find_package(GSL REQUIRED)
    MESSAGE(STATUS "nurex will use GSL library" ${GSL_INCLUDE_DIRS})
    MESSAGE(STATUS "GSL include dir: " ${GSL_INCLUDE_DIRS})
endif()

if(THREADS)
    find_package(Threads REQUIRED)
    set (USE_THREADS ON)
    MESSAGE(STATUS "Nurex will use threads")
endif(THREADS)

include(FetchContent)
FetchContent_Declare(nlohmann_json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(nlohmann_json)

######## create build config files  ###################
configure_file(
                    "${CMAKE_CURRENT_SOURCE_DIR}/nurex/Config.h.in"
                    "${CMAKE_CURRENT_BINARY_DIR}/include/nurex/Config.h"
                    )

#configure_file("${PROJECT_SOURCE_DIR}/init.sh.in"
#               "${PROJECT_BINARY_DIR}/init.sh"
#              )

######################## main build ##################################

if(PYTHON_MODULE)
    set(BUILD_SHARED_LIBS OFF)
    set(TESTS OFF)
    set(CMAKE_CXX_VISIBILITY_PRESET hidden)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
    MESSAGE(STATUS "building python wheel")
endif(PYTHON_MODULE)

file(GLOB SOURCES nurex/*.cpp)
file(GLOB HEADERS nurex/*.h libs/*.h)

set(LIB_SRC
    nurex/libs/spline.cpp
    )
if(GSL_INTEGRATION OR GSL_INTERPOLATION)
LIST( APPEND LIB_SRC libs/gsl_integrator.cpp)
endif()

add_library(nurex ${SOURCES} ${LIB_SRC})

set_target_properties(nurex PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(nurex PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(nurex PRIVATE cxx_std_17)
target_compile_features(nurex INTERFACE cxx_std_17)

target_include_directories(nurex PUBLIC
                                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nurex/libs>
                                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                                )
target_link_libraries(nurex PRIVATE $<BUILD_INTERFACE:nlohmann_json::nlohmann_json>)

if(USE_THREADS)
    target_link_libraries(nurex PUBLIC Threads::Threads)
endif(USE_THREADS)

if (GSL_INTEGRATION OR GSL_INTERPOLATION)
    target_link_libraries(nurex PUBLIC GSL::gsl)
endif ()


add_library(nurex::nurex ALIAS nurex)
########## Sub Directories ###########

if(APPS)
add_subdirectory("bin")
endif(APPS)

if(TESTS)
enable_testing()
add_subdirectory("tests")
endif(TESTS)

if(BENCHMARKS)
add_subdirectory("benchmarks")
endif(BENCHMARKS)

if(PYTHON_MODULE)
add_subdirectory("pymodule")
return()
endif(PYTHON_MODULE)


####### DOCS generation #########
if(DOCS)
find_package(Doxygen REQUIRED)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(docs
                   COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    COMMENT "Generating API documentation with Doxygen" VERBATIM )
endif(DOXYGEN_FOUND)
endif(DOCS)

####### install part #######
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/nurexConfigVersion.cmake
        VERSION ${PACKAGE_VERSION}
        COMPATIBILITY AnyNewerVersion
)

include(GNUInstallDirs)

# Generate targets 
install(TARGETS nurex
        EXPORT nurexTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


# Install header files (TODO, this can be replaced with FILE_SET in cmake 3.23)
install (FILES ${HEADERS} "${CMAKE_CURRENT_BINARY_DIR}/include/nurex/Config.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Generate CMake file into the install tree (rather than export(TARGETS) which installs into build tree)
# Directly generate `Config` instead of `Targets`, as we have no additional configuration
install(EXPORT nurexTargets 
        NAMESPACE nurex:: 
        DESTINATION ${CMAKE_INSTALL_DATADIR}/nurex
        FILE nurexConfig.cmake)
# Install version
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/nurexConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATADIR}/nurex
)

###### packaging #######
set(CPACK_PACKAGE_NAME "nurex")
set(CPACK_PACKAGE_VENDOR "ProchazkaTakechi")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Nurex library - nuclear reaction cross section calculator")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Nurex")
include(CPack)