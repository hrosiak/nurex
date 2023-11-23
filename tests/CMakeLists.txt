include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(NUREX_TESTS 
    #test_functions 
    test_evaporation
    #test_density 
    #test_nucleus 
    #test_nnsigma 
    #test_integration 
    test_gm 
    #test_nucid 
    #test_defaults
    #test_json
    test_profile_function
    #test_barrier
    )

foreach(entry ${NUREX_TESTS})
    add_executable(${entry} ${entry}.cpp)
    target_link_libraries(${entry} nurex)
    add_test(${entry} ${entry})
endforeach(entry in ${NUREX_TESTS})
