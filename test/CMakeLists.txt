add_executable(parallel_test parallel_test.cpp)
target_include_directories(parallel_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(parallel_test libsiga Boost::unit_test_framework)
add_test(NAME parallel_test COMMAND parallel_test)

add_executable(preprocess_test preprocess_test.cpp)
target_include_directories(preprocess_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(preprocess_test libsiga Boost::unit_test_framework)
add_test(NAME preprocess_test COMMAND preprocess_test)

add_executable(index_test index_test.cpp)
target_include_directories(index_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(index_test libsiga Boost::unit_test_framework)
add_test(NAME index_test COMMAND preprocess_test)

add_executable(overlap_test overlap_test.cpp)
target_include_directories(overlap_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(overlap_test libsiga Boost::unit_test_framework)
add_test(NAME overlap_test COMMAND preprocess_test)

add_executable(assemble_test assemble_test.cpp)
target_include_directories(assemble_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(assemble_test libsiga Boost::unit_test_framework)
add_test(NAME assemble_test COMMAND preprocess_test)

add_executable(utils_test utils_test.cpp)
target_include_directories(utils_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(utils_test libsiga Boost::unit_test_framework)
add_test(NAME utils_test COMMAND preprocess_test)
