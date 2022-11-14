cmake_minimum_required(VERSION 3.1...3.16)

set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)

add_executable(srl_projection_test "TestEquirectangularCamera.cpp" "TestPinholeCamera.cpp" "TestRangeImageLidar.cpp" "TestLidar.cpp")
target_include_directories(srl_projection_test PRIVATE)
target_link_libraries(srl_projection_test
  PRIVATE
    SRL::Projection
    gtest_main
    ${CMAKE_THREAD_LIBS_INIT}
)
set_target_properties(${LIB_NAME} PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)
target_compile_options(srl_projection_test PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
	-Wall -pedantic -Wextra>
)
add_test(NAME srl_projection_test COMMAND srl_projection_test)

