
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  test
  test.cc
)

target_link_libraries(
  test
  GTest::gtest_main
  GascoigneStd
)


target_include_directories(test PUBLIC
	${GASCOIGNE_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}
)

configure_file(square.inp ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(square_2.inp ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Results)

include(GoogleTest)
gtest_discover_tests(test)
