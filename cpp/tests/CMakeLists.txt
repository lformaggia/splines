include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.2
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

set(TEST_NAMES
  vec
  tridiagonal
  spline_curve
  arc_length
  spline_surface
  geometry
  nearest_point
)

foreach(t ${TEST_NAMES})
  add_executable(test_${t} test_${t}.cpp)
  target_link_libraries(test_${t} PRIVATE splines Catch2::Catch2WithMain)
  target_compile_options(test_${t} PRIVATE -O2)
  catch_discover_tests(test_${t})
endforeach()
