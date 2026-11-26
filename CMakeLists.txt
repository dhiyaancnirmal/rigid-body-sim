cmake_minimum_required(VERSION 3.18)
project(rbsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(RBSIM_PYTHON "Build the Python module" ON)

# Determinism: no fused multiply-add contraction and no fast-math, so every
# build performs the same IEEE operations in the same order. No -march=native.
set(RBSIM_FP_FLAGS "-ffp-contract=off;-fno-fast-math" CACHE STRING "Floating-point flags for the core")

add_library(rbsim_core STATIC src/shape.cpp src/collide.cpp src/world.cpp)
target_include_directories(rbsim_core PUBLIC include)
target_compile_options(rbsim_core PUBLIC ${RBSIM_FP_FLAGS} -Wall -Wextra -Wpedantic)
set_target_properties(rbsim_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

enable_testing()
add_executable(rbsim_tests tests/cpp/test_main.cpp)
target_link_libraries(rbsim_tests PRIVATE rbsim_core)
add_test(NAME rbsim_tests COMMAND rbsim_tests)

add_executable(bench_step bench/bench_step.cpp)
target_link_libraries(bench_step PRIVATE rbsim_core)
add_executable(bench_integrators bench/bench_integrators.cpp)
target_include_directories(bench_integrators PRIVATE include)
target_compile_options(bench_integrators PRIVATE ${RBSIM_FP_FLAGS})

if(RBSIM_PYTHON)
  find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_rbsim python/bindings.cpp)
  target_link_libraries(_rbsim PRIVATE rbsim_core)
  set_target_properties(_rbsim PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/rbsim)
endif()
