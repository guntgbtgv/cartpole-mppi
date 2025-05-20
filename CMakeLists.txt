cmake_minimum_required(VERSION 3.16)

# Make CMAKE_C_VISIBILITY_PRESET work properly.
set(CMAKE_POLICY_DEFAULT_CMP0063 NEW)
# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
# Default to GLVND if available.
set(CMAKE_POLICY_DEFAULT_CMP0072 NEW)
# Avoid BUILD_SHARED_LIBS getting overridden by an option() in ccd.
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

set(CMAKE_CXX_STANDARD 17)
# This line has to appear before 'PROJECT' in order to be able to disable incremental linking
set(MSVC_INCREMENTAL_DEFAULT ON)

project(mppi_test)

enable_language(C)
enable_language(CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/mujoco/simulate/cmake")
option(MUJOCO_BUILD_EXAMPLES "Build samples for MuJoCo" OFF)
option(MUJOCO_BUILD_TESTS "Build tests for MuJoCo" OFF)

add_subdirectory(${PROJECT_SOURCE_DIR}/mujoco) 

include(SimulateDependencies)

find_package(Eigen3 REQUIRED)

add_executable(${PROJECT_NAME} src/main.cc src/mppi_controller.cpp)


target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/mujoco
  ${EIGEN3_INCLUDE_DIR}
)


target_link_libraries(
  ${PROJECT_NAME}
  libsimulate
  mujoco::mujoco
  glfw
  Threads::Threads
  lodepng
)