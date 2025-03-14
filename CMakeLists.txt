cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(
  ViennaHRLE
  LANGUAGES CXX
  VERSION 0.6.0)

# --------------------------------------------------------------------------------------------------------
# Library options
# --------------------------------------------------------------------------------------------------------

option(VIENNAHRLE_OPENMP "Build OpenMP" OFF)

option(VIENNAHRLE_BUILD_TESTS "Build tests" OFF)
option(VIENNAHRLE_BUILD_EXAMPLES "Build examples" OFF)

option(VIENNAHRLE_ENABLE_SANITIZER "Enable Sanitizer for debug builds" OFF)
option(VIENNAHRLE_USE_SANITIZER "Enable Sanitizer (Requires GCC or Clang)" OFF)

# --------------------------------------------------------------------------------------------------------
# Ensure valid library options
# --------------------------------------------------------------------------------------------------------

if(VIENNAHRLE_BUILD_TESTS)
  set(VIENNAHRLE_OPENMP ON)
endif()

# --------------------------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------------------------

if((VIENNAHRLE_ENABLE_SANITIZER AND CMAKE_BUILD_TYPE STREQUAL "Debug") OR VIENNAHRLE_USE_SANITIZER)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=memory)
    add_link_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=memory)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=leak)
    add_link_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=leak)
  endif()

  message(STATUS "[ViennaHRLE] Using Sanitizer")
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Library
# --------------------------------------------------------------------------------------------------------

add_library(${PROJECT_NAME} INTERFACE)
add_library(ViennaTools::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)
set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES CXX_STANDARD 11
             CXX_EXTENSIONS OFF
             CXX_STANDARD_REQUIRED ON)

# --------------------------------------------------------------------------------------------------------
# Include directories
# --------------------------------------------------------------------------------------------------------

target_include_directories(
  ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/viennahrle>
                            $<INSTALL_INTERFACE:include/viennahrle-${PROJECT_VERSION}>)

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------

include("cmake/cpm.cmake")

CPMAddPackage(
  NAME PackageProject
  VERSION 1.13.0
  GIT_REPOSITORY "https://github.com/TheLartians/PackageProject.cmake")

if(VIENNAHRLE_OPENMP)
  message(STATUS "[ViennaHRLE] Enabling OMP Support")

  find_package(OpenMP REQUIRED)
  target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Examples
# --------------------------------------------------------------------------------------------------------

if(VIENNAHRLE_BUILD_EXAMPLES)
  message(STATUS "[ViennaHRLE] Building Examples")
  add_subdirectory(Examples)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Tests
# --------------------------------------------------------------------------------------------------------

if(VIENNAHRLE_BUILD_TESTS)
  message(STATUS "[ViennaHRLE] Building Tests")

  enable_testing()
  add_subdirectory(Tests)
endif()

# --------------------------------------------------------------------------------------------------------
# Install Target
# --------------------------------------------------------------------------------------------------------

packageProject(
  NAME ${PROJECT_NAME}
  VERSION ${PROJECT_VERSION}
  BINARY_DIR ${PROJECT_BINARY_DIR}
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/viennahrle
  INCLUDE_DESTINATION include/viennahrle-${PROJECT_VERSION})
