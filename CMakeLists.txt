cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(
  ViennaHRLE
  LANGUAGES CXX
  VERSION 0.7.0)

# --------------------------------------------------------------------------------------------------------
# Library options
# --------------------------------------------------------------------------------------------------------

option(VIENNAHRLE_BUILD_TESTS "Build tests" OFF)
option(VIENNAHRLE_BUILD_EXAMPLES "Build examples" OFF)

option(VIENNAHRLE_USE_SANITIZER "Enable Sanitizers" OFF)

# --------------------------------------------------------------------------------------------------------
# Setup Library
# --------------------------------------------------------------------------------------------------------

add_library(${PROJECT_NAME} INTERFACE)
add_library(ViennaTools::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES CXX_STANDARD 17
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
  NAME ViennaCore
  VERSION 1.4.3
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaCore")

CPMAddPackage(
  NAME PackageProject
  VERSION 1.13.0
  GIT_REPOSITORY "https://github.com/TheLartians/PackageProject.cmake")

target_link_libraries(${PROJECT_NAME} INTERFACE ViennaTools::ViennaCore)

# --------------------------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------------------------

if(VIENNAHRLE_USE_SANITIZER)
  message(STATUS "[ViennaHRLE] Using Sanitizer")
  viennacore_enable_sanitizer()
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Examples
# --------------------------------------------------------------------------------------------------------

if(VIENNAHRLE_BUILD_EXAMPLES)
  message(STATUS "[ViennaHRLE] Building Examples")
  add_subdirectory(examples)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Tests
# --------------------------------------------------------------------------------------------------------

if(VIENNAHRLE_BUILD_TESTS)
  message(STATUS "[ViennaHRLE] Building Tests")

  enable_testing()
  add_subdirectory(tests)
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
