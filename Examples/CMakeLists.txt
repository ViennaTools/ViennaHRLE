cmake_minimum_required(VERSION 3.4)

project(ViennaHRLEExamples)

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      IF(NOT ${child} STREQUAL "build")
        LIST(APPEND dirlist ${child})
      ENDIF()
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

set(ViennaHRLE_DIR "${ViennaHRLE_BINARY_DIR}")
set(ViennaHRLE_COMMON_TARGET "ViennaHRLEExamples")

SUBDIRLIST(EXAMPLE_DIRS ${CMAKE_CURRENT_LIST_DIR})

# add all examples
foreach(subdir ${EXAMPLE_DIRS})
  add_subdirectory(${subdir})
endforeach()
