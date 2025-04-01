add_custom_target(ViennaHRLE_Examples ALL)

file(
  GLOB entries
  LIST_DIRECTORIES true
  "*")

foreach(entry ${entries})
  if(NOT IS_DIRECTORY ${entry})
    continue()
  endif()

  message(STATUS "[ViennaHRLE] Adding Example ${entry}")
  add_subdirectory(${entry})
endforeach()
