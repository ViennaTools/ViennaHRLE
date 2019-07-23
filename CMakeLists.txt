cmake_minimum_required(VERSION 3.4)

project(
  "ViennaHRLE"
  VERSION 0.0.1)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)

# set the correct paths for installation
set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}-${PROJECT_VERSION}/")
set(LOCAL_CONFIG_DIR "lib/cmake/${PROJECT_NAME}-${PROJECT_VERSION}")

# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUInstallDirs)
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

# Run tests
#enable_testing()
#add_subdirectory(test)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                               COMPATIBILITY AnyNewerVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
  ${LOCAL_CONFIG_DIR})

# build examples
# Do we want examples built?
option(BUILD_EXAMPLES "Build examples." OFF)
if(BUILD_EXAMPLES)
  add_subdirectory(Examples)
endif()

#################################################
# INSTALL
#################################################
# install config files
# locations are provided by GNUInstallDirs
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${LOCAL_CONFIG_DIR})


install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${LOCAL_CONFIG_DIR})

# install include files
file(GLOB_RECURSE HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*")
install(FILES ${HEADER_FILES} DESTINATION include)
