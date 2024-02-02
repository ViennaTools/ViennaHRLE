cmake_minimum_required(VERSION 3.14)
project(MultiThreading LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaHRLE)

add_dependencies(ViennaHRLE_Examples ${PROJECT_NAME})
