cmake_minimum_required(VERSION 3.4)

project("CrossIterator")

find_package(ViennaHRLE REQUIRED)

add_executable(${PROJECT_NAME} CrossIterator.cpp)
target_link_libraries(${PROJECT_NAME} ViennaHRLE)
