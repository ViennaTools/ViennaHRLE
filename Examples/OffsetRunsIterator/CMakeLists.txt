cmake_minimum_required(VERSION 3.4)

project("OffsetRunsIterator")

find_package(ViennaHRLE REQUIRED)

add_executable(${PROJECT_NAME} OffsetRunsIterator.cpp)
target_link_libraries(${PROJECT_NAME} ViennaHRLE)
