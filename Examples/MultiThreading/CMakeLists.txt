cmake_minimum_required(VERSION 3.4)

project("MultiThreading")

find_package(ViennaHRLE REQUIRED)

add_executable(${PROJECT_NAME} MultiThreading.cpp)
target_link_libraries(${PROJECT_NAME} ViennaHRLE)
