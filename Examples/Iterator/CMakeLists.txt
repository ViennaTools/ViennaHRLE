cmake_minimum_required(VERSION 3.4)

project("Iterator")

find_package(ViennaHRLE REQUIRED)

add_executable(${PROJECT_NAME} Iterator.cpp)
target_link_libraries(${PROJECT_NAME} ViennaHRLE)
