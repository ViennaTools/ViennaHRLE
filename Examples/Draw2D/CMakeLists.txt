cmake_minimum_required(VERSION 3.4)

project("Draw2D")

find_package(ViennaHRLE REQUIRED)

add_executable(${PROJECT_NAME} Draw2D.cpp)
target_link_libraries(${PROJECT_NAME} ViennaHRLE)
