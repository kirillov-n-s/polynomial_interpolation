cmake_minimum_required(VERSION 3.20)
project(laba4)

set(CMAKE_CXX_STANDARD 20)

find_package(Matplot++ REQUIRED)

add_executable(laba4 main.cpp split.h interpol.h format.h)

target_link_libraries(laba4 PUBLIC Matplot++::matplot)
