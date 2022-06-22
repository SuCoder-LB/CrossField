cmake_minimum_required(VERSION 3.21)
project(CrossField)

set(CMAKE_CXX_STANDARD 14)

include_directories(libs/eigen)
include_directories(.)

set(CMAKE_EXE_LINKER_FLAGS "-static")

add_executable(CrossField main.cc cross_field.cc linear_solver.cc)
