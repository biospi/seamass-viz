cmake_minimum_required(VERSION 2.8.12)
project(seamass-viz C CXX)

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
endif(CMAKE_COMPILER_IS_GNUCC)
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif(CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory(viz)
