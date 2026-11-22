cmake_minimum_required(VERSION 3.16)
project(Azimuth CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (MSVC)
    add_compile_options(/O2 /arch:AVX2 /W4)
else()
    add_compile_options(-O3 -march=native -Wall -Wextra -Wpedantic)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(pathfinder main.cpp)

if (WIN32)
    target_link_libraries(pathfinder psapi)
endif()
