cmake_minimum_required(VERSION 3.12)
project(montecarlo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(montecarlo main.cpp)

# If using GCC, ensure it's a recent version
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
        message(FATAL_ERROR "GCC version must be at least 10!")
    endif()
endif()

# If using Clang, ensure it's a recent version
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
        message(FATAL_ERROR "Clang version must be at least 14!")
    endif()
endif()