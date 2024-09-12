cmake_minimum_required(VERSION 3.12)
project(montecarlo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable OpenMP
find_package(OpenMP REQUIRED)

# Add compiler flags for optimization and AVX2 support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mavx2")

# Check if AVX2 is supported
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)

if(COMPILER_SUPPORTS_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
else()
    message(WARNING "AVX2 is not supported by the compiler")
endif()

add_executable(montecarlo main.cpp)

# Link OpenMP
target_link_libraries(montecarlo PRIVATE OpenMP::OpenMP_CXX)

# Compiler version checks
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
        message(FATAL_ERROR "GCC version must be at least 10!")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
        message(FATAL_ERROR "Clang version must be at least 14!")
    endif()
endif()