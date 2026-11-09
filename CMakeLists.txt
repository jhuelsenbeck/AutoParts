cmake_minimum_required(VERSION 3.10)

# On Apple silicon, build a native arm64 binary even when CMake itself is an Intel
# (x86_64) build running under Rosetta; otherwise the program is compiled for x86_64
# and runs translated. Must be set before project(). Override with -DCMAKE_OSX_ARCHITECTURES=...
if(APPLE AND NOT CMAKE_OSX_ARCHITECTURES)
    execute_process(COMMAND sysctl -n hw.optional.arm64
                    OUTPUT_VARIABLE APPLE_SILICON_HOST
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    if(APPLE_SILICON_HOST STREQUAL "1")
        set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "macOS build architecture" FORCE)
    endif()
endif()

project(Autoparts)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

find_package(Eigen3 QUIET)

if(NOT Eigen3_FOUND)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(eigen)
endif()

file(GLOB_RECURSE SOURCES
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

add_executable(autoparts ${SOURCES})

target_include_directories(autoparts PRIVATE
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(autoparts PRIVATE
    Eigen3::Eigen
)

# Tune for the build machine: -mcpu=native on ARM (including Apple silicon),
# -march=native on x86
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mcpu=native" COMPILER_SUPPORTS_MCPU_NATIVE)
check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MCPU_NATIVE)
    target_compile_options(autoparts PRIVATE -mcpu=native)
elseif(COMPILER_SUPPORTS_MARCH_NATIVE)
    target_compile_options(autoparts PRIVATE -march=native)
endif()

target_compile_options(autoparts PRIVATE
    -O3
)