cmake_minimum_required(VERSION 3.16)
project(m2424 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(SEAL_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(SEAL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(SEAL_USE_ZLIB OFF CACHE BOOL "" FORCE)
set(SEAL_USE_ZSTD OFF CACHE BOOL "" FORCE)
set(SEAL_USE_MSGSL OFF CACHE BOOL "" FORCE)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/seal/CMakeLists.txt")
    message(FATAL_ERROR "SEAL submodule not found. Run: git submodule update --init --recursive")
endif()
add_subdirectory(extern/seal)

add_subdirectory(src)
add_subdirectory(apps)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
