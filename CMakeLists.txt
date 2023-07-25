cmake_minimum_required(VERSION 3.14)
project(MPPconverter)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_TEST "Build tests" ON)
option(ENABLE_EXAMPLE "Build examples" OFF)

set(MPPCONVERTER_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

if (ENABLE_TEST)
    enable_testing()
    add_subdirectory(test)
endif()

add_subdirectory(src)

if (ENABLE_EXAMPLE)
    add_subdirectory(example)
endif()