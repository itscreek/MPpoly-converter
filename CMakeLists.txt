cmake_minimum_required(VERSION 3.14)
project(MPPconverter)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_TEST "Build tests" OFF)

set(MPPCONVERTER_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

if (ENABLE_TEST)
    enable_testing()
    add_subdirectory(test)
endif()

add_subdirectory(src)
add_subdirectory(example)