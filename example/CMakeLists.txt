set(HEXL_HINT_DIR $ENV{HOME}/.local/lib/cmake/hexl-1.2.5)

find_package(HEXL 1.2.5
    HINTS ${HEXL_HINT_DIR}
    REQUIRED)
if (NOT TARGET HEXL::hexl)
    message(FATAL_ERROR "TARGET HEXL::hexl not found")
endif()

add_executable(speed_test speed_test.cpp)

target_include_directories(speed_test PRIVATE ${MPPCONVERTER_SRC_DIR})

target_link_libraries(speed_test PRIVATE mppconverter HEXL::hexl)