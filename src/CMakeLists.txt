set(MPPCONVERTER_SRC converter.cpp)

set(GMP_LIB_DIR $ENV{HOME}/.local/lib/libgmp.a)
set(HEXL_HINT_DIR $ENV{HOME}/.local/lib/cmake/hexl-1.2.5)

find_package(HEXL 1.2.5
    HINTS ${HEXL_HINT_DIR}
    REQUIRED)
if (NOT TARGET HEXL::hexl)
    message(FATAL_ERROR "TARGET HEXL::hexl not found")
endif()

add_library(mppconverter STATIC ${MPPCONVERTER_SRC})
target_link_libraries(mppconverter PRIVATE ${GMP_LIB_DIR})