set(RNSTFHE_SRC converter.cpp)

set(LOCAL_LIB_DIR /home/jukia/.local/lib)
set(HEXL_HINT_DIR /home/jukia/.local/lib/cmake/hexl-1.2.5)

find_package(HEXL 1.2.5
    HINTS ${HEXL_HINT_DIR}
    REQUIRED)
if (NOT TARGET HEXL::hexl)
    message(FATAL_ERROR "TARGET HEXL::hexl not found")
endif()

add_library(rnstfhe ${RNSTFHE_SRC})
target_link_libraries(rnstfhe PRIVATE ${LOCAL_LIB_DIR}/libgmp.a)