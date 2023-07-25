set(MPPCONVERTER_SRC converter.cpp)

set(GMP_LIB_DIR $ENV{HOME}/.local/lib/libgmp.a)

add_library(mppconverter STATIC ${MPPCONVERTER_SRC})
target_link_libraries(mppconverter PRIVATE ${GMP_LIB_DIR})