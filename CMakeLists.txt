cmake_minimum_required(VERSION 3.16)
project(UniDesign LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
file(GLOB LIB_SOURCES CONFIGURE_DEPENDS "${SRC_DIR}/*.cpp")
list(REMOVE_ITEM LIB_SOURCES "${SRC_DIR}/Main.cpp")

add_library(unidesign STATIC ${LIB_SOURCES})
set_target_properties(unidesign PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(unidesign
    PUBLIC
        ${SRC_DIR}
)

if (MSVC)
    target_compile_options(unidesign PRIVATE /W0 /O2)
else()
    target_compile_options(unidesign PRIVATE -w -O3 -ffast-math)
endif()

add_executable(UniDesign "${SRC_DIR}/Main.cpp")

target_link_libraries(UniDesign PRIVATE unidesign)

if (MSVC)
    target_compile_options(UniDesign PRIVATE /W0 /O2)
else()
    target_compile_options(UniDesign PRIVATE -w -O3 -ffast-math)
endif()
