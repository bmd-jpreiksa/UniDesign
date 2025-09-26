cmake_minimum_required(VERSION 3.16)
project(UniDesign LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
file(GLOB LIB_SOURCES CONFIGURE_DEPENDS "${SRC_DIR}/*.cpp")
list(REMOVE_ITEM LIB_SOURCES "${SRC_DIR}/Main.cpp")

find_package(pybind11 CONFIG QUIET)
if (NOT pybind11_FOUND)
    message(STATUS "pybind11 not found via CMake package, attempting to use bundled subdirectory")
    set(PYBIND11_SUBMODULE "${CMAKE_CURRENT_SOURCE_DIR}/extbin/pybind11")
    if (EXISTS "${PYBIND11_SUBMODULE}/CMakeLists.txt")
        add_subdirectory("${PYBIND11_SUBMODULE}" EXCLUDE_FROM_ALL)
    else()
        message(FATAL_ERROR "pybind11 could not be located. Install pybind11 or add it as a submodule in extbin/pybind11")
    endif()
endif()

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

set(PYBIND11_TARGET_NAME "unidesign_python")
pybind11_add_module(${PYBIND11_TARGET_NAME} MODULE
    src/python/bindings.cpp
)

set_target_properties(${PYBIND11_TARGET_NAME} PROPERTIES OUTPUT_NAME "_core")
target_link_libraries(${PYBIND11_TARGET_NAME} PRIVATE unidesign)
target_include_directories(${PYBIND11_TARGET_NAME}
    PRIVATE
        ${SRC_DIR}
)

if (MSVC)
    target_compile_options(${PYBIND11_TARGET_NAME} PRIVATE /W0 /O2)
else()
    target_compile_options(${PYBIND11_TARGET_NAME} PRIVATE -w -O3 -ffast-math)
endif()
