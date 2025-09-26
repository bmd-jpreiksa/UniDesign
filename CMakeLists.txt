cmake_minimum_required(VERSION 3.16)
project(UniDesign LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB UNIDESIGN_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
list(REMOVE_ITEM UNIDESIGN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Main.cpp")
list(APPEND UNIDESIGN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/support/ProgramGlobals.cpp")

add_library(unidesign ${UNIDESIGN_SOURCES})

target_include_directories(unidesign
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(unidesign PUBLIC cxx_std_17)

install(TARGETS unidesign
  EXPORT UniDesignTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

install(DIRECTORY include/ DESTINATION include)

install(EXPORT UniDesignTargets
  FILE UniDesignTargets.cmake
  NAMESPACE UniDesign::
  DESTINATION lib/cmake/UniDesign
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/UniDesignConfigVersion.cmake"
  VERSION 1.0.0
  COMPATIBILITY AnyNewerVersion
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/UniDesignConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/UniDesignConfig.cmake"
  @ONLY
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/UniDesignConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/UniDesignConfigVersion.cmake"
  DESTINATION lib/cmake/UniDesign
)

option(UNIDESIGN_BUILD_TESTS "Build UniDesign tests" ON)
if(UNIDESIGN_BUILD_TESTS)
  enable_testing()
  add_executable(StructureHandleTest tests/StructureHandleTest.cpp)
  target_link_libraries(StructureHandleTest PRIVATE unidesign)
  add_test(NAME StructureHandleTest COMMAND StructureHandleTest)
endif()
