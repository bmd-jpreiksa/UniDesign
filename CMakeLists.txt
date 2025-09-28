cmake_minimum_required(VERSION 3.16)
project(UniDesign LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(unidesign STATIC
  src/Atom.cpp
  src/AtomParamsSet.cpp
  src/Chain.cpp
  src/DesignSite.cpp
  src/EnergyFunction.cpp
  src/EnergyMatrix.cpp
  src/ErrorTracker.cpp
  src/EvoAminoName.cpp
  src/EvoGetPhiPsi.cpp
  src/EvoGetSeq2SA.cpp
  src/EvoGetSeq2SS.cpp
  src/EvoSeqAlign.cpp
  src/EvoUtility.cpp
  src/Evolution.cpp
  src/GeometryCalc.cpp
  src/Getopt.cpp
  src/ProgramFunction.cpp
  src/ProgramPreprocess.cpp
  src/ProteinDesign.cpp
  src/Residue.cpp
  src/ResidueTopology.cpp
  src/Rotamer.cpp
  src/RotamerBuilder.cpp
  src/RotamerOptimizer.cpp
  src/Sequence.cpp
  src/SmallMol.cpp
  src/SmallMolEEF1.cpp
  src/SmallMolParAndTopo.cpp
  src/Structure.cpp
  src/Utility.cpp
  src/WeightOpt.cpp
)

target_sources(unidesign
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/unidesign/SequenceHandle.hpp>
)

target_include_directories(unidesign
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

target_compile_definitions(unidesign
  PUBLIC
    $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
)

add_executable(unidesign_cli src/Main.cpp)
target_link_libraries(unidesign_cli PRIVATE unidesign)

if (PROJECT_IS_TOP_LEVEL)
  enable_testing()

  add_executable(sequence_handle_tests tests/SequenceHandleTests.cpp)
  target_link_libraries(sequence_handle_tests PRIVATE unidesign)
  add_test(NAME SequenceHandleTests COMMAND sequence_handle_tests)

  find_package(pybind11 QUIET)
  if (pybind11_FOUND)
    pybind11_add_module(unidesign_python MODULE bindings/python/SequenceBindings.cpp)
    target_link_libraries(unidesign_python PRIVATE unidesign)
    target_include_directories(unidesign_python PRIVATE include src)
    target_compile_definitions(unidesign_python PRIVATE UNIDESIGN_ENABLE_PYBIND11)
  endif()
endif()

