cmake_minimum_required(VERSION 4.1)
project(libsim)


# --------------------------------------------------
# LLVM + OpenMP sur macOS
# --------------------------------------------------
if(APPLE)
    execute_process(
        COMMAND brew --prefix llvm
        OUTPUT_VARIABLE LLVM_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(EXISTS "${LLVM_PREFIX}/bin/clang++")
        message(STATUS "Using LLVM from Homebrew: ${LLVM_PREFIX}")
        set(CMAKE_C_COMPILER   "${LLVM_PREFIX}/bin/clang"   CACHE FILEPATH "" FORCE)
        set(CMAKE_CXX_COMPILER "${LLVM_PREFIX}/bin/clang++" CACHE FILEPATH "" FORCE)
    else()
        message(FATAL_ERROR "LLVM not found. Install with: brew install llvm")
    endif()
endif()

project(galaxy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --------------------------------------------------
# OpenMP
# --------------------------------------------------
find_package(OpenMP REQUIRED)

set(SRCS
    src/simulator.cpp
    )
    
set(HEADERS
    include/simulator.hpp
    )
    

add_library(libsim ${SRCS} ${HEADERS})
target_include_directories(libsim PUBLIC include)
target_include_directories(libsim PUBLIC .)


find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test tests/test_runner.cpp)
find_package(GTest REQUIRED)

target_link_libraries(test PRIVATE libsim GTest::GTest GTest::Main)

enable_testing()
find_package(GTest REQUIRED)

add_test(test libsim)

gtest_discover_tests(test)

target_link_libraries(libsim PRIVATE
    OpenMP::OpenMP_CXX
)
