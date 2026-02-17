cmake_minimum_required(VERSION 3.16)
project(libsim LANGUAGES CXX)

# -------------------
# Library
# -------------------
add_library(libsim
    src/simulator.cpp
    include/simulator.hpp
)

target_include_directories(libsim
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# -------------------
# OpenMP
# -------------------
if(APPLE)
    # Apple clang + Homebrew libomp
    target_compile_options(libsim PUBLIC -Xclang -fopenmp -march=native)
    target_include_directories(libsim PUBLIC /opt/homebrew/opt/libomp/include)
    target_link_libraries(libsim PUBLIC /opt/homebrew/opt/libomp/lib/libomp.dylib)
else()
    find_package(OpenMP REQUIRED)
    target_link_libraries(libsim PUBLIC OpenMP::OpenMP_CXX)
endif()

# -------------------
# Tests
# -------------------
find_package(GTest REQUIRED)

add_executable(test tests/test_runner.cpp)
target_link_libraries(test
    PRIVATE
        libsim
        GTest::gtest
        GTest::gtest_main
)

if(APPLE)
    target_compile_options(test PRIVATE -Xclang -fopenmp)
    target_include_directories(test PRIVATE /opt/homebrew/opt/libomp/include)
    target_link_libraries(test PRIVATE /opt/homebrew/opt/libomp/lib/libomp.dylib)
else()
    target_link_libraries(test PRIVATE OpenMP::OpenMP_CXX)
endif()

enable_testing()
gtest_discover_tests(test)
