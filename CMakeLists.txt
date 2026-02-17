cmake_minimum_required(VERSION 3.20)
project(galaxy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --------------------------------------------------
# Sources
# --------------------------------------------------
add_subdirectory(libsim)

add_executable(main
    main.cpp
)

target_compile_options(main PRIVATE
    $<$<CONFIG:Debug>:-Wall -Wextra -g>
    $<$<CONFIG:Release>:-Wall -Wextra -O3 -march=native -DNDEBUG>
)

target_link_libraries(main PRIVATE
    libsim
)
