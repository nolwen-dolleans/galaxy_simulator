cmake_minimum_required(VERSION 3.20)

# --------------------------------------------------
# Sources
# --------------------------------------------------
add_subdirectory(libsim)

add_executable(main
    main.cpp
)

target_compile_options(main PRIVATE
    $<$<CONFIG:Debug>:-Wall -Wextra -g -fsanitize=address>
    $<$<CONFIG:Release>:-Wall -Wextra -O3 -march=native -DNDEBUG>
)

target_link_libraries(main PRIVATE
    libsim
)
