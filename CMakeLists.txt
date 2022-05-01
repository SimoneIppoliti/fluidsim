cmake_minimum_required(VERSION 3.12)
project(FluidSimulation)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(lib/piksel)

set(SOURCES
    src/main.cpp
    src/tools.cpp
    src/tools.hpp
    src/fluid.cpp
    src/fluid.hpp
    src/app.cpp
)

add_executable(fluidsim ${SOURCES})

target_include_directories(fluidsim
    PUBLIC src
    PUBLIC lib/piksel/src
)

target_link_libraries(fluidsim
    piksel
)