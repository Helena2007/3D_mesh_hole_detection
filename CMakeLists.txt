cmake_minimum_required(VERSION 3.10)
project(HoleDetection)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PCL REQUIRED COMPONENTS
    common
    io
    surface
    search
    kdtree
)

add_executable(HoleDetection
    main.cpp
)

target_include_directories(HoleDetection PRIVATE
    ${PCL_INCLUDE_DIRS}
)

target_compile_definitions(HoleDetection PRIVATE
    ${PCL_DEFINITIONS}
)

target_link_libraries(HoleDetection PRIVATE
    ${PCL_LIBRARIES}
)

