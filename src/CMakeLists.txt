cmake_minimum_required(VERSION 3.16)
project(BarnesHut VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer")
string(APPEND CMAKE_EXE_LINKER_DEBUG " -fsanitize=address -fno-omit-frame-pointer")

find_package(SFML 2.5 COMPONENTS graphics REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(barnes-hut.out main.cpp initial_distr.cpp nbody-sim.cpp gravity.cpp barnes-hut.cpp)
target_link_libraries(barnes-hut.out PRIVATE sfml-graphics)
target_link_libraries(barnes-hut.out PUBLIC OpenMP::OpenMP_CXX)

