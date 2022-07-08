cmake_minimum_required(VERSION 3.13)
project(linkphase3)
set(CMAKE_CXX_STANDARD 11)

FIND_PACKAGE( Boost 1.4 COMPONENTS program_options REQUIRED )

add_executable(linkphase src/linkphase.cpp src/AnimalInfo.cpp  src/HalfsibPhaser.cpp)

target_link_libraries(linkphase Boost::program_options)
