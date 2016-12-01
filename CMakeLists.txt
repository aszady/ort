cmake_minimum_required(VERSION 2.8.4)
project(example)

add_subdirectory(../../cppjson ${CMAKE_CURRENT_BINARY_DIR}/cppjson)
add_subdirectory(../../gogui_core ${CMAKE_CURRENT_BINARY_DIR}/gogui_core)

#set(CMAKE_CXX_COMPILER "clang++-3.8")
set(CMAKE_CXX_FLAGS "-std=c++1y -O2 -static")


include_directories(".")
include_directories("../../cppjson/include")
include_directories("../../gogui_core/include")


add_executable(ort main.cpp)
target_link_libraries(ort gogui_core cppjson)
