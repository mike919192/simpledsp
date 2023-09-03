cmake_minimum_required(VERSION 3.15)
project(testFFT)

find_package(Catch2 REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(${PROJECT_NAME} test/testFFT.cpp)
include_directories(include)
target_link_libraries(${PROJECT_NAME} Catch2::Catch2WithMain)
