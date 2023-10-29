cmake_minimum_required(VERSION 3.15)
project(testFFT)

find_package(Catch2 REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

add_executable(${PROJECT_NAME} test/testFFT.cpp test/testIIR.cpp)
include_directories(include)
target_link_libraries(${PROJECT_NAME} Catch2::Catch2WithMain)
