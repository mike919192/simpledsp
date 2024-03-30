cmake_minimum_required(VERSION 3.15)
project(testFFT)

find_package(Catch2 REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wconversion")
endif()

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE test/testFFT.cpp test/testIIR.cpp)
add_compile_definitions(CATCH_CONFIG_MAIN)
include_directories(include)
target_link_libraries(${PROJECT_NAME} PRIVATE Catch2::Catch2WithMain)
