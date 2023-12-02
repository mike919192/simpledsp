cmake_minimum_required(VERSION 3.15)
project(testFFT)

file(GLOB test CONFIGURE_DEPENDS "test/*.cpp")

find_package(Catch2 REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

add_executable(${PROJECT_NAME} ${test})
add_compile_definitions(CATCH_CONFIG_MAIN)
include_directories(include)
target_link_libraries(${PROJECT_NAME} Catch2::Catch2WithMain)
