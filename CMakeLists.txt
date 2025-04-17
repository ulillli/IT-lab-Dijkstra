cmake_minimum_required(VERSION 3.10)

option(BUILD_SAMPLES ON)

set(PROJECT_NAME Dijkstra)
project(${PROJECT_NAME})

set(MP2_CUSTOM_PROJECT "${PROJECT_NAME}")
set(MP2_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(include)
add_subdirectory(main)


# REPORT
message( STATUS "")
message( STATUS "General configuration for ${PROJECT_NAME}")
message( STATUS "======================================")
message( STATUS "")
message( STATUS "   Configuration: ${CMAKE_BUILD_TYPE}")
message( STATUS "")
