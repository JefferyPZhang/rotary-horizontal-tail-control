# ChatGPT-Generated Makefile

cmake_minimum_required(VERSION 3.15)
project(Horizontal_Rotary_Tail LANGUAGES CXX)

# GLOBAL SETTINGS
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Prefer Release builds
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# DEPENDENCIES
set(EIGEN_DIR "${CMAKE_SOURCE_DIR}/extern/eigen-5.0.0")
include_directories(${EIGEN_DIR})

# SOURCE FILES
set(SOURCES
    src/main.cpp
    src/Aircraft.cpp
    src/Dynamics.cpp
    src/Controllers.cpp
    src/Utils.cpp
)

# EXECUTABLE
add_executable(Horizontal_Rotary_Tail ${SOURCES})

target_include_directories(Horizontal_Rotary_Tail PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

# OUTPUT DIRECTORY
set_target_properties(Horizontal_Rotary_Tail PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/bin"
)
