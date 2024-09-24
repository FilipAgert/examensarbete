# Minimum required version of CMake
cmake_minimum_required(VERSION 3.16)

# Project name and version
project(examensarbete VERSION 0.1.0 LANGUAGES Fortran)

# Metadata
set(PROJECT_LICENSE "license")
set(PROJECT_AUTHOR "Filip Agert")
set(PROJECT_MAINTAINER "filip.agert@gmail.com")
set(PROJECT_COPYRIGHT "Copyright 2024, Filip Agert")

set(FVERSION "-std=f95")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")
# Set Fortran Standard (assuming Fortran 2008 as fpm defaults to this)
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

# Compiler options based on fpm.toml
add_compile_options(-fimplicit-none)  # For implicit-typing = false

# Define project sources (You may want to include subdirectories here)

# Define the main executable source
set(MAIN_SOURCE app/main.f90)

# Collect all other source files in the src directory
file(GLOB_RECURSE PROJECT_SOURCES src/*.f90)

# Combine main source and src sources
set(ALL_SOURCES ${MAIN_SOURCE} ${PROJECT_SOURCES})

# Add the main executable
add_executable(examensarbete ${ALL_SOURCES})


# Add external libraries (stdlib, fsparse)
# These should be cloned manually or through some script
# Use ExternalProject to download dependencies (simpler if projects provide CMake)

# Path to dependencies (stdlib and fsparse directories)
# Replace with your actual relative or absolute path
# Include the stdlib project
set(STDLIB_DIR ${CMAKE_SOURCE_DIR}/build/dependencies/stdlib/src/)
set(FSPARSE_DIR ${CMAKE_SOURCE_DIR}/build/dependencies/fspars/src/)
set(FSPARSE_CONV_DIR ${CMAKE_SOURCE_DIR}/build/dependencies/fspars/src/conversion)
set(FSPARSE_MATVEC__DIR ${CMAKE_SOURCE_DIR}/build/dependencies/fspars/src/matvec/)
set(FSPARSE_DIR ${CMAKE_SOURCE_DIR}/build/dependencies/fspars/src/)


# Collect stdlib source files
file(GLOB_RECURSE STDLIB_SOURCES ${STDLIB_DIR}/*.f90)

# Collect fsparse source files from main directory and specific subdirectories
file(GLOB_RECURSE FSPARSE_SOURCES ${FSPARSE_DIR}/*.f90)
file(GLOB_RECURSE FSPARSE_CONV_SOURCES ${FSPARSE_CONV_DIR}/*.f90)
file(GLOB_RECURSE FSPARSE_MATVEC_SOURCES ${FSPARSE_MATVEC_DIR}/*.f90)

# Combine all fsparse sources
set(ALL_FSPARSE_SOURCES ${FSPARSE_SOURCES} ${FSPARSE_CONV_SOURCES} ${FSPARSE_MATVEC_SOURCES})

# Add libraries for stdlib and fsparse
add_library(stdlib ${STDLIB_SOURCES})
add_library(fsparse ${ALL_FSPARSE_SOURCES})

# Link stdlib and fsparse to the main executable
target_link_libraries(examensarbete stdlib fsparse)

# Collect all test files in the tests directory
file(GLOB_RECURSE TEST_SOURCES tests/*.f90)

# Create an executable for each test (assuming one test per source file)
foreach(test_src ${TEST_SOURCES})
    # Extract the test name from the file name (e.g., test1.f90 -> test1)
    get_filename_component(test_name ${test_src} NAME_WE)

    # Add executable for the test
    add_executable(${test_name} ${test_src})

    # Link the test with dependencies (stdlib, fsparse)
    target_link_libraries(${test_name} stdlib fsparse)

    # Register the test with CTest
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

# Add tests (assuming auto-tests = true)
enable_testing()
add_test(NAME BasicTest COMMAND examensarbete)
