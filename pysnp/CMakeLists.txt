cmake_minimum_required(VERSION 3.10)

project(pysnpprocess) 

set(CMAKE_CXX_STANDARD 17)
set(PYBIND11_PYTHON_VERSION 3.6)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC")

find_package(pybind11 REQUIRED PATHS /usr/local/lib/python3.6/site-packages/pybind11/)

pybind11_add_module(${PROJECT_NAME} pysnp.cpp)

target_compile_definitions(${PROJECT_NAME} PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_include_directories(${PROJECT_NAME} PRIVATE ${PYBIND11_INCLUDE_DIRS})
