cmake_minimum_required (VERSION 3.8)
project ("primer design")

set(CMAKE_CXX_STANDARD 17)
set(PYBIND11_PYTHON_VERSION 3.6)
set(CMAKE_CXX_FLAGS "-w -fPIC")

find_package(pybind11 REQUIRED PATHS /usr/local/lib/python3.6/site-packages/pybind11/)

file(GLOB_RECURSE ALL_SRC "src/*.cpp" "src/*.h")

pybind11_add_module(primer_design ${ALL_SRC})
target_compile_definitions(primer_design PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_include_directories(primer_design PRIVATE ${PYBIND11_INCLUDE_DIRS})



