cmake_minimum_required(VERSION 3.4...3.18)
project(pyhipop)

find_package(pybind11 REQUIRED)
find_package(OpenMP REQUIRED)

pybind11_add_module(cpp 
    hipop/cpp/main.cpp 
    hipop/cpp/graph.cpp
    hipop/cpp/shortest_path.cpp)
set_target_properties(cpp PROPERTIES CXX_STANDARD 17)

find_package(Hipop REQUIRED)
target_link_libraries(cpp PUBLIC hipoplib)