set(HIPOP_SRC
  create.cpp
  graph.cpp
  shortest_path.cpp
)

add_library(hipoplib ${HIPOP_SRC})
set_target_properties(hipoplib PROPERTIES POSITION_INDEPENDENT_CODE ON)  # Needed to link the hipop static lib to the pybind dynamic one

target_link_libraries(hipoplib PUBLIC OpenMP::OpenMP_CXX)

set_property(TARGET hipoplib PROPERTY CXX_STANDARD 17)

target_include_directories( hipoplib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:include>
)