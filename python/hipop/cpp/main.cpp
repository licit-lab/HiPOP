#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <hipop/graph.h>
#include <hipop/shortest_path.h>
#include <hipop/create.h>

namespace py = pybind11;

namespace hipop_wrappers {
    void graph(py::module_ &);
    void shortest_path(py::module_ &);
}

PYBIND11_MODULE(cpp, m) {

    py::module graph = m.def_submodule("graph", "Graph module");
    hipop_wrappers::graph(graph);

    py::module shortest_path = m.def_submodule("shortest_path", "Shortest path module");
    hipop_wrappers::shortest_path(shortest_path);

}
