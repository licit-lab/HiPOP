#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <hipop/graph.h>
#include <hipop/shortest_path.h>
#include <hipop/create.h>

namespace py = pybind11;


void graph(py::module_ &);
void shortest_path(py::module_ &);

PYBIND11_MODULE(cpp, m) {
    
    graph(m);
    shortest_path(m);

}
