#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <hipop/shortest_path.h>

namespace py = pybind11;

namespace hipop_wrappers {

void shortest_path(py::module_ &m) {
    m.def(
        "dijkstra", 
        &dijkstra, 
        py::arg("graph"), 
        py::arg("origin"), 
        py::arg("destination"), 
        py::arg("cost"),
        py::arg("map_label_cost"),
        py::arg("available_labels") = setstring());
    m.def(
        "parallel_dijkstra", 
        &parallelDijkstra, 
        py::arg("graph"), 
        py::arg("origins"), 
        py::arg("destinations"),
        py::arg("map_label_costs"), 
        py::arg("cost"), 
        py::arg("thread_number"), 
        py::arg("available_labels") = std::vector<setstring>());
    m.def("k_shortest_path", &KShortestPath);
    m.def("parallel_k_shortest_path", &parallelKShortestPath);
    m.def("yen_k_shortest_path", &YenKShortestPath);
    m.def("astar_euclidian_dist", &aStarEuclidianDist);
    m.def("compute_path_length", &computePathLength);
}

}
