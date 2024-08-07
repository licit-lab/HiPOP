#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <hipop/shortest_path.h>

namespace py = pybind11;

namespace hipop_wrappers {

void shortest_path(py::module_ &m) {
    m.def(
        "dijkstra", 
        &hipop::dijkstra, 
        py::arg("graph"), 
        py::arg("origin"), 
        py::arg("destination"), 
        py::arg("cost"),
        py::arg("map_label_cost"),
        py::arg("available_labels") = setstring());
    m.def(
        "dijkstra_single_source",
        &hipop::dijkstraSingleSource,
        py::arg("graph"),
        py::arg("origin"),
        py::arg("cost"),
        py::arg("map_label_cost"),
        py::arg("available_labels") = setstring());
    m.def(
      "floyd_warshall",
      &hipop::floydWarshall,
      py::arg("graph"),
      py::arg("cost"),
      py::arg("map_label_cost"),
      py::arg("available_labels") = setstring()
    );
    m.def(
        "parallel_dijkstra", 
        &hipop::parallelDijkstra, 
        py::arg("graph"), 
        py::arg("origins"), 
        py::arg("destinations"),
        py::arg("map_label_costs"), 
        py::arg("cost"), 
        py::arg("thread_number"), 
        py::arg("available_labels") = std::vector<setstring>());
    m.def(
        "parallel_dijkstra_single_source",
        &hipop::parallelDijkstraSingleSource,
        py::arg("graph"),
        py::arg("origins"),
        py::arg("map_label_costs"),
        py::arg("cost"),
        py::arg("thread_number"),
        py::arg("available_labels") = std::vector<setstring>());
    m.def(
        "parallel_dijkstra_heterogeneous_costs",
        &hipop::parallelDijkstraHeterogeneousCosts,
        py::arg("graph"),
        py::arg("origins"),
        py::arg("destinations"),
        py::arg("map_label_costs"),
        py::arg("costs"), 
        py::arg("thread_number"),
        py::arg("available_labels") = std::vector<setstring>());
    m.def("k_shortest_path", &hipop::KShortestPath);
    m.def("parallel_k_shortest_path", &hipop::parallelKShortestPath);
    m.def("yen_k_shortest_path", &hipop::YenKShortestPath);
    m.def("astar_euclidian_dist", &hipop::aStarEuclidianDist);
    m.def("compute_path_length", &hipop::computePathLength);
    m.def("compute_path_cost", &hipop::computePathCost);
    m.def("compute_paths_costs", &hipop::computePathsCosts);
    m.def(
        "parallel_k_intermodal_shortest_path",
        &hipop::parallelKIntermodalShortestPath,
        py::arg("graph"),
        py::arg("origins"),
        py::arg("destinations"),
        py::arg("map_label_costs"),
        py::arg("cost"),
        py::arg("thread_number"),
        py::arg("pair_mandatory_labels"),
        py::arg("max_diff_cost"),
        py::arg("max_dist_in_common"),
        py::arg("cost_multiplier"),
        py::arg("max_retry"),
        py::arg("nb_paths"),
        py::arg("available_labels") = std::vector<setstring>());
}

}
