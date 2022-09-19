#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <hipop/graph.h>
#include <hipop/create.h>

namespace py = pybind11;

namespace hipop_wrappers {

void graph(py::module_ &m) {
    py::class_<Link>(m, "Link")
        .def(py::init<std::string, std::string, std::string, double, std::unordered_map<std::string, double>,  std::string>(),
             py::arg("id"), py::arg("up"), py::arg("down"), py::arg("length"), py::arg("cost"), py::arg("label") = "")
        .def_readonly("id", &Link::mid)
        .def_readonly("upstream", &Link::mupstream)
        .def_readonly("downstream", &Link::mdownstream)
        .def_readonly("costs", &Link::mcosts)
        .def_readonly("label", &Link::mlabel)
        .def_readonly("length", &Link::mlength)
        .def("update_costs", &Link::updateCosts);

    py::class_<Node>(m, "Node")
          .def(py::init<std::string, double, double, std::string, std::unordered_map<std::string, std::set<std::string> > >(),
               py::arg("id"), py::arg("x"), py::arg("y"), py::arg("label") = "",py::arg("exclude_movements") = mapsets())
          .def_readonly("id", &Node::mid)
          .def_readonly("position", &Node::mposition)
          .def_readonly("adj", &Node::madj)
          .def_readonly("radj", &Node::mradj)
          .def_readonly("label", &Node::mlabel)
          .def_readonly("exclude_movements", &Node::mexclude_movements)
          .def("get_exits", &Node::getExits, py::arg("predecessor")="_default")
          .def("get_entrances", &Node::getEntrances, py::arg("predecessor")="_default");

    py::class_<OrientedGraph>(m, "OrientedGraph")
          .def(py::init<>())
          .def_readwrite("nodes", &OrientedGraph::mnodes)
          .def_readwrite("links", &OrientedGraph::mlinks)
          .def("add_node", py::overload_cast<std::string, double, double, std::string, mapsets>(&OrientedGraph::AddNode),
               py::arg("id"), py::arg("x"), py::arg("y"), py::arg("label"), py::arg("exclude_movements") = mapsets())
          .def("add_link", py::overload_cast<std::string, std::string, std::string, double, std::unordered_map<std::string, double>, std::string>(&OrientedGraph::AddLink),
               py::arg("id"), py::arg("up"), py::arg("down"), py::arg("length"), py::arg("costs"), py::arg("label") = "_def")
          .def("get_link", &OrientedGraph::getLink)
          .def("update_link_costs", &OrientedGraph::UpdateLinkCosts);

    m.def("generate_manhattan", &makeManhattan);

    m.def("merge_oriented_graph", &mergeOrientedGraph);
}

}
