from typing import Dict

from hipop.cpp.graph import Link, Node, OrientedGraph, generate_manhattan, merge_oriented_graph, copy_graph


def node_to_dict(node: Node):
    d = {"ID": node.id,
         "X": node.position[0],
         "Y": node.position[1],
         "LABEL": node.label,
         "EXCLUDE_MOVEMENTS": node.exclude_movements}
    return d


def dict_to_node(G:OrientedGraph, d:Dict):
    G.add_node(d["ID"],
               d["X"],
               d["Y"],
               d["LABEL"],
               d["EXCLUDE_MOVEMENTS"])


def link_to_dict(link:Link):
    d = {"ID": link.id,
         "UPSTREAM": link.upstream,
         "DOWNSTREAM": link.downstream,
         "LENGTH": link.length,
         "COSTS": link.costs,
         "LABEL": link.label}
    return d


def dict_to_link(G: OrientedGraph, d: Dict):
    G.add_link(d["ID"],
               d["UPSTREAM"],
               d["DOWNSTREAM"],
               d["LENGTH"],
               d["COSTS"],
               d["LABEL"])


def graph_to_dict(G: OrientedGraph) -> Dict:
    nodes = []
    links = []

    for n in G.nodes.values():
        d = node_to_dict(n)
        nodes.append(d)

    for l in G.links.values():
        d = link_to_dict(l)
        links.append(d)

    return {"NODES": nodes, "LINKS": links}


def dict_to_graph(d: Dict) -> OrientedGraph:
    G = OrientedGraph()

    for n in d["NODES"]:
        dict_to_node(G, n)

    for l in d["LINKS"]:
        dict_to_link(G, l)

    return G