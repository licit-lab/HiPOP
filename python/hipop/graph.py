from os import link
from typing import Dict

from hipop.cpp import Link, Node, OrientedGraph, generate_manhattan, merge_oriented_graph


def graph_to_dict(G: OrientedGraph) -> Dict:
    nodes = []
    links = []

    for n in G.nodes.values():
        d = {"ID": n.id,
             "X": n.position[0],
             "Y": n.position[1],
             "LABEL": n.label,
             "EXCLUDE_MOVEMENTS": n.exclude_movements}
        nodes.append(d)
    
    for l in G.links.values():
        d = {"ID": l.id,
             "UPSTREAM": l.upstream,
             "DOWNSTREAM": l.downstream,
             "LENGTH": l.length,
             "COSTS": l.costs,
             "LABEL": l.label}
        links.append(d)
    
    return {"NODES": nodes, "LINKS":links}


def dict_to_graph(d: Dict) -> OrientedGraph:
    G = OrientedGraph()

    for n in d["NODES"]:
        G.add_node(n["ID"],
                   n["X"],
                   n["Y"],
                   n["LABEL"],
                   n["EXCLUDE_MOVEMENTS"])
    
    for l in d["LINKS"]:
        G.add_link(l["ID"],
                   l["UPSTREAM"],
                   l["DOWNSTREAM"],
                   l["LENGTH"],
                   l["COSTS"],
                   l["LABEL"])

    return G