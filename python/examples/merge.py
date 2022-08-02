from importlib.resources import path
from hipop.graph import OrientedGraph
from hipop.cpp import merge_oriented_graph

G = OrientedGraph()

G.add_node("0", 0, 0, 'N0')
G.add_node("1", 1, 0, 'N0')
G.add_node("2", 1, 1, 'N0')
G.add_node("3", 0, 1, 'N0')

G.add_link("0_1", "0", "1", 1, {"time": 10})
G.add_link("0_3", "0", "3", 1, {"time": 11})
G.add_link("0_2", "0", "2", 1.44, {"time": 10})
G.add_link("1_2", "1", "2", 4, {"time": 10})
G.add_link("3_2", "3", "2", 5, {"time": 9})
G.add_link("1_3", "1", "3", 1, {"time": 10})

mergedG = merge_oriented_graph([G])