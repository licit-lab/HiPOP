from importlib.resources import path
from hipop.graph import OrientedGraph
from hipop.shortest_path import parallel_k_shortest_path
from pprint import pprint

G = OrientedGraph()

G.add_node("0", 0, 0., "ed")
G.add_node("1", 1, 0., "ed")
G.add_node("2", 1, 1., "ed")
G.add_node("3", 0, 1., "ed")

G.add_link("0_1", "0", "1", 1, {'PV':{'travel_time':15}})
G.add_link("0_3", "0", "3", 1,{'PV':{'travel_time':12}})
G.add_link("0_2", "0", "2", 1.44,{'PV':{'travel_time':10}})
G.add_link("1_2", "1", "2", 4,{'PV':{'travel_time':11}})
G.add_link("3_2", "3", "2", 5, {'PV':{'travel_time':20}})
G.add_link("1_3", "1", "3", 1, {'PV':{'travel_time':8}})

N = int(1e6)
origins = ["0" for _ in range(N)]
destinations = ["2" for _ in range(N)]
labels = [set() for _ in range(N)]

paths = parallel_k_shortest_path(G, origins, destinations, "time", labels, 0, 100, 3, 8)
pprint(paths[-30:])
                                                     