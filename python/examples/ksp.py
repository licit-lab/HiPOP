from importlib.resources import path
from hipop.graph import OrientedGraph
from hipop.shortest_path import parallel_k_shortest_path, compute_path_length 
from pprint import pprint

G = OrientedGraph()

G.add_node("0", 0, 0.,'pv_layer')
G.add_node("1", 1, 0.,'pv_layer')
G.add_node("2", 1, 1.,'pv_layer')
G.add_node("3", 0, 1.,'pv_layer')

G.add_link("0_1", "0", "1", 1, {'PV':{'travel_time':15}},'pv_layer')
G.add_link("0_3", "0", "3", 1,{'PV':{'travel_time':12}},'pv_layer')
G.add_link("0_2", "0", "2", 1.44,{'PV':{'travel_time':100}},'pv_layer')
G.add_link("1_2", "1", "2", 4,{'M':{'travel_time':11}},'metro')
G.add_link("3_2", "3", "2", 5, {'PV':{'travel_time':20}},'pv_layer')
G.add_link("1_3", "1", "3", 1, {'PV':{'travel_time':8}},'pv_layer')

N = int(1e6)
N=3
origins = ["0" for _ in range(N)]
destinations = ["2" for _ in range(N)]

# Usable layer
layers = [{'pv_layer'} for _ in range(N)]

# Usable services (for each layer)
s=dict()
s['pv_layer']='PV'
services = [s for _ in range(N)]

# Number of k shortest paths to calculate
kpath=2

# Number of threads to use
nthread=8

# Minimal distance difference between the k-rd shortest path and the shortest path
min_dist=0

# Maximal distance difference between the k-rd shortest path and the shortest path
max_dist=5

paths = parallel_k_shortest_path(G, origins, destinations, 'travel_time',services, layers,min_dist,max_dist,kpath,nthread)
pprint(paths)

print('Lengths: ',compute_path_length(G,paths[0][0][0]),compute_path_length(G,paths[0][1][0]))

                                                     