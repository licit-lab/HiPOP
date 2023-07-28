from hipop.graph import generate_manhattan
from hipop.shortest_path import dijkstra, parallel_dijkstra

from time import time

g = generate_manhattan(100, 10)

s=dict()
s['']='PersonalCar'

print(dijkstra(g, "NORTH_0", "EAST_0", "length",s))

N = int(3000)

origins = ["NORTH_0"]*N
dests = ["EAST_0"]*N

s=dict()
s['']='PersonalCar'
services = [s for _ in range(N)]

print("Launch")
start = time()
res = parallel_dijkstra(g, origins, dests, services, "length", 8)
end = time()
print("Done", f"[{end-start} s]")

print(res[0])
