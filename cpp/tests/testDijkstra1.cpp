#include "helpers.h"

#include <hipop/graph.h>
#include <hipop/shortest_path.h>


int testDijkstra1(int argc, char *argv[])
{   
    hipop::OrientedGraph G;

    G.AddNode("0", 0, 0);
    G.AddNode("1", 1, 0);
    G.AddNode("2", 1, 1);
    G.AddNode("3", 0, 1);

    G.AddLink("0_1", "0", "1", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("1_2", "1", "2", 1, {{"PersonalVehicle", {{"time", 13}}}}, "CAR");
    G.AddLink("0_3", "0", "3", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("3_2", "3", "2", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");

    auto path = hipop::dijkstra(G, "0", "2", "time", {{"CAR", "PersonalVehicle"}});

    assertTrue(path.second==24, "Path cost not equal to 24");
    assertTrue(path.first==std::vector<std::string>{"0", "3", "2"}, "Path nodes not equal to 0, 3, 2");

    return 0;
}

