#include "helpers.h"

#include <hipop/graph.h>
#include <hipop/shortest_path.h>



int testparallelKShortestPath(int argc, char *argv[])
{
    OrientedGraph G;

    G.AddNode("0", 0, 0);
    G.AddNode("1", 1, 1);
    G.AddNode("2", 1, -1);
    G.AddNode("3", 2, 0);
    G.AddNode("4", 2, 1);

    G.AddLink("0_1", "0", "1", 1, {{"time", 14}});
    G.AddLink("1_3", "1", "3", 1, {{"time", 12}});
    G.AddLink("0_2", "0", "2", 1, {{"time", 12}});
    G.AddLink("2_3", "2", "3", 1, {{"time", 12}});
    G.AddLink("0_3", "0", "3", 1, {{"time", 12}});
    G.AddLink("0_4", "0", "4", 11, {{"time", 14}});
    G.AddLink("4_3", "4", "3", 11, {{"time", 12}});



    std::vector<std::string> origins = {"0", "0", "0", "0"}; 
    std::vector<std::string> destinations = {"3", "3", "3", "3"}; 

    auto paths = parallelKShortestPath(G, origins, destinations, "time", {}, 0, 10, 4, 4);


    return 0;
}
