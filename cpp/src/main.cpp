#include <iostream>
#include <hipop/graph.h>
#include <hipop/shortest_path.h>
//include "/Users/manon/Documents/mnms/HiPOP/cpp/include/hipop/graph.h"
//include "/Users/manon/Documents/mnms/HiPOP/cpp/include/hipop/shortest_path.h"

int main() {

    hipop::OrientedGraph G;
    
    G.AddNode("0", 0, 0);
    G.AddNode("1", 1, 0);
    G.AddNode("2", 1, 1);
    G.AddNode("3", 0, 1);

    G.AddLink("0_1", "0", "1", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("1_2", "1", "2", 1, {{"PersonalVehicle", {{"time", 13}}}}, "CAR");
    G.AddLink("0_3", "0", "3", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("3_2", "3", "2", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");

    /*
    auto res = hipop::multiDestDijkstra(G, "0", "time", {{"CAR", "PersonalVehicle"}});
    */

   std::cout << "Hello World!";
   return 0;
}

