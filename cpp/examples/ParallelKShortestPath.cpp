#include <hipop/graph.h>
#include <hipop/shortest_path.h>

int main()
{
    hipop::OrientedGraph G;

    G.AddNode("0", 0, 0);
    G.AddNode("1", 1, 1);
    G.AddNode("2", 1, -1);
    G.AddNode("3", 2, 0);
    G.AddNode("4", 2, 1);

    G.AddLink("0_1", "0", "1", 1, {{"PersonalVehicle", {{"time", 14}}}}, "CAR");
    G.AddLink("1_3", "1", "3", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("0_2", "0", "2", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("2_3", "2", "3", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("0_3", "0", "3", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("0_4", "0", "4", 11, {{"PersonalVehicle", {{"time", 14}}}}, "CAR");
    G.AddLink("4_3", "4", "3", 11, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");



    std::vector<std::string> origins = {"0", "0", "0", "0"}; 
    std::vector<std::string> destinations = {"3", "3", "3", "3"};
    std::vector<std::unordered_map<std::string, std:: string> > vecMapLabelCosts = {
        {{"CAR", "PersonalVehicle"}}, 
        {{"CAR", "PersonalVehicle"}}, 
        {{"CAR", "PersonalVehicle"}}, 
        {{"CAR", "PersonalVehicle"}}
    };
    
    //auto paths = hipop::parallelKShortestPath(G, origins, destinations, "time", vecMapLabelCosts, {}, 0, 10, 4, 4);
    auto paths = hipop::parallelDijkstra(G, origins, destinations, vecMapLabelCosts, "time", 4, {});

    return EXIT_SUCCESS;
}