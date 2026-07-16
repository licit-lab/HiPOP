#include <hipop/graph.h>
#include <hipop/shortest_path.h>

int main()
{
    hipop::OrientedGraph G;

    /*G.AddNode("0", 0, 0);
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
    };*/
    
    //auto paths = hipop::parallelKShortestPath(G, origins, destinations, "time", vecMapLabelCosts, {}, 0, 10, 4, 4);

    /*
    G.AddNode("0", 0, 0);
    G.AddNode("1", 1, 0);
    G.AddNode("2", 1, 1);
    G.AddNode("3", 0, 1);

    G.AddLink("0_1", "0", "1", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("1_2", "1", "2", 1, {{"PersonalVehicle", {{"time", 13}}}}, "CAR");
    G.AddLink("0_3", "0", "3", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("3_2", "3", "2", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    */

    //auto path = hipop::dijkstra(G, "0", "2", "time", {{"CAR", "PersonalVehicle"}});

    G.AddNode("0", 0, 0);
    G.AddNode("1", 1, 0);
    G.AddNode("2", 1, 1);
    G.AddNode("3", 0, 1);

    G.AddLink("0_1", "0", "1", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("1_2", "1", "2", 1, {{"PersonalVehicle", {{"time", 13}}}}, "CAR");
    G.AddLink("0_3", "0", "3", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");
    G.AddLink("3_2", "3", "2", 1, {{"PersonalVehicle", {{"time", 12}}}}, "CAR");

    std::vector<std::string> origins = {"0", "1", "2", "3"}; 
    std::vector<std::unordered_map<std::string, std:: string> > vecMapLabelCosts = {
        {{"CAR", "PersonalVehicle"}}, 
        {{"CAR", "PersonalVehicle"}}, 
        {{"CAR", "PersonalVehicle"}}, 
        {{"CAR", "PersonalVehicle"}}
    };

    auto res1 = hipop::multiDestDijkstra(G, "0", "time", {{"CAR", "PersonalVehicle"}});
    auto res2 = hipop::parallelMultiDestDijkstra(G, origins, "time", vecMapLabelCosts, 4);

    return EXIT_SUCCESS;
}
