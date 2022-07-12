#include "helpers.h"

#include <hipop/graph.h>

#include <unordered_map>
#include <set>
#include <string>
#include <memory>
#include <iostream>



int testGraph(int argc, char *argv[])
{
    OrientedGraph G;
    G.AddNode("a", 0, 0);

    std::unordered_map<std::string, std::set<std::string> > excludeMovements;
    excludeMovements["a"] = {"c"};
    G.AddNode("b", 2, 5, "", excludeMovements);
    
    Node *newNode = new Node("c", 12., 43.);
    G.AddNode(newNode);
    
    G.AddNode("d", 435, 345);
    G.AddLink("a_b", "a", "b", 12, {{"time", 12}});
    G.AddLink("b_c", "b", "c", 12, {{"time", 12}});
    G.AddLink("b_d", "b", "d", 12, {{"time", 12}});

    std::vector<Link*> exits = G.mnodes["b"]->getExits("a");
    
    assertTrue(exits.size()==1, "Exits does not return one link");
    assertTrue(exits[0]->mdownstream=="d", "Node should be d");



    return 0;
}
