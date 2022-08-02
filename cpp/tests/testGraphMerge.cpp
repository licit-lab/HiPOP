#include "helpers.h"

#include <hipop/graph.h>

#include <unordered_map>
#include <set>
#include <string>
#include <memory>
#include <iostream>



int testGraphMerge(int argc, char *argv[])
{
    OrientedGraph* G1 = new OrientedGraph();
    G1->AddNode("a", 0, 0);
    G1->AddNode("b", 2, 5, "", {{"a", {"c"}}});
    G1->AddNode("c", 12, 43);
    G1->AddNode("d", 435, 345);

    G1->AddLink("a_b", "a", "b", 12, {{"time", 12}});
    G1->AddLink("b_c", "b", "c", 12, {{"time", 12}});
    G1->AddLink("b_d", "b", "d", 12, {{"time", 12}});


    OrientedGraph* G2 = new OrientedGraph();

    G2->AddNode("f", 39, 3);
    G2->AddNode("y", 42, 0);

    G2->AddLink("f_y", "f", "y", 22, {{"time", 22}});
    
    
    OrientedGraph* G3 = new OrientedGraph();

    G3->AddNode("h", 39, 3);

    OrientedGraph* mergeG = mergeOrientedGraph({G1, G2, G3});


    assertTrue(mergeG->mnodes.size()==7, "Merge graph does not have 7 nodes");
    assertTrue(mergeG->mlinks.size()==4, "Merge graph does not have 4 links");


    return 0;
}
