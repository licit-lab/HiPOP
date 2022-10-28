#include "helpers.h"

#include <hipop/graph.h>


int testCopyGraph(int argc, char *argv[])
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

    hipop::OrientedGraph *copyG = hipop::copyGraph(G);

    G.mlinks["0_1"]->mcosts["PersonalVehicle"]["time"] = 42;

    assertTrue(
        G.mlinks["0_1"]->mcosts["PersonalVehicle"]["time"] != copyG->mlinks["0_1"]->mcosts["PersonalVehicle"]["time"],
        "Cost of original and copy should be different");

    assertTrue(
        G.mlinks["0_3"]->mcosts["PersonalVehicle"]["time"] == copyG->mlinks["0_3"]->mcosts["PersonalVehicle"]["time"],
        "Cost of original and copy should be the same");

    assertTrue(
        G.mlinks["0_3"] != copyG->mlinks["0_3"],
        "Pointers on links should be different");

    return 0;
}