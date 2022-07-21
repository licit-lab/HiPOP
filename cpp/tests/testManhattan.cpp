#include "helpers.h"

#include <hipop/graph.h>
#include <hipop/create.h>

#include <unordered_map>
#include <set>
#include <string>
#include <memory>
#include <iostream>
#include <unistd.h>

int testManhattan(int argc, char *argv[])
{
    OrientedGraph *G = makeManhattan(500, 10);
    return 0;
}
