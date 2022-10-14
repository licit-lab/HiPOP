#include "helpers.h"

#include <hipop/graph.h>
#include <hipop/create.h>


int testManhattan(int argc, char *argv[])
{
    OrientedGraph *G = makeManhattan(100, 10);
    return 0;
}
