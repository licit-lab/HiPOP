#include "hipop/graph.h"


#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <algorithm>
#include <queue>
#include <numeric>
#include <limits>


OrientedGraph::~OrientedGraph() {
    // std::cout<<"Del Graph "<<this<<std::endl;
    // std::cout<<"Links Graph "<<&mlinks<<std::endl;
    for (auto iter : mlinks) {
        // std::cout<<" Del Link "<<iter.first<<" "<<iter.second<<std::endl;
        // printf("Del Link %p %s \n", iter.second, iter.first.c_str());
        delete iter.second;
        iter.second = NULL;
    }
    
    for (auto iter : mnodes) {
        // for(auto &keyVal:iter.second->madj) {
        //     keyVal.second = NULL;
        // }

        // for(auto &keyVal:iter.second->mradj) {
        //     keyVal.second = NULL;
        // }
        // std::cout<<" Del Node "<<iter.first<<" "<<iter.second<<std::endl;
        delete iter.second;
        iter.second = NULL;
    }

    mnodes.clear();
    mlinks.clear();
}


void OrientedGraph::AddNode(std::string _id, double x, double y, std::string label, mapsets excludeMovements) {
    Node *new_node = new Node(_id, x, y, label, excludeMovements);
    mnodes[_id] = new_node;
};

void OrientedGraph::AddNode(Node* n) {
    mnodes[n->mid] = n;
};


void OrientedGraph::AddLink(std::string _id, std::string _up, std::string _down, double length, mapcosts _costs, std::string label) {
    Link *new_link = new Link(_id, _up, _down, length, _costs, label);
    mnodes[_up]->madj.emplace(_down, new_link);
    mnodes[_down]->mradj.emplace(_up, new_link);

    mlinks.emplace(_id, new_link);
};


void OrientedGraph::AddLink(Link *l) {
    mnodes[l->mupstream]->madj[l->mdownstream] = l;
    mnodes[l->mdownstream]->mradj[l->mupstream] = l;

    mlinks[l->mid] = l;
};


void OrientedGraph::UpdateLinkCosts(std::string lid, mapcosts _costs) {
    mlinks[lid]->updateCosts(_costs);
}


void OrientedGraph::ShowNodes() {
    for(const auto &elem: mnodes) {
        std::cout << "Node(" << elem.first << ", [" << elem.second->mposition[0] << ",\t" << elem.second->mposition[1] << "])\n";
    }
}

void OrientedGraph::ShowLinks() {
    for(const auto &elem: mlinks) {
        std::cout << "Link(" << elem.first << ", " << elem.second->mupstream << ", " << elem.second->mdownstream << ")\n";
    }
}



OrientedGraph* copyGraph(const OrientedGraph &G) {
    OrientedGraph* newGraph = new OrientedGraph();

    for(const auto &keyVal: G.mnodes) {
        // Node *n = new Node();

        newGraph->AddNode(keyVal.second->mid,
                           keyVal.second->mposition[0],
                           keyVal.second->mposition[1],
                           keyVal.second->mlabel,
                           keyVal.second->mexclude_movements);
    }


    for(const auto &keyVal: G.mlinks) {
        // Link *l = new Link(*keyVal.second);
        newGraph->AddLink(keyVal.second->mid,
                          keyVal.second->mupstream,
                          keyVal.second->mdownstream,
                          keyVal.second->mlength,
                          keyVal.second->mcosts,
                          keyVal.second->mlabel);
    }

    return newGraph;
}



OrientedGraph* mergeOrientedGraph(std::vector<const OrientedGraph*> allGraphs){
    OrientedGraph *newGraph = new OrientedGraph();

    for(auto G:allGraphs) {
        for(const auto &keyValNodes:G->mnodes) {
            mapsets excludeMovements;
            for(const auto &keyVal: keyValNodes.second->mexclude_movements) {
                setstring copy;
                for(const auto &s: keyVal.second) {
                    copy.insert(s.c_str());
                }
                excludeMovements[keyVal.first] = copy;
            }
            newGraph->AddNode(keyValNodes.first, keyValNodes.second->mposition[0], keyValNodes.second->mposition[1], keyValNodes.second->mlabel, excludeMovements);
        }

        for(const auto &keyVal: G->mlinks) {
            mapcosts costs;

            // for(const auto &keyVal: keyVal.second->mcosts) {
            //     costs[keyVal.first] = keyVal.second;
            // }
            for(const auto &keyMapCost: keyVal.second->mcosts) {
                for(const auto &keyVal: keyMapCost.second) {
                    costs[keyMapCost.first][keyVal.first] = keyVal.second;
            }
        }

            newGraph->AddLink(keyVal.first, 
                    keyVal.second->mupstream,
                    keyVal.second->mdownstream,
                    keyVal.second->mlength,
                    costs,
                    keyVal.second->mlabel);

        }
    }



    return newGraph;
}
