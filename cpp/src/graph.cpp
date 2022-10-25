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


namespace hipop
{
    /**
     * @brief Destroy the Oriented Graph:: Oriented Graph object
     * 
     */
    OrientedGraph::~OrientedGraph() {
        for (auto iter : mlinks) {
            delete iter.second;
            iter.second = NULL;
        }
        
        for (auto iter : mnodes) {
            delete iter.second;
            iter.second = NULL;
        }

        mnodes.clear();
        mlinks.clear();
    }

    /**
     * @brief Create and Add a new Node to the OrientedGraph
     * 
     * @param _id The id of the Node
     * @param x The x coordinate of the Node
     * @param y The y coordinate of the Node
     * @param label The optional label associated to the Node
     * @param excludeMovements The map of exclude movements with adjacent Nodes
     */
    void OrientedGraph::AddNode(std::string _id, double x, double y, std::string label, mapsets excludeMovements) {
        Node *new_node = new Node(_id, x, y, label, excludeMovements);
        mnodes[_id] = new_node;
    };


    /**
     * @brief Add an exisiting Node to the OrientedGraph
     * 
     * @param n 
     */
    void OrientedGraph::AddNode(Node* n) {
        mnodes[n->mid] = n;
    };


    /**
     * @brief Create and add a new Link to the OrientedGraph
     * 
     * @param _id The id of the Link
     * @param _up The id of the upstream Node of the Link
     * @param _down The id of the downstream Node of the Link
     * @param length The length of the Link
     * @param _costs The costs of the Link
     * @param label The optional label of the Link
     */
    void OrientedGraph::AddLink(std::string _id, std::string _up, std::string _down, double length, mapcosts _costs, std::string label) {
        Link *new_link = new Link(_id, _up, _down, length, _costs, label);
        mnodes[_up]->madj.emplace(_down, new_link);
        mnodes[_down]->mradj.emplace(_up, new_link);

        mlinks.emplace(_id, new_link);
    };

    /**
     * @brief Add an existing Link to the OrientedGraph
     * 
     * @param l The Link to add
     */
    void OrientedGraph::AddLink(Link *l) {
        mnodes[l->mupstream]->madj[l->mdownstream] = l;
        mnodes[l->mdownstream]->mradj[l->mupstream] = l;

        mlinks[l->mid] = l;
    };

    /**
     * @brief Update a Link costs
     * 
     * @param lid The id of the Link to update
     * @param _costs The new costs
     */
    void OrientedGraph::UpdateLinkCosts(std::string lid, mapcosts _costs) {
        mlinks[lid]->updateCosts(_costs);
    }

    /**
     * @brief Print the Nodes informations
     * 
     */
    void OrientedGraph::ShowNodes() {
        for(const auto &elem: mnodes) {
            std::cout << "Node(" << elem.first << ", [" << elem.second->mposition[0] << ",\t" << elem.second->mposition[1] << "])\n";
        }
    }

    /**
     * @brief Print the Link informations
     * 
     */
    void OrientedGraph::ShowLinks() {
        for(const auto &elem: mlinks) {
            std::cout << "Link(" << elem.first << ", " << elem.second->mupstream << ", " << elem.second->mdownstream << ")\n";
        }
    }


    /**
     * @brief Make a deep copy of an OrientedGraph
     * 
     * @param G The graph to copy
     * @return OrientedGraph* The copy
     */
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


    /**
     * @brief Merge multiple OrientedGraph together into one
     * 
     * @param allGraphs Vector of OrientedGraph to merge
     * @return OrientedGraph* The result of the merge
     */
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
        
} // namespace hipop


