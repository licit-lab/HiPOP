#include <string>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <set>
#include <vector>
#include <array>

#include <iostream>

#pragma once

typedef std::set<std::string> setstring;
typedef std::vector<std::string> vecstring;
typedef std::unordered_map<std::string, std::set<std::string> > mapsets;
typedef std::unordered_map<std::string, std::unordered_map<std::string, double> > mapcosts;

// first key : uplink, second key : downlink, third key : cost label, value : cost
typedef std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, double> > > mapnodecosts;


namespace hipop
{
    class Link {
    public:
        std::string mid;
        std::string mupstream;
        std::string mdownstream;
        mapcosts mcosts;
        std::string mlabel;
        double mlength;

        Link(std::string _id, std::string _up, std::string _down, double length, mapcosts _costs, std::string label = "") {
            mid = _id.c_str();
            mlabel = label.c_str();
            mupstream = _up.c_str();
            mdownstream = _down.c_str();
            mcosts = _costs;
            mlength = length;
        }

        Link(const Link &other) {
            mid = other.mid.c_str();
            mlabel = other.mlabel.c_str();
            mupstream = other.mupstream.c_str();
            mdownstream = other.mdownstream.c_str();
            mlength = other.mlength;

            for(const auto &keyVal: other.mcosts) {
                mcosts[keyVal.first] = keyVal.second;
            }


        }

        void updateCosts(mapcosts costs) {
            for(const auto &keyMapCost: costs) {
                for(const auto &keyVal: keyMapCost.second) {
                    mcosts[keyMapCost.first][keyVal.first] = keyVal.second;
                }
            }
        }
    };


    class Node {
    public:
        std::string mid;
        std::array<double, 2> mposition;
        std::unordered_map<std::string, Link* > madj;
        std::unordered_map<std::string, Link* > mradj;
        mapnodecosts mcosts;
        std::string mlabel;

        mapsets mexclude_movements;

        Node(std::string _id, double x, double y, std::string label = "", mapsets exclude_movements = {}, mapnodecosts costs = {}) {
            mid = _id.c_str();
            mposition[0] = x;
            mposition[1] = y;
            mexclude_movements = exclude_movements;
            mcosts = costs;
            mlabel = label.c_str();
        }

        Node(const Node &other) {
            mid = other.mid.c_str();
            mposition[0] = other.mposition[0];
            mposition[1] = other.mposition[1];

            for(const auto &keyVal: other.madj) {
                Link *l = new Link(*keyVal.second);
                madj[keyVal.first] = l;
            }

            for(const auto &keyVal: other.mradj) {
                Link *l = new Link(*keyVal.second);
                mradj[keyVal.first] = l;
            }

            for(const auto &keyVal: other.mexclude_movements) {
                setstring copy;
                for(const auto &s: keyVal.second) {
                    copy.insert(s.c_str());
                }

                mexclude_movements[keyVal.first] = copy;
            }

            mcosts = other.mcosts;
        }

        void updateCosts(mapnodecosts costs) {
            for(const auto &upLink: costs) {
                for(const auto &downLink: upLink.second) {
                    for(const auto &costLabel: downLink.second) {
                        mcosts[upLink.first][downLink.first][costLabel.first] = costLabel.second;
                    }
                }
            }
        }

        std::vector<Link*> getExits(std::string predecessor = "_default") {
            std::vector<Link*> res;
            for(const auto &l: madj) {
                std::string neighbor = l.second->mdownstream;
                if(mexclude_movements.find(predecessor) == mexclude_movements.end() || mexclude_movements[predecessor].find(neighbor) == mexclude_movements[predecessor].end()) {
                    res.push_back(l.second);
                }
            }
            return res;
        }

        std::vector<Link*> getEntrances(std::string predecessor) {
            std::vector<Link*> res;
            for(const auto &l: mradj) {
                std::string neighbor = l.second->mupstream;
                if(mexclude_movements[predecessor].find(neighbor) == mexclude_movements[predecessor].end()) {
                    res.push_back(l.second);
                }
            }
            return res;
        }

        double getCost(const std::string & upLink, const std::string & downLink, const std::string & cost) {
            auto iterUp = mcosts.find(upLink);
            if (iterUp != mcosts.end()) {
                auto iterDown = iterUp->second.find(downLink);
                if (iterDown != iterUp->second.end()) {
                    auto iterCost = iterDown->second.find(cost);
                    if (iterCost != iterDown->second.end()) {
                        return iterCost->second;
                    }
                }
            }
            return 0;
        }
    };


    class OrientedGraph {
    public:
        std::unordered_map<std::string, Node* > mnodes;
        std::unordered_map<std::string, Link* > mlinks;
        void AddNode(std::string _id, double x, double y, std::string label = "", mapsets excludeMovements = {}, mapnodecosts costs = {});
        void AddNode(Node *n);
        void AddLink(std::string _id, std::string _up, std::string _down, double length, mapcosts _costs, std::string label = "");
        void AddLink(Link* l);
        void UpdateLinkCosts(std::string lid, mapcosts _costs);
        double getLength(std::string _up, std::string _down);

        void ShowNodes();
        void ShowLinks();

        Link* getLink(std::string  _id) {
            return mlinks[_id];
        }

        OrientedGraph() {};
        OrientedGraph(const OrientedGraph &other) {
            for(const auto &keyVal: other.mnodes) {
                Node *newNode = new Node(*keyVal.second);
                mnodes[keyVal.first] = newNode;
            }

            for(const auto &keyVal: other.mlinks) {
                Link *newLink = new Link(*keyVal.second);
                mlinks[keyVal.first] = newLink;
            }
        }
        ~OrientedGraph();

    };


    OrientedGraph* copyGraph(const OrientedGraph &G);


    OrientedGraph* mergeOrientedGraph(std::vector<const OrientedGraph*> allGraphs);
    
} // namespace hipop


