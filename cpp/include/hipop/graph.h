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
        std::string mlabel;

        mapsets mexclude_movements;

        Node(std::string _id, double x, double y, std::string label = "", mapsets exclude_movements = {}) {
            mid = _id.c_str();
            mposition[0] = x;
            mposition[1] = y;
            mexclude_movements = exclude_movements;
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

    };


    class OrientedGraph {
    public:
        std::unordered_map<std::string, Node* > mnodes;
        std::unordered_map<std::string, Link* > mlinks;
        void AddNode(std::string _id, double x, double y, std::string label = "", mapsets excludeMovements = {});
        void AddNode(Node *n);
        void AddLink(std::string _id, std::string _up, std::string _down, double length, mapcosts _costs, std::string label = "");
        void AddLink(Link* l);
        void DeleteLink(std::string _id);
        void DeleteAllLinksToNode(std::string _id);
        void UpdateLinkCosts(std::string lid, mapcosts _costs);
        void UpdateCosts(std::unordered_map<std::string, mapcosts> maplinkcosts);
        double getLength(std::string _up, std::string _down);
        std::vector<std::string> GetLinksWithoutCost(std::string cost, const std::unordered_map<std::string, std::string> &mapLabelCost);

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


