#include "hipop/graph.h"

#include <vector>
#include <string>
#include <utility>
#include <functional>

#pragma once

typedef std::set<std::string> setstring;
typedef std::pair<std::vector<std::string>, double> pathCost;

double computePathLength(OrientedGraph &G, const std::vector<std::string> &path);

pathCost dijkstra(const OrientedGraph &G, const std::string &origin, const std::string &destination, const std::string &cost, setstring accessibleLabels = {});
pathCost aStar(const OrientedGraph &G, const std::string &origin, const std::string &destination, const std::string &cost, const setstring &accessibleLabels, std::function<double(const Node *, const Node *)> heuristic);
pathCost aStarEuclidianDist(const OrientedGraph &G, const std::string &origin, const std::string &destination, const std::string &cost, const setstring &accessibleLabels);

std::vector<pathCost> parallelDijkstra(const OrientedGraph &G, std::vector<std::string> origins, std::vector<std::string> destinations, std::string cost, int threadNumber, std::vector<setstring> vecAvailableLabels = {});

std::vector<pathCost> YenKShortestPath(OrientedGraph &G, std::string origin, std::string destination, std::string cost, setstring accessibleLabels, int kPath);
std::vector<pathCost> KShortestPath(OrientedGraph &G, const std::string &origin, const std::string &destination, const std::string &cost, setstring accessibleLabels, double minDist, double maxDist, int kPath);

std::vector<std::vector<pathCost>> parallelKShortestPath(OrientedGraph &G, const std::vector<std::string> &origins, const std::vector<std::string> &destinations, const std::string &cost, const std::vector<setstring> accessibleLabels, double minDist, double maxDist, int kPath, int threadNumber);
