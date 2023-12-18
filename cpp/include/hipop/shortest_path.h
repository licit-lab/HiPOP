#include "hipop/graph.h"

#include <vector>
#include <string>
#include <utility>
#include <functional>

#pragma once

typedef std::set<std::string> setstring;
typedef std::pair<std::vector<std::string>, double> pathCost;

namespace hipop
{
    double computePathLength(OrientedGraph &G, const std::vector<std::string> &path);

    double computePathCost(OrientedGraph &G, const std::vector<std::string> &path, std::string cost, const std::unordered_map<std::string, std::string> mapLabelCost);

    std::vector<std::vector<double>> computePathsCosts(OrientedGraph &G,
        const std::vector<std::vector<std::vector<std::string>>> &paths,
        const std::string &cost,
        const std::unordered_map<std::string, std::string> mapLabelCost,
        int threadNumber);

    pathCost dijkstra(
        const OrientedGraph &G, 
        const std::string &origin, 
        const std::string &destination, 
        const std::string &cost, 
        const std::unordered_map<std::string, std::string> &mapLabelCost, 
        setstring accessibleLabels = {});
    pathCost aStar(
        const OrientedGraph &G, 
        const std::string &origin, 
        const std::string &destination, 
        const std::string &cost,
        const std::unordered_map<std::string, std::string> &mapLabelCost, 
        const setstring &accessibleLabels, 
        std::function<double(const Node *, const Node *)> heuristic);
    pathCost aStarEuclidianDist(
        const OrientedGraph &G, 
        const std::string &origin, 
        const std::string &destination, 
        const std::string &cost,
        const std::unordered_map<std::string, std::string> &mapLabelCost, 
        const setstring &accessibleLabels);

    std::vector<pathCost> parallelDijkstra(
        const OrientedGraph &G, 
        std::vector<std::string> origins, 
        std::vector<std::string> destinations, 
        std::vector<std::unordered_map<std::string, std::string> > vecMapLabelCosts, 
        std::string cost, 
        int threadNumber, 
        std::vector<setstring> vecAvailableLabels = {});

    std::vector<pathCost> YenKShortestPath(
        OrientedGraph &G, 
        std::string origin, 
        std::string destination, 
        std::string cost, 
        setstring accessibleLabels,
        const std::unordered_map<std::string, std::string> &mapLabelCost,
        int kPath);
    std::vector<pathCost> KShortestPath(
        OrientedGraph &G, 
        const std::string &origin, 
        const std::string &destination, 
        const std::string &cost, 
        setstring accessibleLabels,
        const std::unordered_map<std::string, std::string> &mapLabelCost,
        double maxDiffCost,
        double maxDistInCommon,
        int costMultiplier,
        int maxRetry,
        int kPath,
        bool intermodal);

    std::vector<std::vector<pathCost>> parallelKShortestPath(
        OrientedGraph &G, 
        const std::vector<std::string> &origins, 
        const std::vector<std::string> &destinations, 
        const std::string &cost,
        const std::vector<std::unordered_map<std::string, std::string> > vecMapLabelCosts,
        const std::vector<setstring> accessibleLabels, 
        double maxDiffCost,
        double maxDistInCommon,
        int costMultiplier,
        int maxRetry,
        const std::vector<int> &kPaths, 
        int threadNumber);

    std::vector<std::vector<pathCost>> parallelKIntermodalShortestPath(
        const OrientedGraph &G,
        std::vector<std::string> origins,
        std::vector<std::string> destinations,
        std::vector<std::unordered_map<std::string, std::string> > vecMapLabelCosts,
        std::string cost,
        int threadNumber,
        std::pair<std::unordered_set<std::string>, std::unordered_set<std::string>> pairMandatoryLabels,
        double maxDiffCost,
        double maxDistInCommon,
        int costMultiplier,
        int maxRetry,
        std::vector<int> kPaths,
        std::vector<setstring> vecAvailableLabels = {});

} // namespace hipop
