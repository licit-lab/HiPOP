#include "hipop/graph.h"
#include "hipop/shortest_path.h"

#include <omp.h>

#include <vector>
#include <queue>
#include <unordered_map>
#include <string>
#include <numeric>
#include <limits>
#include <algorithm>
#include <functional>
#include <iostream>
#include <cmath>

typedef std::pair<double, std::string> QueueItem;
typedef std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem>> PriorityQueue;


namespace hipop
{
    /**
     * @brief Compute the shortest path between origin and destination using the Dijkstra algorithm
     * 
     * @param G The OrientedGrah used for the shortest path
     * @param origin The origin 
     * @param destination The destination
     * @param cost The costs to consider in the shortest path algorithm
     * @param mapLabelCost The type of cost map to choose on each label (mulitple set of costs can be defined on a Link)
     * @param accessibleLabels The set of accessible label
     * @return pathCost The list of Nodes defining the shortest path and the associated cost
     */
    pathCost dijkstra(
        const OrientedGraph &G, 
        const std::string &origin, 
        const std::string &destination, 
        const std::string &cost, 
        const std::unordered_map<std::string, std::string> &mapLabelCost, 
        setstring accessibleLabels)
    {
        pathCost path;

        PriorityQueue pq;

        std::unordered_map<std::string, double> dist;
        std::unordered_map<std::string, std::string> prev;
        prev.reserve(G.mnodes.size());
        dist.reserve(G.mnodes.size());
        double inf = std::numeric_limits<double>::infinity();
        for (const auto keyVal : G.mnodes)
        {
            dist[keyVal.first] = inf;
        }
        pq.push(make_pair(0, origin));
        dist[origin] = 0;

        path.second = inf;
        prev[origin] = "";

        if (origin==destination) {
        path.second = 0;
        return path;
        }

        while (!pq.empty())
        {
            QueueItem current = pq.top();
            pq.pop();
            std::string u = current.second;

            if (u == destination)
            {
                std::string v = prev[u];
                path.first.push_back(u);

                while (v != origin)
                {
                    path.first.push_back(v);
                    v = prev[v];
                }

                path.first.push_back(v);
                std::reverse(path.first.begin(), path.first.end());
                path.second = dist[destination];
                return path;
            }

            try
            {
                for (const auto link : G.mnodes.at(u)->getExits(prev[u]))
                {
                    if (accessibleLabels.empty() || accessibleLabels.find(link->mlabel) != accessibleLabels.end())
                    {
                        std::string neighbor = link->mdownstream;
                        double new_dist = dist[u] + link->mcosts[mapLabelCost.at(link->mlabel)][cost];

                        if (dist[neighbor] > new_dist)
                        {
                            dist[neighbor] = new_dist;
                            pq.push(QueueItem(new_dist, neighbor));
                            prev[neighbor] = u;
                        }
                    }
                }
            }
            catch(const std::out_of_range& e)
            {
                std::cerr <<  "The node " << u << " does not belong to the graph \n";
            }
            
            
        }
        return path;
    }

    /**
     * @brief Batch computation of shortest with openmp using the Dijkstra algorithm
     * 
     * @param G The OrientedGrah used for the shortest paths
     * @param origins The vector of origins
     * @param destinations The vector of destinations
     * @param vecMapLabelCosts The vector of type of cost map to choose on each label
     * @param cost The cost to consider in the shortest path algorithm
     * @param threadNumber The number of thread for openmp
     * @param vecAvailableLabels The vector of available labels
     * @return std::vector<pathCost> The vector of computed shortest path
     */
    std::vector<pathCost> parallelDijkstra(
        const OrientedGraph &G, 
        std::vector<std::string> origins, 
        std::vector<std::string> destinations,
        std::vector<std::unordered_map<std::string, std::string> > vecMapLabelCosts,
        std::string cost, 
        int threadNumber, 
        std::vector<setstring> vecAvailableLabels)
    {
        omp_set_num_threads(threadNumber);

        int nbPath = origins.size();
        std::vector<pathCost> res(nbPath);

        #pragma omp parallel for shared(res, vecAvailableLabels, vecMapLabelCosts) schedule(dynamic)
        for (int i = 0; i < nbPath; i++)
        {
            if (vecAvailableLabels.empty())
            {
                res[i] = dijkstra(G, origins[i], destinations[i], cost, vecMapLabelCosts[i], {});
            }
            else
            {
                res[i] = dijkstra(G, origins[i], destinations[i], cost, vecMapLabelCosts[i], vecAvailableLabels[i]);
            }
        }

        return res;
    }

    typedef std::unordered_map<std::string, mapcosts> linkMapCosts;

    /**
     * @brief Increase the cost in a OrientedGraph for a path
     * 
     * @param G The OrientedGraph on which the increase occurs
     * @param path The path where the costs should be increased
     * @param initial_costs The intial cost of the links to save
     * @param costMultiplier The multiplier to apply to the links costs of the path
     */
    void increaseCostsFromPath(OrientedGraph &G, const std::vector<std::string> &path, linkMapCosts &initial_costs, double costMultiplier)
    {

        for (size_t i = 0; i < path.size() - 1; i++)
        {
            Link *link = G.mnodes[path[i]]->madj[path[i + 1]];
            if (initial_costs.find(link->mid) == initial_costs.end())
            {
                for (auto &keyMapCost : link->mcosts)
                {
                    for(auto &keyVal: keyMapCost.second) {
                        initial_costs[link->mid][keyMapCost.first][keyVal.first] = keyVal.second; // save initial cost
                        keyVal.second *= costMultiplier;
                    }
                }
            }
            else
            {
                for(auto &keyMapCost: link->mcosts) {
                    for(auto &keyVal: keyMapCost.second) {
                        keyVal.second *= costMultiplier;
                    }
                }
            }
        }
    }

    /**
     * @brief Increase the cost in a OrientedGraph for an intermodal path
     *
     * @param G The OrientedGraph on which the increase occurs
     * @param path The path where the costs should be increased
     * @param initial_costs The intial cost of the links to save
     * @param costMultiplier The multiplier to apply to the links costs of the path
     */
    void increaseCostsFromIntermodalPath(OrientedGraph &G, const std::vector<std::string> &path, linkMapCosts &initial_costs, double costMultiplier)
    {

        for (size_t i = 0; i < path.size() - 1; i++)
        {
            Link *link = G.mnodes[path[i]]->madj[path[i + 1]];
            std::string decoded_link_id = "";

            if (link->mid.size() > 5 && (link->mid.compare(link->mid.size() - 5, 5, "_TWIN") == 0 || link->mid.compare(link->mid.size() - 5, 5, "_TRPL") == 0))
            {
              decoded_link_id = link->mid.substr(0, link->mid.size() - 5);
            }
            else if (link->mid.size() > 17 && (link->mid.compare(link->mid.size() - 17, 17, "_ORIGINAL_TO_TWIN") == 0))
            {
                decoded_link_id = link->mid.substr(0, link->mid.size() - 17);
            }
            else if (link->mid.size() > 13 && (link->mid.compare(link->mid.size() - 13, 13, "_ORIGINAL_TO_TWIN") == 0))
            {
                decoded_link_id = link->mid.substr(0, link->mid.size() - 13);
            }
            else
            {
                decoded_link_id = link->mid;
            }

            std::vector<std::string> corresponding_links_ids = {decoded_link_id, decoded_link_id + "_TWIN", decoded_link_id + "_TRPL", decoded_link_id + "_ORIGINAL_TO_TWIN", decoded_link_id + "_TWIN_TO_TRPL"};

            if (initial_costs.find(decoded_link_id) == initial_costs.end())
            {
                // Increase cost of all corresponding links and save initial costs
                for (int j = 0; j < 5; ++j)
                {
                    std::string corresponding_link_id = corresponding_links_ids[j];
                    if (G.mlinks.find(corresponding_link_id) != G.mlinks.end())
                    {
                        Link *corresponding_link = G.mlinks[corresponding_link_id];
                        for (auto &keyMapCost : corresponding_link->mcosts)
                        {
                            for(auto &keyVal: keyMapCost.second) {
                                initial_costs[corresponding_link_id][keyMapCost.first][keyVal.first] = keyVal.second;
                                keyVal.second *= costMultiplier;
                            }
                        }
                    }
                }
            }
            else
            {
                // Only increase cost of all corresponding links
                for (int j = 0; j < 5; ++j)
                {
                    std::string corresponding_link_id = corresponding_links_ids[j];
                    if (G.mlinks.find(corresponding_link_id) != G.mlinks.end())
                    {
                        Link *corresponding_link = G.mlinks[corresponding_link_id];
                        for (auto &keyMapCost : corresponding_link->mcosts)
                        {
                            for(auto &keyVal: keyMapCost.second) {
                                keyVal.second *= costMultiplier;
                            }
                        }
                    }
                }
            }
        }
    }


    /**
     * @brief Compute the length of a path
     * 
     * @param G The OrientedGraph on which the path is computed
     * @param path The path
     * @return double The length of the path
     */
    double computePathLength(OrientedGraph &G, const std::vector<std::string> &path)
    {
        double length = 0;

        for (size_t i = 0; i < path.size() - 1; i++)
        {
            Link *link = G.mnodes[path[i]]->madj[path[i + 1]];
            length += link->mlength;
        }

        return length;
    }

    /**
     * @brief Compute the relative distances in common between path and paths
     *
     * @param G The OrientedGraph
     * @param path The path to compare to paths
     * @param paths The paths
     * @return std::vector<double> The relative distances in common
     */
    std::vector<double> computeRelativeDistancesInCommon(OrientedGraph &G, const std::vector<std::string> &path, std::vector<pathCost> &paths)
    {
        //assert path.size() > 0;
        int nbPaths = paths.size();
        std::vector<double> relDists(nbPaths);

        for (int i = 0; i < nbPaths; i++)
        {
            std::vector<std::string> compared_p = paths[i].first;
            //assert compared_p.size() > 0;
            double path_length = computePathLength(G, path);
            double compared_p_length = computePathLength(G, compared_p);
            int compared_p_size = compared_p.size();
            std::vector<std::string> compared_p_links(compared_p_size);
            for (int j = 0; j < compared_p_size - 1; j++)
            {
                Link *link = G.mnodes[compared_p[j]]->madj[compared_p[j + 1]];
                compared_p_links[j] = link->mid;
            }
            double commonDist = 0;
            for (int j = 0; j < path.size() - 1; j++)
            {
                Link *link = G.mnodes[path[j]]->madj[path[j + 1]];
                if (std::find(compared_p_links.begin(), compared_p_links.end(), link->mid) != compared_p_links.end())
                {
                  commonDist += link->mlength;
                }
            }
            relDists[i] = commonDist / fmax(path_length, compared_p_length);
        }

        return relDists;
    }

    /**
     * @brief Batch computation of paths costs
     *
     * @param G The OrientedGraph on which the path is computed
     * @param paths The paths grouped in batches
     * @param cost The cost to consider
     * @param mapLabelCost The type of cost map to choose on each label
     * @param threadNumber Number of threads to use
     * @return std::vector<double> The total costs of the paths
     */
    std::vector<std::vector<double>> computePathsCosts(
        OrientedGraph &G,
        const std::vector<std::vector<std::vector<std::string>>> &paths,
        const std::string &cost,
        const std::unordered_map<std::string, std::string> mapLabelCost,
        int threadNumber)
    {
        omp_set_num_threads(threadNumber);
        int nbBatches = paths.size();

        std::vector<std::vector<double>> res(nbBatches);
        OrientedGraph *privateG;

        #pragma omp parallel shared(res, G, paths, cost, mapLabelCost) private(privateG)
        {
            privateG = copyGraph(G);

            #pragma omp for
            for (int i = 0; i < nbBatches; i++)
            {
                int nbPaths = paths[i].size();
                std::vector<double> res_(nbPaths);
                for (int j = 0; j < nbPaths; j++)
                {
                  res_[j] = computePathCost(*privateG, paths[i][j], cost, mapLabelCost);
                }
                res[i] = res_;
            }

            // Not sure if the omp critical is necessary
            #pragma omp critical
            {
                delete privateG;
            }

        }

        return res;
    }

    /**
     * @brief Compute the total cost of a path
     * 
     * @param G The OrientedGraph on which the path is computed
     * @param path The path
     * @param cost The cost to consider 
     * @param mapLabelCost The type of cost map to choose on each label
     * @return double The total cost of the path
     */
    double computePathCost(OrientedGraph &G, const std::vector<std::string> &path, std::string cost, const std::unordered_map<std::string, std::string> mapLabelCost)
    {
        double c = 0;

        if (path.size() == 0)
        {
          return c;
        }
        else
        {
          for (size_t i = 0; i < path.size() - 1; i++)
          {
              Link *link = G.mnodes[path[i]]->madj[path[i + 1]];
              c += link->mcosts[mapLabelCost.at(link->mlabel)][cost];
          }
          return c;
        }
    }

    /**
     * @brief Compute the total cost of a path using map for effective cost on some links of the graph.
     *
     * @param G The OrientedGraph on which the path is computed
     * @param path The path
     * @param cost The cost to consider
     * @param mapLabelCost The type of cost map to choose on each label
     * @param initialCosts The effective costs values to use for some links
     * @return double The total cost of the path
     */
    double computePathCostWithInitialCostsDict(OrientedGraph &G, const std::vector<std::string> &path, std::string cost, const std::unordered_map<std::string, std::string> mapLabelCost, linkMapCosts initialCosts)
    {
        double c = 0;

        if (path.size() == 0)
        {
          return c;
        }
        else
        {
          for (size_t i = 0; i < path.size() - 1; i++)
          {
              Link *link = G.mnodes[path[i]]->madj[path[i + 1]];
              if (initialCosts.find(link->mid) != initialCosts.end())
              {
                  c += initialCosts[link->mid][mapLabelCost.at(link->mlabel)][cost];
              }
              else
              {
                  c += link->mcosts[mapLabelCost.at(link->mlabel)][cost];
              }
          }
          return c;
        }
    }

    /**
     * @brief Print a path
     * 
     * @param path The path to print
     */
    void showPath(pathCost path)
    {
        std::cout << path.second << " [";
        for (const auto &p : path.first)
        {
            std::cout << p << ", ";
        }
        std::cout << "]" << std::endl;
    }

    /**
     * @brief Print a path
     *
     * @param path The path to print
     */
    void showPathNodes(std::vector<std::string> path)
    {
        std::cout << " [";
        for (const auto &p : path)
        {
            std::cout << p << ", ";
        }
        std::cout << "]" << std::endl;
    }


    /**
     * @brief Function that decodes an intermodal path computed on a tripled graph
     *
     * @param path The path to decode
     * @return decodedPath The decoded path
     */
    std::vector<std::string> decodeIntermodalPath(std::vector<std::string> path)
    {
        int pathSize = path.size();
        std::vector<std::string> decodedPath(pathSize);
        for (int k = 0; k < pathSize; k++) {
            if (path[k].size() > 5 && (path[k].compare(path[k].size() - 5, 5, "_TWIN") == 0 || path[k].compare(path[k].size() - 5, 5, "_TRPL") == 0)) {
                decodedPath[k] = path[k].substr(0, path[k].size() - 5);
            }
            else
            {
                decodedPath[k] = path[k];
            }
        }
        return decodedPath;
    }

    /**
     * @brief Compute K shortest paths for an origin/destination. The first path is normally computed, then we increase the costs of the first path and perform a Dijkstra again.
     * The last step is performed until either we reach K accepted shortest paths or the last computed shortest path do not meet the requirements in terms of maxDiffCost and maxDistInCommon.
     *
     * @param G The OrientedGraph on which we compute the paths
     * @param origin The origin
     * @param destination The destination
     * @param cost The cost to consider
     * @param accessibleLabels The set of accessible label
     * @param mapLabelCost The type of cost map to choose on each label
     * @param maxDiffCost The maximal difference between the cost of the first computed
     *        shortest path and the cost of the next ones, expressed as a percentage
     *        (e.g. 0.1 means that the cost of the next path should be less than 101%
     *        of the cost of the first computed shortest path)
     * @param maxDistInCommon The maximal distance in common between the first shortest
     *        path found and the next ones, expressed as a percentage (e.g. 0.6 means that
     *        the next path should have less than 60% common distance with the first computed
     *        shortest path to be accepted)
     * @param costMultiplier The multiplier applied to the links costs of an
     *        accepted shortest path
     * @param maxRetry Maximum number of times we retry to find an acceptable shorest path
     * @param kPath The number of paths to compute
     * @param intermodal Specifies we search k intermodal shortest paths on a tripled graph
     * @return std::vector<pathCost> The vector of k computed paths
     */
     std::vector<pathCost> KShortestPath(
        OrientedGraph &G,
        const std::string &origin,
        const std::string &destination,
        const std::string &cost,
        setstring accessibleLabels,
        const std::unordered_map<std::string, std::string> &mapLabelCost,
        double maxDiffCost,
        double maxDistInCommon,
        double costMultiplier,
        int maxRetry,
        int kPath,
        bool intermodal)
    {
        //assert (maxDiffCost >= 0);
        //assert (maxDistInCommon >= 0 and maxDistInCommon <= 1);

        std::vector<pathCost> paths;
        linkMapCosts initial_costs;

        pathCost firstPath = dijkstra(G, origin, destination, cost, mapLabelCost, accessibleLabels);
        paths.push_back(firstPath);

        if (firstPath.first.empty()) // no path found
        {
            return paths;
        }

        if (intermodal)
        {
            increaseCostsFromIntermodalPath(G, firstPath.first, initial_costs, costMultiplier);
        }
        else
        {
            increaseCostsFromPath(G, firstPath.first, initial_costs, costMultiplier);
        }

        int pathCounter = 1, retry = 0;

        while (pathCounter < kPath && retry < maxRetry )
        {
            pathCost newPath = dijkstra(G, origin, destination, cost, mapLabelCost, accessibleLabels);
            newPath.second = computePathCostWithInitialCostsDict(G, newPath.first, cost, mapLabelCost, initial_costs);

            if (newPath.first.empty())
            {
                break; // no path found
            }

            // Check conditions to accept this new path
            std::vector<double> relDistancesInCommon = computeRelativeDistancesInCommon(G, newPath.first, paths); // relative distances in common between newPath and the already found ones
            bool maxDistInCommonChecked = (std::all_of(relDistancesInCommon.cbegin(), relDistancesInCommon.cend(), [maxDistInCommon](double rd){ return rd  <= maxDistInCommon; }));
            bool maxDiffCostChecked = ( (newPath.second - firstPath.second) / firstPath.second <= maxDiffCost);
            bool isNew = true;
            if (intermodal)
            {
                isNew = (std::all_of(paths.cbegin(), paths.cend(), [newPath](pathCost p){ return decodeIntermodalPath(p.first) != decodeIntermodalPath(newPath.first); }));
            }
            else
            {
                isNew = (std::all_of(paths.cbegin(), paths.cend(), [newPath](pathCost p){ return p.first != newPath.first; }));
            }

            if (maxDistInCommonChecked && maxDiffCostChecked && isNew)
            {
                paths.push_back(newPath);
                retry = 0;
                pathCounter += 1;
            }
            else
            {
                retry += 1;
            }

            if (intermodal)
            {
                increaseCostsFromIntermodalPath(G, newPath.first, initial_costs, costMultiplier);
            }
            else
            {
                increaseCostsFromPath(G, newPath.first, initial_costs, costMultiplier);
            }
        }

        // Reset initial costs /!\ Take care if this function is called in parallel !!
        for (const auto &keyVal : initial_costs)
        {
            G.mlinks[keyVal.first]->mcosts = keyVal.second;
        }

        return paths;
    }

    /**
     * @brief Compute K shortest path using the Yen algorithm
     * 
     * @param G The OrientedGraph on which we compute the paths
     * @param origin The origin
     * @param destination The destination
     * @param cost The cost to consider
     * @param accessibleLabels The set of accessible label
     * @param mapLabelCost The type of cost map to choose on each label
     * @param kPath The number of path to compute
     * @return std::vector<pathCost> The vector of computed paths
     */
    std::vector<pathCost> YenKShortestPath(
        OrientedGraph &G, 
        std::string origin, 
        std::string destination, 
        std::string cost, 
        setstring accessibleLabels,
        const std::unordered_map<std::string, std::string> &mapLabelCost, 
        int kPath)
    {
        std::vector<pathCost> A;
        std::vector<pathCost> B;
        A.push_back(dijkstra(G, origin, destination, cost, mapLabelCost, accessibleLabels));

        double inf = std::numeric_limits<double>::infinity();

        for (size_t k = 1; k < kPath; k++)
        {
            for (size_t i = 0; i < A[k - 1].first.size() - 2; i++)
            {
                std::unordered_map<std::string, std::pair<std::string, double> > initial_costs;
                std::string spurNode = A[k - 1].first[i];
                pathCost rootPath;
                rootPath.second = 0;
                rootPath.first.insert(rootPath.first.begin(), A[k - 1].first.begin(), A[k - 1].first.begin() + i + 1);

                for (int j = 0; j < rootPath.first.size() - 1; j++)
                {
                    Link *l = G.mnodes[rootPath.first[j]]->madj[rootPath.first[j + 1]];
                    rootPath.second += l->mcosts[mapLabelCost.at(l->mlabel)][cost];
                }

                for (const pathCost &pc : A)
                {
                    if (std::equal(pc.first.begin(), pc.first.begin() + i, rootPath.first.begin()))
                    {
                        Link *l = G.mnodes[pc.first[i]]->madj[pc.first[i + 1]];

                        if (initial_costs.find(l->mid) == initial_costs.end())
                        {
                            initial_costs[l->mid] = {mapLabelCost.at(l->mlabel), l->mcosts[mapLabelCost.at(l->mlabel)][cost]};
                        }
                        l->mcosts[mapLabelCost.at(l->mlabel)][cost] = inf;
                    }
                }

                pathCost spurPath = dijkstra(G, spurNode, destination, cost, mapLabelCost, accessibleLabels);
                pathCost totalPath;
                totalPath.first = rootPath.first;

                totalPath.first.insert(totalPath.first.end(), spurPath.first.begin() + 1, spurPath.first.end());
                totalPath.second = rootPath.second + spurPath.second;

                for (const auto &keyVal : initial_costs)
                {
                    G.mlinks[keyVal.first]->mcosts[keyVal.second.first][cost] = keyVal.second.second;
                }

                bool toAdd = true;
                for (const auto &prevPath : B)
                {
                    if (totalPath.first == prevPath.first)
                    {
                        toAdd = false;
                        break;
                    }
                }

                if (toAdd)
                {
                    B.push_back(totalPath);
                }
            }

            if (B.empty())
            {
                break;
            }

            std::sort(B.begin(), B.end(), [](pathCost a, pathCost b)
                    { return a.second < b.second; });
            A.push_back(B[0]);
            B.erase(B.begin());
        }

        return A;
    }

    /**
     * @brief Batch computation of K shortest paths using openmp, each thread has its own deep copy of the OrientedGraph to ensure that the increase of the cost do not collapse with the other threads
     * 
     * @param G The OrientedGraph on which we compute the paths
     * @param origins The origins
     * @param destinations The destinations
     * @param cost The cost to consider
     * @param vecMapLabelCosts The vector of type of cost map to choose on each label
     * @param accessibleLabels The vector set of accessible label
     * @param maxDiffCost The maximal difference between the cost of the first computed
     *                    shortest path and the cost of the next ones
     * @param maxDistInCommon The maximal distance in common between the first shortest
     *                        path found and the next ones
     * @param costMultiplier The multiplier applied to the links costs of an
     *                       accepted shortest path
     * @param maxRetry Maximum number of times we retry to find an acceptable shorest path
     * @param kPaths The number of paths to compute
     * @param threadNumber Number of threads to use
     * @return std::vector<std::vector<pathCost>>
     */
    std::vector<std::vector<pathCost>> parallelKShortestPath(
        OrientedGraph &G, 
        const std::vector<std::string> &origins, 
        const std::vector<std::string> &destinations, 
        const std::string &cost,
        const std::vector<std::unordered_map<std::string, std::string> > vecMapLabelCosts,
        const std::vector<setstring> accessibleLabels,
        double maxDiffCost,
        double maxDistInCommon,
        double costMultiplier,
        int maxRetry,
        const std::vector<int> &kPaths,
        int threadNumber)
    {
        omp_set_num_threads(threadNumber);
        int nbOD = origins.size();

        std::vector<std::vector<pathCost>> res(nbOD);
        OrientedGraph *privateG;

        #pragma omp parallel shared(res, accessibleLabels, G, vecMapLabelCosts, origins, destinations, kPaths) private(privateG)
        {
            privateG = copyGraph(G);

            #pragma omp for
            for (int i = 0; i < nbOD; i++)
            {
                if (accessibleLabels.empty())
                {
                    res[i] = KShortestPath(*privateG, origins[i], destinations[i], cost, {}, vecMapLabelCosts[i], maxDiffCost, maxDistInCommon, costMultiplier, maxRetry, kPaths[i], false);
                }
                else
                {
                    res[i] = KShortestPath(*privateG, origins[i], destinations[i], cost, accessibleLabels[i], vecMapLabelCosts[i], maxDiffCost, maxDistInCommon, costMultiplier, maxRetry, kPaths[i], false);
                }
            }
            
            // Not sure if the omp critical is necessary
            #pragma omp critical
            {
                delete privateG;
            }

        }
        
        return res;
    }

    /**
     * @brief Compute the shortest between origin and destination using the A* algorithm
     * 
     * @param G The OrientedGraph on which we compute the path
     * @param origin The origin
     * @param destination The destination
     * @param cost The cost to consider in the shortest path algoritm
     * @param mapLabelCost The type of cost map to choose on each label
     * @param accessibleLabels The set of accessible label
     * @param heuristic An heuristic to speed up the shortest path computation
     * @return pathCost The computed path
     */
    pathCost aStar(
        const OrientedGraph &G, 
        const std::string &origin, 
        const std::string &destination, 
        const std::string &cost, 
        const std::unordered_map<std::string, std::string> &mapLabelCost,
        const setstring &accessibleLabels,
        std::function<double(const Node *, const Node *)> heuristic)
    {
        pathCost path;

        PriorityQueue pq;

        std::unordered_map<std::string, double> dist;
        std::unordered_map<std::string, std::string> prev;
        prev.reserve(G.mnodes.size());
        dist.reserve(G.mnodes.size());
        double inf = std::numeric_limits<double>::infinity();
        for (const auto keyVal : G.mnodes)
        {
            dist[keyVal.first] = inf;
        }
        pq.push(make_pair(0, origin));
        dist[origin] = 0;

        path.second = inf;
        prev[origin] = "";

        if (origin == destination){
        path.second = 0;
        return path;
        }

        while (!pq.empty())
        {
            std::string u = pq.top().second;
            pq.pop();

            if (u == destination)
            {
                std::string v = prev[u];
                path.first.push_back(u);

                while (v != origin)
                {
                    path.first.push_back(v);
                    v = prev[v];
                }

                path.first.push_back(v);
                std::reverse(path.first.begin(), path.first.end());
                path.second = dist[destination];
                return path;
            }

            for (const auto link : G.mnodes.at(u)->getExits(prev[u]))
            {
                if (accessibleLabels.empty() || accessibleLabels.find(link->mlabel) != accessibleLabels.end())
                {
                    std::string neighbor = link->mdownstream;
                    double tentative_score = dist[u] + link->mcosts[mapLabelCost.at(link->mlabel)][cost];

                    if (tentative_score < dist[neighbor])
                    {
                        dist[neighbor] = tentative_score;
                        pq.push(QueueItem(tentative_score + heuristic(G.mnodes.at(u), G.mnodes.at(destination)), neighbor));
                        prev[neighbor] = u;
                    }
                }
            }
        }
        return path;
    }

    /**
     * @brief A simple heuristic for the A* based on the euclidian distance between two nodes
     * 
     * @param current The current Node
     * @param dest The destination Node
     * @return double The distance between current and dest
     */
    double euclidianDist(const Node *current, const Node *dest)
    {
        return std::sqrt(std::pow(dest->mposition[0] - current->mposition[0], 2) + std::pow(dest->mposition[1] - current->mposition[1], 2));
    }

    /**
     * @brief A* algorithm with an euclidian distance as heuristic
     * 
     * @param G The OrientedGraph on which we compute the path
     * @param origin The origin
     * @param destination The destination
     * @param cost The cost to consider in the shortest path algoritm 
     * @param mapLabelCost The type of cost map to choose on each label
     * @param accessibleLabels The set of accessible label
     * @return pathCost The computed path
     */
    pathCost aStarEuclidianDist(
        const OrientedGraph &G, 
        const std::string &origin, 
        const std::string &destination, 
        const std::string &cost,
        const std::unordered_map<std::string, std::string> &mapLabelCost, 
        const setstring &accessibleLabels)
    {
        return aStar(G, origin, destination, cost, mapLabelCost, accessibleLabels, euclidianDist);
    }

    /**
     * @brief Batch computation of intermodal shortest paths with openmp using the
     *        a doubled graph to be sure that the path found is intermodal and Dijkstra algorithm
     *
     * @param G The OrientedGrah used for the shortest paths
     * @param origins The vector of origins
     * @param destinations The vector of destinations
     * @param vecMapLabelCosts The vector of type of cost map to choose on each label
     * @param cost The cost to consider in the shortest path algorithm
     * @param threadNumber The number of thread for openmp
     * @param pairMandatoryLabels The pair of labels groups the shortest paths must contain
     * @param kPaths The number of paths to compute
     * @param maxDiffCost The maximal difference between the cost of the first computed
     *                    shortest path and the cost of the next ones
     * @param maxDistInCommon The maximal distance in common between the first shortest
     *                        path found and the next ones
     * @param costMultiplier The multiplier applied to the links costs of an
     *                       accepted shortest path
     * @param maxRetry Maximum number of times we retry to find an acceptable shorest path
     * @param vecAvailableLabels The vector of available labels
     * @return std::vector<std::vector<pathCost>> The vector of computed shortest path
     */
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
        double costMultiplier,
        int maxRetry,
        std::vector<int> kPaths,
        std::vector<setstring> vecAvailableLabels)
    {
        // Create doubled graph two ways
        OrientedGraph *doubledG1 = new OrientedGraph(); // pass first on first elem of pairMandatoryLabels
        OrientedGraph *doubledG2 = new OrientedGraph(); // pass first on second elem of pairMandatoryLabels
        for(const auto &keyVal: G.mnodes) {
            doubledG1->AddNode(keyVal.second->mid,
                            keyVal.second->mposition[0],
                            keyVal.second->mposition[1],
                            keyVal.second->mlabel,
                            keyVal.second->mexclude_movements);
            doubledG1->AddNode(keyVal.second->mid + "_TWIN",
                            keyVal.second->mposition[0],
                            keyVal.second->mposition[1],
                            keyVal.second->mlabel,
                            keyVal.second->mexclude_movements);
            doubledG1->AddNode(keyVal.second->mid + "_TRPL",
                            keyVal.second->mposition[0],
                            keyVal.second->mposition[1],
                            keyVal.second->mlabel,
                            keyVal.second->mexclude_movements);
            doubledG2->AddNode(keyVal.second->mid,
                            keyVal.second->mposition[0],
                            keyVal.second->mposition[1],
                            keyVal.second->mlabel,
                            keyVal.second->mexclude_movements);
            doubledG2->AddNode(keyVal.second->mid + "_TWIN",
                            keyVal.second->mposition[0],
                            keyVal.second->mposition[1],
                            keyVal.second->mlabel,
                            keyVal.second->mexclude_movements);
            doubledG2->AddNode(keyVal.second->mid + "_TRPL",
                            keyVal.second->mposition[0],
                            keyVal.second->mposition[1],
                            keyVal.second->mlabel,
                            keyVal.second->mexclude_movements);
        }
        for(const auto &keyVal: G.mlinks) {
          doubledG1->AddLink(keyVal.second->mid,
                          keyVal.second->mupstream,
                          keyVal.second->mdownstream,
                          keyVal.second->mlength,
                          keyVal.second->mcosts,
                          keyVal.second->mlabel);
          doubledG1->AddLink(keyVal.second->mid + "_TWIN",
                          keyVal.second->mupstream + "_TWIN",
                          keyVal.second->mdownstream + "_TWIN",
                          keyVal.second->mlength,
                          keyVal.second->mcosts,
                          keyVal.second->mlabel);
          doubledG1->AddLink(keyVal.second->mid + "_TRPL",
                          keyVal.second->mupstream + "_TRPL",
                          keyVal.second->mdownstream + "_TRPL",
                          keyVal.second->mlength,
                          keyVal.second->mcosts,
                          keyVal.second->mlabel);
          doubledG2->AddLink(keyVal.second->mid,
                          keyVal.second->mupstream,
                          keyVal.second->mdownstream,
                          keyVal.second->mlength,
                          keyVal.second->mcosts,
                          keyVal.second->mlabel);
          doubledG2->AddLink(keyVal.second->mid + "_TWIN",
                          keyVal.second->mupstream + "_TWIN",
                          keyVal.second->mdownstream + "_TWIN",
                          keyVal.second->mlength,
                          keyVal.second->mcosts,
                          keyVal.second->mlabel);
          doubledG2->AddLink(keyVal.second->mid + "_TRPL",
                          keyVal.second->mupstream + "_TRPL",
                          keyVal.second->mdownstream + "_TRPL",
                          keyVal.second->mlength,
                          keyVal.second->mcosts,
                          keyVal.second->mlabel);
          // Add link on G1 between original and twin graph only if it corresponds to a mandatory label of the first group
          bool original_to_twin_G1 = pairMandatoryLabels.first.count(keyVal.second->mlabel);
          if (original_to_twin_G1) {
            doubledG1->AddLink(keyVal.second->mid + "_ORIGINAL_TO_TWIN",
                            keyVal.second->mupstream,
                            keyVal.second->mdownstream + "_TWIN",
                            keyVal.second->mlength,
                            keyVal.second->mcosts,
                            keyVal.second->mlabel);
          }
          // Add link on G2 between original and twin graph only if it corresponds to a mandatory label of the second group
          bool original_to_twin_G2 = pairMandatoryLabels.second.count(keyVal.second->mlabel);
          if (original_to_twin_G2) {
            doubledG2->AddLink(keyVal.second->mid + "_ORIGINAL_TO_TWIN",
                            keyVal.second->mupstream,
                            keyVal.second->mdownstream + "_TWIN",
                            keyVal.second->mlength,
                            keyVal.second->mcosts,
                            keyVal.second->mlabel);
          }
          // Add link on G1 between twin and triple graph only if it corresponds to a mandatory label of the second group
          bool twin_to_trpl_G1 = pairMandatoryLabels.second.count(keyVal.second->mlabel);
          if (twin_to_trpl_G1) {
            doubledG1->AddLink(keyVal.second->mid + "_TWIN_TO_TRPL",
                            keyVal.second->mupstream + "_TWIN",
                            keyVal.second->mdownstream + "_TRPL",
                            keyVal.second->mlength,
                            keyVal.second->mcosts,
                            keyVal.second->mlabel);
          }
          // Add link on G2 between twin and triple graph only if it corresponds to a mandatory label of the first group
          bool twin_to_trpl_G2 = pairMandatoryLabels.first.count(keyVal.second->mlabel);
          if (twin_to_trpl_G2) {
            doubledG2->AddLink(keyVal.second->mid + "_TWIN_TO_TRPL",
                            keyVal.second->mupstream + "_TWIN",
                            keyVal.second->mdownstream + "_TRPL",
                            keyVal.second->mlength,
                            keyVal.second->mcosts,
                            keyVal.second->mlabel);
          }
        }

        // Set destinations as nodes of the trpl graph
        std::vector<std::string> destinationsTwin;
        for (int i=0; i < destinations.size(); i++) {
          destinationsTwin.push_back(destinations[i] + "_TRPL");
        }

        // Launch dijkstra algo for each OD in parallel
        omp_set_num_threads(threadNumber);

        int nbOD = origins.size();
        std::vector<std::vector<pathCost>> res(nbOD);
        OrientedGraph *privateDoubledG1;
        OrientedGraph *privateDoubledG2;

        #pragma omp parallel shared(res, vecAvailableLabels, vecMapLabelCosts, origins, destinationsTwin, kPaths, doubledG1, doubledG2) private(privateDoubledG1, privateDoubledG2)
        {
          privateDoubledG1 = copyGraph(*doubledG1);
          privateDoubledG2 = copyGraph(*doubledG2);

          #pragma omp for
          for (int i = 0; i < nbOD; i++)
          {
            std::vector<pathCost> resPath1;
            std::vector<pathCost> resPath2;

            if (vecAvailableLabels.empty())
            {
                // Look for shortest paths on G1
                resPath1 = KShortestPath(*privateDoubledG1, origins[i], destinationsTwin[i], cost, {}, vecMapLabelCosts[i], maxDiffCost, maxDistInCommon, costMultiplier, maxRetry, kPaths[i], true);
                // Look for shortest paths on G2
                resPath2 = KShortestPath(*privateDoubledG2, origins[i], destinationsTwin[i], cost, {}, vecMapLabelCosts[i], maxDiffCost, maxDistInCommon, costMultiplier, maxRetry, kPaths[i], true);
            }
            else
            {
                resPath1 = KShortestPath(*privateDoubledG1, origins[i], destinationsTwin[i], cost, vecAvailableLabels[i], vecMapLabelCosts[i], maxDiffCost, maxDistInCommon, costMultiplier, maxRetry, kPaths[i], true);
                resPath2 = KShortestPath(*privateDoubledG2, origins[i], destinationsTwin[i], cost, vecAvailableLabels[i], vecMapLabelCosts[i], maxDiffCost, maxDistInCommon, costMultiplier, maxRetry, kPaths[i], true);
            }
            // Concat resPath1 and resPath2
            resPath1.insert(resPath1.end(), resPath2.begin(), resPath2.end());
            // Decode paths
            for (int j = 0; j < resPath1.size(); j++){
              if (resPath1[j].first.size() > 0) {
                for (int k = 0; k < resPath1[j].first.size(); k++) {
                  if (resPath1[j].first[k].size() > 5 && (resPath1[j].first[k].compare(resPath1[j].first[k].size() - 5, 5, "_TWIN") == 0 || resPath1[j].first[k].compare(resPath1[j].first[k].size() - 5, 5, "_TRPL") == 0)) {
                    resPath1[j].first[k] = resPath1[j].first[k].substr(0, resPath1[j].first[k].size() - 5);
                  }
                }
              }
            }
            // Keep only unique paths
            sort( resPath1.begin(), resPath1.end() );
            resPath1.erase(std::unique( resPath1.begin(), resPath1.end() ), resPath1.end() );

            // Keep the k best paths found
            std::sort(resPath1.begin(), resPath1.end(), [](pathCost a, pathCost b)
              { return a.second < b.second; });
            if (resPath1.size() >= kPaths[i])
            {
                std::vector<pathCost> resPath(resPath1.begin(), resPath1.begin() + kPaths[i]);
                res[i] = resPath;
            }
            else
            {
                res[i] = resPath1;
            }


        }

        #pragma omp critical
        {
            if (privateDoubledG1) delete(privateDoubledG1);
            if (privateDoubledG2) delete(privateDoubledG2);
        }
      }
      if (doubledG1) delete(doubledG1);
      if (doubledG2) delete(doubledG2);
      return res;
    }

} // namespace hipop
