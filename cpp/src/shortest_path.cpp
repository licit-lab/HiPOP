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
     */
    void increaseCostsFromPath(OrientedGraph &G, const std::vector<std::string> &path, linkMapCosts &initial_costs)
    {

        for (size_t i = 0; i < path.size() - 1; i++)
        {
            Link *link = G.mnodes[path[i]]->madj[path[i + 1]];
            if (initial_costs.find(link->mid) == initial_costs.end())
            {
                for (auto &keyMapCost : link->mcosts)
                {
                    for(auto &keyVal: keyMapCost.second) {
                        initial_costs[link->mid][keyMapCost.first][keyVal.first] = keyVal.second;
                        keyVal.second *= 10;
                    }
                }
            }
            else
            {
                for(auto &keyMapCost: link->mcosts) {
                    for(auto &keyVal: keyMapCost.second) {
                        keyVal.second *= 10;
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
     * @brief Compute K shortest path for an origin/destination. The first path is normally computed, then we increase the costs of the first path and perform a Dijkstra again.
     * The last step is performed until either we reach K accepted shortest paths or the last computed shortest path do not meet the requirement minDist <= kPathLength - firstPathLength <= maxDist
     * 
     * @param G The OrientedGraph on which we compute the paths
     * @param origin The origin
     * @param destination The destination
     * @param cost The cost to consider
     * @param accessibleLabels The set of accessible label
     * @param mapLabelCost The type of cost map to choose on each label
     * @param minDist The minimal distance difference
     * @param maxDist The maximal distance difference
     * @param kPath The number of path to compute
     * @return std::vector<pathCost> The vector of k computed paths
     */
    std::vector<pathCost> KShortestPath(
        OrientedGraph &G, 
        const std::string &origin, 
        const std::string &destination, 
        const std::string &cost, 
        setstring accessibleLabels,
        const std::unordered_map<std::string, std::string> &mapLabelCost,
        double minDist, 
        double maxDist, 
        int kPath)
    {
        std::vector<pathCost> paths;
        linkMapCosts initial_costs;

        pathCost firstPath = dijkstra(G, origin, destination, cost, mapLabelCost, accessibleLabels);
        paths.push_back(firstPath);

        if (firstPath.first.empty())
        {
            return paths;
        }

        double firstPathLength = computePathLength(G, firstPath.first);

        increaseCostsFromPath(G, firstPath.first, initial_costs);

        int pathCounter = 1, retry = 0;

        while (pathCounter < kPath && retry < 10)
        {
            pathCost newPath = dijkstra(G, origin, destination, cost, mapLabelCost, accessibleLabels);

            if (newPath.first.empty())
            {
                break;
            }

            increaseCostsFromPath(G, newPath.first, initial_costs);
            double newPathLength = computePathLength(G, newPath.first);

            double diffPathLength = newPathLength - firstPathLength;

            if (minDist <= diffPathLength && diffPathLength <= maxDist)
            {
                bool isNew = true;
                for (const auto &p : paths)
                {
                    if (p.first == newPath.first)
                    {
                        isNew = false;
                        break;
                    }
                }

                if (isNew)
                {
                    paths.push_back(newPath);
                    retry = 0;
                    pathCounter += 1;
                }
                else
                {
                    retry += 1;
                }
            }

            else
            {
                retry += 1;
            }
        }

        for (const auto &keyVal : initial_costs)
        {
            G.mlinks[keyVal.first]->mcosts = keyVal.second;
        }

        for (auto &p : paths)
        {
            p.second = computePathCost(G, p.first, cost, mapLabelCost);
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
     * @param minDist The minimal distance difference
     * @param maxDist The maximal distance difference
     * @param kPath The number of path to compute
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
        double minDist, 
        double maxDist, 
        int kPath, 
        int threadNumber)
    {
        omp_set_num_threads(threadNumber);
        int nbOD = origins.size();

        std::vector<std::vector<pathCost>> res(nbOD);
        OrientedGraph *privateG;

        #pragma omp parallel shared(res, accessibleLabels, G, vecMapLabelCosts, origins, destinations) private(privateG)
        {
            privateG = copyGraph(G);

            #pragma omp for
            for (int i = 0; i < nbOD; i++)
            {
                if (accessibleLabels.empty())
                {
                    res[i] = KShortestPath(*privateG, origins[i], destinations[i], cost, {}, vecMapLabelCosts[i], minDist, maxDist, kPath);
                }
                else
                {
                    res[i] = KShortestPath(*privateG, origins[i], destinations[i], cost, accessibleLabels[i], vecMapLabelCosts[i], minDist, maxDist, kPath);
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
     * @param vecAvailableLabels The vector of available labels
     * @param pairMandatoryLabels The pair of labels groups the shortest paths must contain
     * @return std::vector<pathCost> The vector of computed shortest path
     */
    std::vector<pathCost> parallelIntermodalDijkstra(
        const OrientedGraph &G,
        std::vector<std::string> origins,
        std::vector<std::string> destinations,
        std::vector<std::unordered_map<std::string, std::string> > vecMapLabelCosts,
        std::string cost,
        int threadNumber,
        std::pair<std::unordered_set<std::string>, std::unordered_set<std::string>> pairMandatoryLabels,
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

        int nbPath = origins.size();
        std::vector<pathCost> res(nbPath);

        #pragma omp parallel for shared(res, vecAvailableLabels, vecMapLabelCosts) schedule(dynamic)
        for (int i = 0; i < nbPath; i++)
        {
            pathCost resPath1;
            pathCost resPath2;
            if (vecAvailableLabels.empty())
            {
                // Look for shortest path on G1
                resPath1 = dijkstra(*doubledG1, origins[i], destinationsTwin[i], cost, vecMapLabelCosts[i], {});
                // Look for shortest path on G2
                resPath2 = dijkstra(*doubledG2, origins[i], destinationsTwin[i], cost, vecMapLabelCosts[i], {});
            }
            else
            {
                resPath1 = dijkstra(*doubledG1, origins[i], destinationsTwin[i], cost, vecMapLabelCosts[i], vecAvailableLabels[i]);
                resPath2 = dijkstra(*doubledG2, origins[i], destinationsTwin[i], cost, vecMapLabelCosts[i], vecAvailableLabels[i]);
            }
            // Keep best of the two paths
            pathCost resPath;
            if (resPath1.second <= resPath2.second) {
               resPath = resPath1;
            }
            else {
               resPath = resPath2;
            }
            // Decode path
            if (resPath.first.size() > 0) {
              for (int k = 0; k < resPath.first.size(); k++) {
                if (resPath.first[k].size() > 5 && (resPath.first[k].compare(resPath.first[k].size() - 5, 5, "_TWIN") == 0 || resPath.first[k].compare(resPath.first[k].size() - 5, 5, "_TRPL") == 0)) {
                  resPath.first[k] = resPath.first[k].substr(0, resPath.first[k].size() - 5);
                }
              }
            }
            res[i] = resPath;
        }

        return res;
    }

    /**
     * @brief Compute the shortest path between origin and every possible destination using the Dijkstra algorithm
     * 
     * @param G The OrientedGrah used for the shortest path
     * @param origin The origin 
     * @param cost The costs to consider in the shortest path algorithm
     * @param mapLabelCost The type of cost map to choose on each label (mulitple set of costs can be defined on a Link)
     * @param accessibleLabels The set of accessible label
     * @return std::unordered_map<std::string, pathCost> The map of shortest paths for each destination
     */
    std::unordered_map<std::string, pathCost> multiDestDijkstra(
        const OrientedGraph &G, 
        const std::string &origin, 
        const std::string &cost, 
        const std::unordered_map<std::string, std::string> &mapLabelCost, 
        setstring accessibleLabels)
    {
        int nbPath = G.mnodes.size();
        //std::vector<pathCost> res(nbPath-1);
        std::unordered_map<std::string, pathCost> res(nbPath-1);

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

        prev[origin] = "";
        // Explore full graph
        while (!pq.empty())
        {
            QueueItem current = pq.top();
            pq.pop();
            std::string u = current.second;

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
        //for (const auto keyVal: prev){
        //    std::cout << keyVal.first << ": " << keyVal.second <<"\n";
        //}
        // Compute all paths
        std::unordered_map<std::string, std::vector<std::string>> paths = restorePaths(origin, prev);
        // Extract path results
        int i = 0;
        for (const auto u : G.mnodes) //for path in paths and lists 
        {
            if (u.first != origin)
            {   
                // initialize current path
                pathCost path;
                // initialise path cost
                path.second = dist[u.first];
                // get path
                path.first = paths[u.first];
                res[u.first] = path;
                i++;
            }
        }
        return res;
    }

    /**
     * @brief Based on the exploration of the graph within multiDestDijkstra, use the predecessor array prev to reconstruct efficiently all shortest paths.
     * 
     * @param rootNode The root node or origin of the multiDestDijkstra computation
     * @param prev The predecessor array
     * @return std::unordered_map<std::string, std::vector<std::string>> The map of paths to each destination
     */
    std::unordered_map<std::string, std::vector<std::string>> restorePaths(
        const std::string &rootNode,
        std::unordered_map<std::string, std::string> &prev
    ){
        int nbPath = prev.size();
        std::unordered_map<std::string, std::vector<std::string>> paths(nbPath-1);
        paths.reserve(nbPath-1);
        for(const auto &currentNode: prev) {
            std::string leave = currentNode.first;
            if (leave!=rootNode){
                getRootPath(leave, rootNode, prev, paths);
            }
        }
        return paths;
    }

    /**
     * @brief Recursive function to retrieve the shortest path between the root node and the currently considered node
     * 
     * @param currentNode The current destination considered
     * @param rootNode The root node or origin of the multiDestDijkstra computation
     * @param prev The predecessor array
     * @param paths The map of paths to each destination
     * @return std::vector<std::string> the shortest path between rootNode and currentNode
     */
    std::vector<std::string> getRootPath(
        const std::string &currentNode, 
        const std::string &rootNode,
        std::unordered_map<std::string, std::string> &prev,
        std::unordered_map<std::string, std::vector<std::string>> &paths)
    {
        std::vector<std::string> pathCurrent;
        pathCurrent.reserve(prev.size());
        if(paths.find(currentNode) != paths.end()){
            // If path is known already
            pathCurrent = paths[currentNode];
        } else {
            //Get previous node 
            std::string previousNode = prev[currentNode];
            if (previousNode == rootNode){
                //Initialize path
                pathCurrent.push_back(rootNode);
            } else {
                std::vector<std::string> pathPrevious = getRootPath(previousNode, rootNode, prev, paths);
                //Create current path from pathPrevious
                pathCurrent = pathPrevious;
            }
            //Append current node at the end
            pathCurrent.push_back(currentNode);
            //Save path
            paths[currentNode] = pathCurrent;
        }
        return pathCurrent;
    }

    /**
     * @brief Batch computation of Dijkstra paths for all destination nodes using openmp, each thread has its own deep copy of the OrientedGraph to ensure that the increase of the cost do not collapse with the other threads
     * 
     * @param G The OrientedGraph on which we compute the paths
     * @param origins The origins
     * @param destinations The destinations
     * @param cost The cost to consider
     * @param vecMapLabelCosts The vector of type of cost map to choose on each label
     * @param accessibleLabels The vector set of accessible label
     * @param minDist The minimal distance difference
     * @param maxDist The maximal distance difference
     * @param kPath The number of path to compute
     * @param threadNumber Number of threads to use
     * @return std::vector<std::vector<pathCost>> 
     */
    std::unordered_map<std::string, std::unordered_map<std::string, pathCost>> parallelMultiDestDijkstra(
        const OrientedGraph &G, 
        const std::vector<std::string> &origins, 
        const std::string &cost,
        const std::vector<std::unordered_map<std::string, std::string> > vecMapLabelCosts,
        int threadNumber,
        const std::vector<setstring> accessibleLabels)
    {
        omp_set_num_threads(threadNumber);
        int nbOD = origins.size();

        std::unordered_map<std::string, std::unordered_map<std::string, pathCost>> res(nbOD);
        res.reserve(nbOD);

        OrientedGraph *privateG;

        #pragma omp parallel shared(res, accessibleLabels, G, vecMapLabelCosts, origins) private(privateG)
        {
            privateG = copyGraph(G);

            /*
            #pragma omp for
            for (int i = 0; i < nbOD; i++)
            {
                if (accessibleLabels.empty())
                {
                    res[origins[i]] = multiDestDijkstra(*privateG, origins[i], cost, vecMapLabelCosts[i], {});
                }
                else
                {
                    res[origins[i]] = multiDestDijkstra(*privateG, origins[i], cost, vecMapLabelCosts[i], accessibleLabels[i]);
                }
            }
            */


            
            int q = 50;
            int max_int = nbOD/q;
            #pragma omp for
            //for (int i = 0; i < nbOD; i++)
            for (int k = 0; k<=max_int; k+=1)
            {
                for (int i=k*q; i<(k+1)*q; i++){
                    if (i<nbOD){
                        if (accessibleLabels.empty())
                        {
                            res[origins[i]] = multiDestDijkstra(*privateG, origins[i], cost, vecMapLabelCosts[i], {});
                        }
                        else
                        {
                            res[origins[i]] = multiDestDijkstra(*privateG, origins[i], cost, vecMapLabelCosts[i], accessibleLabels[i]);
                        }
                    }
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


} // namespace hipop
