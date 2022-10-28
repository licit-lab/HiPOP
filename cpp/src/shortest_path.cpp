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

        for (size_t i = 0; i < path.size() - 1; i++)
        {
            Link *link = G.mnodes[path[i]]->madj[path[i + 1]];
            c += link->mcosts[mapLabelCost.at(link->mlabel)][cost];
        }

        return c;
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
    
} // namespace hipop



