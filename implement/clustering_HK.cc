#include "common.h"
#include <ctime>
#include <iomanip>
#include <map>
#include <climits>
#include <random>
#include <algorithm>

vector<Point> cities;
int n;
string dataset_name;
string edge_weight_type;

map<pair<int, int>, double> clusterMemo;
vector<int> currentCluster;

struct Cluster {
    vector<int> cityIndices;
    Point centroid;
};

double distancePoints(Point a, Point b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

double heldKarpRecursive(int mask, int pos) {
    int clusterSize = currentCluster.size();
    if (mask == (1 << clusterSize) - 1) {
        return distance(currentCluster[pos], currentCluster[0]);
    }
    
    pair<int, int> state = {mask, pos};
    if (clusterMemo.find(state) != clusterMemo.end()) {
        return clusterMemo[state];
    }
    
    double result = 1e9;
    for (int next = 0; next < clusterSize; next++) {
        if (mask & (1 << next)) continue;
        
        double cost = distance(currentCluster[pos], currentCluster[next]) + heldKarpRecursive(mask | (1 << next), next);
        result = min(result, cost);
    }
    
    clusterMemo[state] = result;
    return result;
}

double heldKarpCluster(vector<int>& cluster, vector<int>& tour) {
    int clusterSize = cluster.size();
    if (clusterSize == 1) {
        tour = {cluster[0], cluster[0]};
        return 0.0;
    }
    
    clusterMemo.clear();
    currentCluster = cluster;
    
    double bestCost = heldKarpRecursive(1, 0);
    
    tour.clear();
    int mask = 1;
    int pos = 0;
    tour.push_back(cluster[0]);
    
    while (mask != (1 << clusterSize) - 1) {
        int nextPos = -1;
        double minCost = 1e9;
        
        for (int next = 0; next < clusterSize; next++) {
            if (mask & (1 << next)) continue;
            
            double cost = distance(cluster[pos], cluster[next]) + heldKarpRecursive(mask | (1 << next), next);
            if (cost < minCost) {
                minCost = cost;
                nextPos = next;
            }
        }
        
        tour.push_back(cluster[nextPos]);
        mask |= (1 << nextPos);
        pos = nextPos;
    }
    
    tour.push_back(cluster[0]);
    return bestCost;
}

vector<Cluster> kMeansPlusPlus(int k) {
    vector<Cluster> clusters(k);
    vector<Point> centroids;
    
    random_device rd;
    mt19937 gen(12345);
    
    centroids.push_back(cities[uniform_int_distribution<int>(0, n-1)(gen)]);
    
    for (int i = 1; i < k; i++) {
        vector<double> distances(n);
        double totalDist = 0.0;
        
        for (int j = 0; j < n; j++) {
            double minDist = 1e9;
            for (Point& centroid : centroids) {
                minDist = min(minDist, distancePoints(cities[j], centroid));
            }
            distances[j] = minDist * minDist;
            totalDist += distances[j];
        }
        
        double r = uniform_real_distribution<double>(0.0, totalDist)(gen);
        double cumSum = 0.0;
        int selected = 0;
        
        for (int j = 0; j < n; j++) {
            cumSum += distances[j];
            if (cumSum >= r) {
                selected = j;
                break;
            }
        }
        
        centroids.push_back(cities[selected]);
    }
    
    for (int iter = 0; iter < 100; iter++) {
        for (int i = 0; i < k; i++) {
            clusters[i].cityIndices.clear();
        }
        
        for (int i = 0; i < n; i++) {
            int bestCluster = 0;
            double minDist = distancePoints(cities[i], centroids[0]);
            
            for (int j = 1; j < k; j++) {
                double dist = distancePoints(cities[i], centroids[j]);
                if (dist < minDist) {
                    minDist = dist;
                    bestCluster = j;
                }
            }
            
            clusters[bestCluster].cityIndices.push_back(i);
        }
        
        bool converged = true;
        for (int i = 0; i < k; i++) {
            if (clusters[i].cityIndices.empty()) continue;
            
            Point newCentroid = {0.0, 0.0};
            for (int idx : clusters[i].cityIndices) {
                newCentroid.x += cities[idx].x;
                newCentroid.y += cities[idx].y;
            }
            newCentroid.x /= clusters[i].cityIndices.size();
            newCentroid.y /= clusters[i].cityIndices.size();
            
            if (distancePoints(centroids[i], newCentroid) > 0.001) {
                converged = false;
            }
            centroids[i] = newCentroid;
            clusters[i].centroid = newCentroid;
        }
        
        if (converged) break;
    }
    
    vector<Cluster> validClusters;
    for (int i = 0; i < k; i++) {
        if (!clusters[i].cityIndices.empty()) {
            validClusters.push_back(clusters[i]);
        }
    }
    
    for (int i = 0; i < validClusters.size(); i++) {
        while (validClusters[i].cityIndices.size() > 20) {
            vector<int>& clusterCities = validClusters[i].cityIndices;
            int clusterSize = clusterCities.size();
            int splitSize = clusterSize / 2;
            
            vector<vector<double>> distances(clusterSize, vector<double>(clusterSize));
            for (int j = 0; j < clusterSize; j++) {
                for (int k = 0; k < clusterSize; k++) {
                    distances[j][k] = distance(clusterCities[j], clusterCities[k]);
                }
            }
            
            vector<bool> assigned(clusterSize, false);
            vector<int> group1, group2;
            
            int seed1 = 0, seed2 = 0;
            double maxDist = 0;
            for (int j = 0; j < clusterSize; j++) {
                for (int k = j + 1; k < clusterSize; k++) {
                    if (distances[j][k] > maxDist) {
                        maxDist = distances[j][k];
                        seed1 = j;
                        seed2 = k;
                    }
                }
            }
            
            group1.push_back(clusterCities[seed1]);
            group2.push_back(clusterCities[seed2]);
            assigned[seed1] = assigned[seed2] = true;
            
            for (int remaining = clusterSize - 2; remaining > 0; remaining--) {
                if (group1.size() >= splitSize) {
                    for (int j = 0; j < clusterSize; j++) {
                        if (!assigned[j]) {
                            group2.push_back(clusterCities[j]);
                            assigned[j] = true;
                            break;
                        }
                    }
                } else if (group2.size() >= splitSize) {
                    for (int j = 0; j < clusterSize; j++) {
                        if (!assigned[j]) {
                            group1.push_back(clusterCities[j]);
                            assigned[j] = true;
                            break;
                        }
                    }
                } else {
                    int bestCity = -1;
                    double bestScore = -1e9;
                    bool assignToGroup1 = false;
                    
                    for (int j = 0; j < clusterSize; j++) {
                        if (assigned[j]) continue;
                        
                        double dist1 = 1e9, dist2 = 1e9;
                        for (int cityIdx : group1) {
                            for (int k = 0; k < clusterSize; k++) {
                                if (clusterCities[k] == cityIdx) {
                                    dist1 = min(dist1, distances[j][k]);
                                    break;
                                }
                            }
                        }
                        for (int cityIdx : group2) {
                            for (int k = 0; k < clusterSize; k++) {
                                if (clusterCities[k] == cityIdx) {
                                    dist2 = min(dist2, distances[j][k]);
                                    break;
                                }
                            }
                        }
                        
                        double score1 = dist2 - dist1;
                        double score2 = dist1 - dist2;
                        
                        if (score1 > bestScore) {
                            bestScore = score1;
                            bestCity = j;
                            assignToGroup1 = true;
                        }
                        if (score2 > bestScore) {
                            bestScore = score2;
                            bestCity = j;
                            assignToGroup1 = false;
                        }
                    }
                    
                    if (assignToGroup1) {
                        group1.push_back(clusterCities[bestCity]);
                    } else {
                        group2.push_back(clusterCities[bestCity]);
                    }
                    assigned[bestCity] = true;
                }
            }
            
            Cluster newCluster;
            newCluster.cityIndices = group2;
            
            Point newCentroid = {0.0, 0.0};
            for (int idx : newCluster.cityIndices) {
                newCentroid.x += cities[idx].x;
                newCentroid.y += cities[idx].y;
            }
            newCentroid.x /= newCluster.cityIndices.size();
            newCentroid.y /= newCluster.cityIndices.size();
            newCluster.centroid = newCentroid;
            
            validClusters[i].cityIndices = group1;
            Point updatedCentroid = {0.0, 0.0};
            for (int idx : group1) {
                updatedCentroid.x += cities[idx].x;
                updatedCentroid.y += cities[idx].y;
            }
            updatedCentroid.x /= group1.size();
            updatedCentroid.y /= group1.size();
            validClusters[i].centroid = updatedCentroid;
            
            validClusters.push_back(newCluster);
        }
    }
    
    return validClusters;
}

vector<int> getClusterOrder(vector<Cluster>& clusters) {
    int numClusters = clusters.size();
    if (numClusters == 1) return {0};
    
    vector<Point> clusterPoints;
    for (int i = 0; i < numClusters; i++) {
        clusterPoints.push_back(clusters[i].centroid);
    }
    
    vector<Point> originalCities = cities;
    int originalN = n;
    cities = clusterPoints;
    n = numClusters;
    
    vector<bool> inMST(n, false);
    vector<double> minEdge(n, 1e9);
    vector<int> parent(n, -1);
    
    minEdge[0] = 0;
    vector<pair<int, int>> mstEdges;
    
    for (int step = 0; step < n; step++) {
        int u = -1;
        for (int v = 0; v < n; v++) {
            if (!inMST[v] && (u == -1 || minEdge[v] < minEdge[u])) {
                u = v;
            }
        }
        
        inMST[u] = true;
        
        if (parent[u] != -1) {
            mstEdges.push_back({parent[u], u});
        }
        
        for (int v = 0; v < n; v++) {
            if (!inMST[v]) {
                double weight = distancePoints(clusterPoints[u], clusterPoints[v]);
                if (weight < minEdge[v]) {
                    minEdge[v] = weight;
                    parent[v] = u;
                }
            }
        }
    }
    
    vector<vector<int>> children(n);
    for (auto edge : mstEdges) {
        children[edge.first].push_back(edge.second);
    }
    
    vector<int> visitedOrder;
    vector<int> stack;
    stack.push_back(0);
    
    while (!stack.empty()) {
        int current = stack.back();
        stack.pop_back();
        visitedOrder.push_back(current);
        
        for (int i = children[current].size() - 1; i >= 0; i--) {
            stack.push_back(children[current][i]);
        }
    }
    
    cities = originalCities;
    n = originalN;
    
    return visitedOrder;
}

vector<int> mergeCycles(vector<int>& cycleA, vector<int>& cycleB) {
    double minCost = 1e9;
    vector<int> bestMergedCycle;
    
    for (int nodeA : cycleA) {
        if (nodeA == cycleA.back()) continue;
        
        for (int nodeB : cycleB) {
            if (nodeB == cycleB.back()) continue;
            
            double connectionCost = distance(nodeA, nodeB);
            
            int posA = 0;
            for (int i = 0; i < cycleA.size() - 1; i++) {
                if (cycleA[i] == nodeA) {
                    posA = i;
                    break;
                }
            }
            
            int posB = 0;
            for (int i = 0; i < cycleB.size() - 1; i++) {
                if (cycleB[i] == nodeB) {
                    posB = i;
                    break;
                }
            }
            
            int prevA = (posA - 1 + cycleA.size() - 1) % (cycleA.size() - 1);
            int nextA = (posA + 1) % (cycleA.size() - 1);
            int prevB = (posB - 1 + cycleB.size() - 1) % (cycleB.size() - 1);
            int nextB = (posB + 1) % (cycleB.size() - 1);
            
            double removedCostA = distance(cycleA[prevA], nodeA) + distance(nodeA, cycleA[nextA]);
            double removedCostB = distance(cycleB[prevB], nodeB) + distance(nodeB, cycleB[nextB]);
            
            double addedCostA = distance(cycleA[prevA], cycleA[nextA]);
            double addedCostB = distance(cycleB[prevB], cycleB[nextB]);
            
            double totalCost = connectionCost + addedCostA + addedCostB - removedCostA - removedCostB;
            
            if (totalCost < minCost) {
                minCost = totalCost;
                
                bestMergedCycle.clear();
                
                for (int i = 0; i < cycleA.size() - 1; i++) {
                    if (cycleA[i] != nodeA) {
                        bestMergedCycle.push_back(cycleA[i]);
                    }
                }
                
                for (int i = 0; i < cycleB.size() - 1; i++) {
                    if (cycleB[i] != nodeB) {
                        bestMergedCycle.push_back(cycleB[i]);
                    }
                }
                
                bestMergedCycle.push_back(bestMergedCycle[0]);
            }
        }
    }
    
    return bestMergedCycle;
}

void saveResult(double tourCost, double executionTime) {
    filesystem::create_directory("results");
    
    ofstream file("results/results_clustering.csv", ios::app);
    
    file.seekp(0, ios::end);
    if (file.tellp() == 0) {
        file << "Dataset,Cities,Tour_Cost,Execution_Time\n";
    }
    
    file << fixed << setprecision(4);
    file << dataset_name << "," << n << "," << tourCost << "," << executionTime << "\n";
    
    file.close();
}

int main() {
    vector<string> files = {"dataset/ulysses16.tsp","dataset/a280.tsp","dataset/xql662","dataset/kz9976.tsp"};
    
    for (string filename : files) {
        cout << "Processing " << filename << "..." << endl;
        
        if (!readTSPFile(filename)) {
            cout << "Failed to read " << filename << endl;
            continue;
        }
        
        clock_t start = clock();
        
        int k = max(1, n / 16);
        vector<Cluster> clusters = kMeansPlusPlus(k);
        
        cout << "Clustering completed" << endl;
        
        vector<vector<int>> clusterTours(clusters.size());
        
        for (int i = 0; i < clusters.size(); i++) {
            cout << "Solving cluster " << (i+1) << "/" << clusters.size() 
                 << " (size: " << clusters[i].cityIndices.size() << ")" << endl;
            
            vector<int> tour;
            heldKarpCluster(clusters[i].cityIndices, tour);
            clusterTours[i] = tour;
        }
        
        vector<int> clusterOrder = getClusterOrder(clusters);
        
        vector<int> finalTour = clusterTours[clusterOrder[0]];
        
        for (int i = 1; i < clusterOrder.size(); i++) {
            int currentClusterIdx = clusterOrder[i];
            finalTour = mergeCycles(finalTour, clusterTours[currentClusterIdx]);
        }
        
        double tourCost = calculateTourCost(finalTour);
        
        clock_t end = clock();
        double time = double(end - start) / CLOCKS_PER_SEC;
        
        cout << fixed << setprecision(2);
        cout << dataset_name << ": " << tourCost << " (time: " << time << "s)" << endl;
        
        saveResult(tourCost, time);
    }
    
    cout << "Done!" << endl;
    return 0;
}