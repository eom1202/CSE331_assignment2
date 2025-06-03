#include "common.h"
#include <ctime>
#include <iomanip>

vector<Point> cities;
int n;
string dataset_name;
string edge_weight_type;

vector<pair<int, int>> buildMST(int root) {
    vector<bool> inMST(n, false);
    vector<double> minEdge(n, 1e9);
    vector<int> parent(n, -1);
    
    minEdge[root] = 0;
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
                double weight = distance(u, v);
                if (weight < minEdge[v]) {
                    minEdge[v] = weight;
                    parent[v] = u;
                }
            }
        }
    }
    
    return mstEdges;
}

vector<int> preorderWalk(vector<pair<int, int>>& mstEdges, int root) {
    vector<vector<int>> children(n);
    for (auto edge : mstEdges) {
        children[edge.first].push_back(edge.second);
    }
    
    vector<int> visitedOrder;
    vector<int> stack;
    stack.push_back(root);
    
    while (!stack.empty()) {
        int current = stack.back();
        stack.pop_back();
        visitedOrder.push_back(current);
        
        for (int i = children[current].size() - 1; i >= 0; i--) {
            stack.push_back(children[current][i]);
        }
    }
    
    return visitedOrder;
}

double calculateMSTCost(vector<pair<int, int>>& mstEdges) {
    double mstCost = 0.0;
    for (auto edge : mstEdges) {
        mstCost += distance(edge.first, edge.second);
    }
    return mstCost;
}

void saveResult(double tourCost, double mstCost, double executionTime) {
    filesystem::create_directory("results");
    
    ofstream file("results/results_mst.csv", ios::app);
    
    file.seekp(0, ios::end);
    if (file.tellp() == 0) {
        file << "Dataset,Cities,Tour_Cost,MST_Cost,Approximation_Ratio,Execution_Time\n";
    }
    
    file << fixed << setprecision(4);
    file << dataset_name << "," << n << "," << tourCost << "," << mstCost << "," 
         << (tourCost / mstCost) << "," << executionTime << "\n";
    
    file.close();
}

int main() {
    vector<string> files = {"dataset/ulysses16.tsp", "dataset/a280.tsp", "dataset/xql662.tsp", "dataset/kz9976.tsp", "dataset/mona-lisa100K.tsp"};
    
    for (string filename : files) {
        cout << "Processing " << filename << "..." << endl;
        
        if (!readTSPFile(filename)) {
            cout << "Failed to read " << filename << endl;
            continue;
        }
        
        clock_t start = clock();
        
        vector<pair<int, int>> mstEdges = buildMST(0);
        double mstCost = calculateMSTCost(mstEdges);
        
        vector<int> tour = preorderWalk(mstEdges, 0);
        tour.push_back(tour[0]);
        double tourCost = calculateTourCost(tour);
        
        clock_t end = clock();
        double time = double(end - start) / CLOCKS_PER_SEC;
        
        cout << fixed << setprecision(2);
        cout << dataset_name << ": " << tourCost << " (ratio: " << (tourCost/mstCost) << ")" << endl;
        
        saveResult(tourCost, mstCost, time);
    }
    
    cout << "Done!" << endl;
    return 0;
}