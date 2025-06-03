#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <filesystem>

using namespace std;

struct Point {
    double x, y;
};

vector<Point> cities;
int n;
string dataset_name;

double distance(int i, int j) {
    if (i == j) return 0.0;
    double dx = cities[i].x - cities[j].x;
    double dy = cities[i].y - cities[j].y;
    return sqrt(dx * dx + dy * dy);
}

bool readTSPFile(string filename) {
    ifstream file(filename);
    if (!file.is_open()) return false;
    
    string line;
    cities.clear();
    
    while (getline(file, line)) {
        if (line.find("NAME") == 0) {
            stringstream ss(line);
            string temp, name;
            getline(ss, temp, ':');
            getline(ss, name);
            dataset_name = name;
            while (dataset_name[0] == ' ') dataset_name.erase(0, 1);
        } 
        else if (line.find("DIMENSION") == 0) {
            stringstream ss(line);
            string temp;
            getline(ss, temp, ':');
            ss >> n;
        } 
        else if (line.find("NODE_COORD_SECTION") == 0) {
            break;
        }
    }
    
    while (getline(file, line) && line != "EOF") {
        if (!line.empty()) {
            stringstream ss(line);
            int id;
            double x, y;
            ss >> id >> x >> y;
            cities.push_back({x, y});
        }
    }
    
    return cities.size() == n;
}

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

double calculateTourCost(vector<int>& tourPath) {
    double totalCost = 0.0;
    for (int i = 0; i < tourPath.size() - 1; i++) {
        totalCost += distance(tourPath[i], tourPath[i + 1]);
    }
    return totalCost;
}

double calculateMSTCost(vector<pair<int, int>>& mstEdges) {
    double mstCost = 0.0;
    for (auto edge : mstEdges) {
        mstCost += distance(edge.first, edge.second);
    }
    return mstCost;
}

void saveResult(double tourCost, double mstCost, double executionTime) {
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
    vector<string> files = {"dataset/a280.tsp", "dataset/xql662.tsp", "dataset/kz9976.tsp", "dataset/mona-lisa100K.tsp"};
    
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
    
    cout << "Done! Results saved to results.csv" << endl;
    return 0;
}