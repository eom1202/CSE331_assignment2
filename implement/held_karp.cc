#include "common.h"
#include <ctime>
#include <iomanip>
#include <map>
#include <climits>

vector<Point> cities;
int n;
string dataset_name;
string edge_weight_type;
map<pair<int, int>, double> memo;

double heldKarp(int mask, int pos) {
    if (mask == (1 << n) - 1) {
        return distance(pos, 0);
    }
    
    pair<int, int> state = {mask, pos};
    if (memo.find(state) != memo.end()) {
        return memo[state];
    }
    
    double result = 1e9;
    for (int next = 0; next < n; next++) {
        if (mask & (1 << next)) continue;
        
        double cost = distance(pos, next) + heldKarp(mask | (1 << next), next);
        result = min(result, cost);
    }
    
    memo[state] = result;
    return result;
}

vector<int> reconstructPath(int start) {
    vector<int> path;
    int mask = 1 << start;
    int pos = start;
    path.push_back(start);
    
    while (mask != (1 << n) - 1) {
        int nextCity = -1;
        double minCost = 1e9;
        
        for (int next = 0; next < n; next++) {
            if (mask & (1 << next)) continue;
            
            double cost = distance(pos, next) + heldKarp(mask | (1 << next), next);
            if (cost < minCost) {
                minCost = cost;
                nextCity = next;
            }
        }
        
        path.push_back(nextCity);
        mask |= (1 << nextCity);
        pos = nextCity;
    }
    
    path.push_back(start);
    return path;
}

void saveResult(double tourCost, double executionTime) {
    filesystem::create_directory("results");
    
    ofstream file("results/results_held_karp.csv", ios::app);
    
    file.seekp(0, ios::end);
    if (file.tellp() == 0) {
        file << "Dataset,Cities,Tour_Cost,Execution_Time\n";
    }
    
    file << fixed << setprecision(4);
    file << dataset_name << "," << n << "," << tourCost << "," << executionTime << "\n";
    
    file.close();
}

int main() {
    vector<string> files = {"dataset/ulysses16.tsp"};
    
    for (string filename : files) {
        cout << "Processing " << filename << "..." << endl;
        
        if (!readTSPFile(filename)) {
            cout << "Failed to read " << filename << endl;
            continue;
        }
        
        if (n > 20) {
            cout << "Skipping " << dataset_name << " (too large: " << n << " cities)" << endl;
            continue;
        }
        
        memo.clear();
        
        clock_t start = clock();
        
        double tourCost = heldKarp(1, 0);
        vector<int> tour = reconstructPath(0);
        
        clock_t end = clock();
        double time = double(end - start) / CLOCKS_PER_SEC;
        
        cout << fixed << setprecision(2);
        cout << dataset_name << ": " << tourCost << " (time: " << time << "s)" << endl;
        
        cout << "Tour order: ";
        for (int i = 0; i < tour.size(); i++) {
            cout << (tour[i] + 1);
            if (i < tour.size() - 1) cout << ", ";
        }
        cout << endl;        
        saveResult(tourCost, time);
    }
    
    cout << "Done!" << endl;
    return 0;
}