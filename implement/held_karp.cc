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
clock_t startTime;
const double TIME_LIMIT = 1200.0;

bool isTimeExceeded() {
    clock_t current = clock();
    double elapsed = double(current - startTime) / CLOCKS_PER_SEC;
    return elapsed > TIME_LIMIT;
}

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
        
        if (isTimeExceeded()) {
            break;
        }
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
        if (isTimeExceeded()) {
            break;
        }
        
        int nextCity = -1;
        double minCost = 1e9;
        
        for (int next = 0; next < n; next++) {
            if (mask & (1 << next)) continue;
            
            double cost = distance(pos, next) + heldKarp(mask | (1 << next), next);
            if (cost < minCost) {
                minCost = cost;
                nextCity = next;
            }
            
            if (isTimeExceeded()) {
                break;
            }
        }
        
        if (nextCity == -1) break;
        
        path.push_back(nextCity);
        mask |= (1 << nextCity);
        pos = nextCity;
    }
    
    path.push_back(start);
    return path;
}

void saveResult(double tourCost, double executionTime, vector<int>& tour, bool timeout = false) {
    filesystem::create_directory("results");
    
    ofstream file("results/results_held_karp.csv", ios::app);
    
    file.seekp(0, ios::end);
    if (file.tellp() == 0) {
        file << "Dataset,Cities,Tour_Cost,Execution_Time,Tour_Order,Status\n";
    }
    
    string tourOrder = "";
    if (!timeout && tour.size() > 0) {
        for (int i = 0; i < tour.size(); i++) {
            tourOrder += to_string(tour[i] + 1);
            if (i < tour.size() - 1) tourOrder += " ";
        }
    } else {
        tourOrder = "TIMEOUT";
    }
    
    file << fixed << setprecision(4);
    file << dataset_name << "," << n << "," << tourCost << "," << executionTime << "," << tourOrder;
    
    if (timeout) {
        file << ",TIMEOUT";
    } else {
        file << ",COMPLETED";
    }
    
    file << "\n";
    
    file.close();
}

int main() {
    vector<string> files = {"dataset/ulysses16.tsp", "dataset/a280.tsp"};
    
    for (string filename : files) {
        cout << "Processing " << filename << "..." << endl;
        
        if (!readTSPFile(filename)) {
            cout << "Failed to read " << filename << endl;
            continue;
        }
        
        memo.clear();
        
        startTime = clock();
        
        double tourCost = heldKarp(1, 0);
        
        clock_t end = clock();
        double time = double(end - startTime) / CLOCKS_PER_SEC;
        
        if (time > TIME_LIMIT || tourCost >= 1e9) {
            cout << fixed << setprecision(2);
            cout << dataset_name << ": TIMEOUT after " << time << "s" << endl;
            
            vector<int> emptyTour;
            saveResult(-1, time, emptyTour, true);
        } else {
            vector<int> tour = reconstructPath(0);
            
            cout << fixed << setprecision(2);
            cout << dataset_name << ": " << tourCost << " (time: " << time << "s)" << endl;
            
            cout << "Tour order: ";
            for (int i = 0; i < tour.size(); i++) {
                cout << (tour[i] + 1);
                if (i < tour.size() - 1) cout << ", ";
            }
            cout << endl;        
            saveResult(tourCost, time, tour);
        }
    }
    
    cout << "Done!" << endl;
    return 0;
}