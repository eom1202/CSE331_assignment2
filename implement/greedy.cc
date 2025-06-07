#include "common.h"
#include <ctime>
#include <iomanip>

vector<Point> cities;
int n;
string dataset_name;
string edge_weight_type;

vector<int> greedyTSP(int start) {
    vector<bool> visited(n, false);
    vector<int> tour;
    
    int current = start;
    visited[current] = true;
    tour.push_back(current);
    
    for (int step = 1; step < n; step++) {
        int nearest = -1;
        double minDist = 1e9;
        
        for (int next = 0; next < n; next++) {
            if (!visited[next]) {
                double dist = distance(current, next);
                if (dist < minDist) {
                    minDist = dist;
                    nearest = next;
                }
            }
        }
        
        visited[nearest] = true;
        tour.push_back(nearest);
        current = nearest;
    }
    
    tour.push_back(start);
    return tour;
}

void saveResult(double tourCost, double executionTime) {
    filesystem::create_directory("results");
    
    ofstream file("results/results_greedy.csv", ios::app);
    
    file.seekp(0, ios::end);
    if (file.tellp() == 0) {
        file << "Dataset,Cities,Tour_Cost,Execution_Time\n";
    }
    
    file << fixed << setprecision(4);
    file << dataset_name << "," << n << "," << tourCost << "," << executionTime << "\n";
    
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
        
        vector<int> tour = greedyTSP(0);
        double tourCost = calculateTourCost(tour);
        
        clock_t end = clock();
        double time = double(end - start) / CLOCKS_PER_SEC;
        
        cout << fixed << setprecision(2);
        cout << dataset_name << ": " << tourCost << " (time: " << time << "s)" << endl;
        
        saveResult(tourCost, time);
    }
    
    cout << "Done!" << endl;
    return 0;
}