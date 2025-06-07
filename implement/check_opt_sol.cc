#include "common.h"

vector<Point> cities;
int n;
string dataset_name;
string edge_weight_type;

int main() {
    vector<string> files = {"a280", "xql662"};
    
    for (string name : files) {
        if (!readTSPFile("dataset/" + name + ".tsp")) {
            continue;
        }
        
        ifstream tourFile("dataset_opt/" + name + ".opt.tour");
        vector<int> tour;
        string line;
        bool reading = false;
        
        while (getline(tourFile, line)) {
            if (line.find("TOUR_SECTION") != string::npos) {
                reading = true;
                continue;
            }
            if (line == "-1" || line == "EOF") break;
            
            if (reading && !line.empty()) {
                string trimmed = trim(line);
                if (!trimmed.empty() && trimmed != "EOF") {
                    try {
                        int city = stoi(trimmed);
                        if (city > 0) {
                            tour.push_back(city - 1);
                        }
                    } catch (...) {
                        continue;
                    }
                }
            }
        }
        
        if (!tour.empty()) {
            tour.push_back(tour[0]);
        }
        
        double total = 0;
        double totalInt = 0;
        cout << "Tour has " << tour.size() << " cities (including return)" << endl;
        
        for (int i = 0; i < tour.size() - 1; i++) {
            double dist = distance(tour[i], tour[i + 1]);
            int intDist = (int)(dist + 0.5);
            total += dist;
            totalInt += intDist;
            if (i < 5) {
                cout << "City " << tour[i]+1 << " to " << tour[i+1]+1 << ": " << dist << " (int: " << intDist << ")" << endl;
            }
        }
        
        cout << name << " optimal tour length (double): " << total << endl;
        cout << name << " optimal tour length (int): " << totalInt << endl;
    }
    
    return 0;
}