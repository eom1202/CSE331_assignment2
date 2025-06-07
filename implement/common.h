#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <filesystem>

using namespace std;

struct Point {
    double x, y;
};

extern vector<Point> cities;
extern int n;
extern string dataset_name;
extern string edge_weight_type;

double distance(int i, int j) {
    if (i == j) return 0.0;
    
    if (edge_weight_type == "GEO") {
        double PI = 3.141592653589793;
        double RRR = 6378.388;
        
        double lat1_deg = (int)cities[i].x + (cities[i].x - (int)cities[i].x) * 100.0 / 60.0;
        double lon1_deg = (int)cities[i].y + (cities[i].y - (int)cities[i].y) * 100.0 / 60.0;
        double lat2_deg = (int)cities[j].x + (cities[j].x - (int)cities[j].x) * 100.0 / 60.0;
        double lon2_deg = (int)cities[j].y + (cities[j].y - (int)cities[j].y) * 100.0 / 60.0;
        
        double lat1 = lat1_deg * PI / 180.0;
        double lon1 = lon1_deg * PI / 180.0;
        double lat2 = lat2_deg * PI / 180.0;
        double lon2 = lon2_deg * PI / 180.0;
        
        double q1 = cos(lon1 - lon2);
        double q2 = cos(lat1 - lat2);
        double q3 = cos(lat1 + lat2);
        
        return (int)(RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    } else {
        double dx = cities[i].x - cities[j].x;
        double dy = cities[i].y - cities[j].y;
        return sqrt(dx * dx + dy * dy);
    }
}

string trim(string str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == string::npos) return "";
    size_t end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

bool readTSPFile(string filename) {
    ifstream file(filename);
    if (!file.is_open()) return false;
    
    string line;
    cities.clear();
    edge_weight_type = "EUC_2D";
    dataset_name = "";
    
    while (getline(file, line)) {
        if (line.find("NAME") == 0) {
            size_t colonPos = line.find(':');
            if (colonPos != string::npos) {
                dataset_name = trim(line.substr(colonPos + 1));
            }
        } 
        else if (line.find("DIMENSION") == 0) {
            size_t colonPos = line.find(':');
            if (colonPos != string::npos) {
                string numStr = line.substr(colonPos + 1);
                stringstream ss(numStr);
                ss >> n;
            }
        }
        else if (line.find("EDGE_WEIGHT_TYPE") == 0) {
            size_t colonPos = line.find(':');
            if (colonPos != string::npos) {
                string typeStr = line.substr(colonPos + 1);
                stringstream ss(typeStr);
                ss >> edge_weight_type;
            }
        }
        else if (line.find("NODE_COORD_SECTION") == 0) {
            break;
        }
    }
    
    while (getline(file, line) && line.find("EOF") == string::npos) {
        if (!line.empty()) {
            line = trim(line);
            
            if (!line.empty()) {
                stringstream ss(line);
                int id;
                double x, y;
                if (ss >> id >> x >> y) {
                    cities.push_back({x, y});
                }
            }
        }
    }
    
    return cities.size() == n;
}

double calculateTourCost(vector<int>& tourPath) {
    double totalCost = 0.0;
    for (int i = 0; i < tourPath.size() - 1; i++) {
        totalCost += distance(tourPath[i], tourPath[i + 1]);
    }
    return totalCost;
}

#endif