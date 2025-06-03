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

bool readTSPFile(string filename) {
    ifstream file(filename);
    if (!file.is_open()) return false;
    
    string line;
    cities.clear();
    edge_weight_type = "EUC_2D";
    
    while (getline(file, line)) {
        if (line.find("NAME") == 0) {
            stringstream ss(line);
            string temp, name;
            getline(ss, temp, ':');
            getline(ss, name);
            dataset_name = name;
            while (!dataset_name.empty() && dataset_name[0] == ' ') dataset_name.erase(0, 1);
        } 
        else if (line.find("DIMENSION") == 0) {
            stringstream ss(line);
            string temp;
            getline(ss, temp, ':');
            ss >> n;
        }
        else if (line.find("EDGE_WEIGHT_TYPE") == 0) {
            stringstream ss(line);
            string temp;
            getline(ss, temp, ':');
            ss >> edge_weight_type;
        }
        else if (line.find("NODE_COORD_SECTION") == 0) {
            break;
        }
    }
    
    while (getline(file, line) && line != "EOF") {
        if (!line.empty()) {
            while (!line.empty() && (line[0] == ' ' || line[0] == '\t')) {
                line.erase(0, 1);
            }
            
            if (!line.empty()) {
                stringstream ss(line);
                int id;
                double x, y;
                ss >> id >> x >> y;
                cities.push_back({x, y});
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