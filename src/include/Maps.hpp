#pragma once

#include <vector>
#include <array>
#include <queue>
#include <limits>
#include <cmath>
#include <algorithm>

class Maps
{
public:
    Maps(
        int grid_width,
        int grid_height,
        const std::vector<std::vector<double>> &distances_costs_map,
        const std::vector<std::vector<int>> &poi_map,
        const std::vector<std::vector<double>> &demand_map,
        const std::vector<std::vector<double>> &land_rental_cost_map);

    int get_width() const { return grid_width; }
    int get_height() const { return grid_height; }

    int get_poi_at(int i, int j) const { return poi_map[i][j]; }
    double get_demand_at(int i, int j) const { return demand_map[i][j]; }
    double get_rental_cost_at(int i, int j) const { return land_rental_cost_map[i][j]; }
    double get_distance(int m, int n, int i, int j) const { return distances_from_to[m][n][i][j]; }

    bool isPassable(int i, int j) const { return distances_costs_map[i][j] != INF_VAL; }

private:
    static constexpr double INF_VAL = -1.0;
    static constexpr std::array<std::pair<int, int>, 4> DIRECTIONS = {{ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }};

    double calculate_distance(int m, int n, int i, int j) const;
    
    std::pair<double, std::pair<int, int>> find_nearest_valid(int start_x, int start_y) const;

    int grid_width;
    int grid_height;

    std::vector<std::vector<double>> distances_costs_map;
    std::vector<std::vector<int>> poi_map;
    std::vector<std::vector<double>> demand_map;
    std::vector<std::vector<double>> land_rental_cost_map;

    std::vector<std::vector<std::vector<std::vector<double>>>> distances_from_to;
};