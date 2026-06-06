#pragma once

#include <vector>
#include <array>
#include <queue>
#include <limits>
#include <cmath>
#include <algorithm>

/**
 * @brief Represents physical, logistical and energetic infrustructure of the problem.
 * 
 */
class Maps
{
public:
    explicit Maps(
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

    /**
     * @brief Checks whether cell (i, j) can be passee through.
     * 
     * @param i Horizontal coordinate.
     * @param j Vertical coordinate.
     * @return Information, whether cell can be passed through.
     */
    bool isPassable(int i, int j) const { return distances_costs_map[i][j] != INF_VAL; }

private:
    static constexpr double INF_VAL = -1.0;
    static constexpr std::array<std::pair<int, int>, 4> DIRECTIONS = {{ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }};

    /**
     * @brief Calculates distance from (m, n) to the (i, j) cell. It uses Dijkstra algorithm
     * to calculate this distance. It avoids obstacles (fields with distance_cost = INF). It 
     * starts from the first valid point on the map. It uses `find_nearest_valid` method to
     * achieve that.
     * 
     * @param m Horizontal start coordinate.
     * @param n Vertical start coordinate.
     * @param i Horizontal end coordinate.
     * @param j Vertical end coordinate.
     * @return Length of the path betweeen two points.
     */
    double calculate_distance(int m, int n, int i, int j) const;
    
    /**
     * @brief Finds coordinates of the first valid point, from which program can start calculating 
     * distance. It uses BFS algorithm to find first cell, where distance_cost != INF.
     * 
     * @param start_x Horizontal search start coordinate.
     * @param start_y Vertical search start coordinate.
     * @return std::pair<double, std::pair<int, int>>  Coordinates of the first valid point.
     */
    std::pair<double, std::pair<int, int>> find_nearest_valid(int start_x, int start_y) const;

    int grid_width;
    int grid_height;

    /**
     * @brief Illustrate problem physical infrastructure. It contains information about
     * obstacles (INF) and the cost of travelling from one point to another.
     * Shape: [grid_width, grid_height].
     * 
     */
    std::vector<std::vector<double>> distances_costs_map;

    /**
     * @brief Represents locations of the potential stations. It's the binary vector,
     * where 1 means that station can be places and 0, that it cannot be placed there.
     * Shape: [grid_width, grid_height].
     * 
     */
    std::vector<std::vector<int>> poi_map;

    /**
     * @brief Represents map of the energetic demand from the electric vehicle users.
     * Demands should usually, be placed in the location of households of the cars owners.
     * Shape: [grid_width, grid_height].
     * 
     */
    std::vector<std::vector<double>> demand_map;

    /**
     * @brief Represents cost of renting specific cell, for every point on the map.
     * Shape: [grid_width, grid_height].
     * 
     */
    std::vector<std::vector<double>> land_rental_cost_map;


    /**
     * @brief It contains information about the distance from every point on the map
     * to other point of the map. It `uses calculate_distance` to achieve that.
     * Shape: [grid_width, grid_height, grid_width, grid_height].
     * 
     */
    std::vector<std::vector<std::vector<std::vector<double>>>> distances_from_to;
};