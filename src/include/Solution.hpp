#pragma once

#include <vector>

/**
 * @brief Represents individual solution of the EVCSP.
 * 
 */
struct Solution {
    /**
     * @brief Locations of the potential stations. It's the
     * binary vector of shape [grid_width, grid_height].
     * 
     */
    std::vector<std::vector<int>> station_location;

    /**
     * @brief Locations of the l2 stations. It's the
     * integer vector of shape [grid_width, grid_height].
     */
    std::vector<std::vector<int>> l2_station_location;

    /**
     * @brief Locations of the l3 stations. It's the
     * integer vector of shape [grid_width, grid_height].
     */
    std::vector<std::vector<int>> l3_station_location;

    /**
     * @brief Represents user demand allocations map. It maps demands from cell
     * (m, n) to the station located at (i, j). It is of shape:
     * [grid_width, grid_height, grid_width, grid_height].
     * 
     */
    std::vector<std::vector<std::vector<std::vector<double>>>> demand_allocation_map;

    /**
     * @brief Represents cost of the individual solution.
     */
    double total_cost;
};