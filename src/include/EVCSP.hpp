#pragma once

#include <array>
#include <vector>
#include <memory>
#include <limits>

#include "gurobi_c++.h"

#include "Problem.hpp"

/**
 * @brief Electric Vehicle Charging Station Problem.
 * 
 */
class EVCSP : public Problem
{
    public:
        EVCSP(
            int grid_width,
            int grid_height,
            int max_stations_per_cell,
            double budget
        );

        void operator()(
            std::vector<std::vector<double>>&& _distances_costs_map,
            std::vector<std::vector<int>>&& _poi_map,
            std::vector<std::vector<double>>&& _demand_map,
            std::vector<std::vector<double>>&& _land_rental_cost_map,
            std::pair<double, double>&& _stations_powers,
            std::pair<double, double>&& _initial_costs,
            std::pair<double, double>&& _maintenance_costs
        );

        void build_variables() override;

        /**
         * @brief For this particular problem following criterion will be used:
         * 
         * ## Criterion
         * 
         * `min cost = operator_cost + users_cost`
         * 
         * Which, can be expressed, as:
         * 
         * ### Operator Cost
         * 
         * `operator_cost = investment_cost + land_rental_cost + maintenance_cost`
         * 
         *  ### Investment Cost
         * 
         * `investment_cost = 
         *      sum(sum(initial_costs.first * sum(l2_station_location) + initial_costs.second * sum(l3_station_location)))`
         * 
         *  ### Land Rental Cost
         * 
         * `land_rental_cost = sum(sum(land_rental_cost_map[i, j] * station_location[i, j]))`
         * 
         *  ### Maintenance Cost
         * 
         * `maintenance_cost = 
         *      sum(sum(maintenance_costs.first * sum(l2_station_location) + maintenance_costs.second * sum(l3_station_location)))`
         * 
         * ### Users Cost
         * 
         * `users_cost = users_travel_costs = 
         *      sum(sum(demand_map[m,n] * demand_allocation_map[m,n,i,j]*distances_from_to[m,n,i,j]))`
         * 
         */
        void build_criterion() override;

        /**
         * @brief Following constraints must be true, for this specific model:
         * 
         * ### Demand Satisfaction Constraint
         * 
         * All users energy demands must be fulfilled.
         * 
         * `sum(sum(demand_allocation_map[m, n, i, j])) == 1 for all (m, n).`
         * 
         * ### Station Construction Constraint
         * 
         * There is a limit of the nmbers of station, that can be build in one place.
         * 
         * `sum(l2_station_location) + sum(l3_station_location) <= max_stations_per_cell*station_location[i, j] for all (i,j).`
         * 
         * ### POI Constraint
         * 
         * Stations can only be build in the point of intrests.
         * 
         * `poi_map[i,j] == 0 => station_location[i,j] = 0 for all (i,j)`
         * 
         * ### Capacity Constraint
         * 
         * It ensures that each station, would not be overloaded.
         * 
         * * `sum(sum(demand_map[m,n] * demand_allocation_map[m,n,i,j])) <= 
         *  (l2_station_location[i,j] * stations_powers.first * 24) + 
         *  (l3_station_location[i,j] * stations_powers.second * 24) for all (i,j).`
         * 
         * ### Budget Constraint
         * 
         * Developers cannot build sth they do not have budet for.
         * 
         * `sum(sum(initial_costs.fist * l2_station_location[i,j] + initial_costs.second * l3_station_location[i,j]
         *      + land_rental_cost_map[i,j]*station_location[i,j])) <= budget`
         * 
         * ### Obstacle Constraint
         * 
         * Users cannot charge their car using station, for which they do not have path to.
         * 
         * `distances_from_to[m][n][i][j] == INF => demand_allocation_map[m][n][i][j] == 0.0`
         * 
         */
        void build_constraints() override;

        void print_solution();
        void print_demand_distribution();

    private:
        GRBEnv env;
        GRBModel model;

        int grid_width;
        int grid_height;
        int max_stations_per_cell;
        double budget;

        /**
         * @brief Map of distances from current cell (i,j) to one of it's neighbours.
         * Shape: [grid_width, grid_height].
         * 
         */
        std::vector<std::vector<double>> distances_costs_map;

        /**
         * @brief Map of Points Of Intrest. It represents map of poins on,
         * where there is a possibility of creating a station.
         * Shape: [grid_width, grid_height].
         * 
         */
        std::vector<std::vector<int>> poi_map;

        /**
         * @brief Map of user energy-demand per each point, in kWh.
         * Shape: [grid_width, grid_height].
         * 
         */
        std::vector<std::vector<double>> demand_map;

        /**
         * @brief Map of rental cost of each potential point on
         * the grid.
         * Shape: [grid_width, grid_height].
         * 
         */
        std::vector<std::vector<double>> land_rental_cost_map;

        std::pair<double, double> stations_powers;
        std::pair<double, double> initial_costs;
        std::pair<double, double> maintenance_costs;

        /**
         * @brief Value of the distance from (m,n) cell, to (i,j) cell.
         * Shape: [grid_width, grid_height, grid_width, grid_height].
         */
        std::vector<std::vector<std::vector<std::vector<double>>>> distances_from_to;


        /**
         * @brief Calculates distance from the point (m,n) to point (i,j).
         * 
         * @details Uses Dijsktra algorithm, to find shortest path. 
         * 
         * @param m x coordinate of origin point.
         * @param n y coordinate of origin point.
         * @param i x coordinate of destination point.
         * @param j y coordinate of destinatio point.
         * @return double distance.
         */
        double calculate_distance(int m, int n, int i, int j) const;

        /**
         * @brief Binary station location map. Shape: [grid_width, grid_height].
         * 
         */
        std::vector<std::vector<GRBVar>> station_location;

        /**
         * @brief Integer L2 station location map. Shape: [grid_width, grid_height].
         * 
         */
        std::vector<std::vector<GRBVar>> l2_station_location;

        /**
         * @brief Integer L3 station location map. Shape: [grid_width, grid_height].
         * 
         */
        std::vector<std::vector<GRBVar>> l3_station_location;

        /**
         * @brief  Allocation map of demands from cell (m,n) into station in (i,j).
         * Shape: [grid_width, grid_height, grid_width, grid_height].
         * 
         */
        std::vector<std::vector<std::vector<std::vector<GRBVar>>>> demand_allocation_map;

        static constexpr char* MODEL_NAME = "Electric Vehicle Charging Station Problem";
        static constexpr std::array<std::pair<int, int>, 4> DIRECTIONS = {{
            {1, 0},
            {-1, 0},
            {0, 1},
            {0, -1}
        }};
};