#pragma once

#include <array>
#include <vector>
#include <memory>
#include <limits>

#include "gurobi_c++.h"

#include "Model.hpp"
#include "GurobiModel.hpp"
#include "Solution.hpp"
#include "SolutionPrinter.hpp"
#include "Maps.hpp"

/**
 * @brief Electric Vehicle Charging Station Problem.
 * */
class LinearProgrammingModel : public GurobiModel, public Model
{
    public:
        explicit LinearProgrammingModel(
            int grid_width,
            int grid_height,
            int max_stations_per_cell,
            double budget
        );

        explicit LinearProgrammingModel(
            int grid_width,
            int grid_height,
            int max_stations_per_cell,
            double budget,
            double mip_gap
        );

        void operator()(
            const Maps &_map_data,
            std::pair<double, double> _stations_powers,
            std::pair<double, double> _initial_costs,
            std::pair<double, double> _maintenance_costs
        ) override;

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

        Solution get_solution() const override;

        std::unique_ptr<SolutionPrinter> printer;

    private:
        GRBEnv env;
        GRBModel model;

        int grid_width;
        int grid_height;
        int max_stations_per_cell;
        double budget;

        const Maps* map_data;

        std::pair<double, double> stations_powers;
        std::pair<double, double> initial_costs;
        std::pair<double, double> maintenance_costs;

        /**
         * @brief Binary station location map. Shape: [grid_width, grid_height].
         * */
        std::vector<std::vector<GRBVar>> station_location;

        /**
         * @brief Integer L2 station location map. Shape: [grid_width, grid_height].
         * */
        std::vector<std::vector<GRBVar>> l2_station_location;

        /**
         * @brief Integer L3 station location map. Shape: [grid_width, grid_height].
         * */
        std::vector<std::vector<GRBVar>> l3_station_location;

        /**
         * @brief  Allocation map of demands from cell (m,n) into station in (i,j).
         * Shape: [grid_width, grid_height, grid_width, grid_height].
         * */
        std::vector<std::vector<std::vector<std::vector<GRBVar>>>> demand_allocation_map;

        static constexpr double EPSILON = 1e-4; 
        static constexpr double DEFAULT_MIP_GAP = 0.01;
        static constexpr char *MODEL_NAME = "Electric Vehicle Charging Station Problem";
};