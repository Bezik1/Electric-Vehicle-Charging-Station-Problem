#pragma once

#include <vector>
#include <memory>

#include "gurobi_c++.h"

#include "Problem.hpp"

class EVCSP : public Problem
{
    public:
        EVCSP(
            double _mean_station_energy_demand,
            int _num_stations,
            int _map_width,
            int _map_height);
        
        void operator()(
            std::vector<std::vector<int>>&& _connectivity_map,
            std::vector<std::vector<double>>&& _activity_map,
            std::vector<std::vector<double>>&& _attractiveness_map,
            std::vector<std::vector<double>>&& _energy_capacity_map
        );

        void build_variables() override;
        void build_criterion() override;
        void build_constraints() override;
        void print_map();

    private:
        inline static const std::string MODEL_NAME =
            "Electric-Vehicle-Charging-Station-Problem";

        static constexpr double COMMERCIAL_POTENTIAL_COEFF = 1.0;
        static constexpr double NETWORK_SECURITY_COEFF = 1.5;
        static constexpr double MINIMAL_ENERGY_DEMAND = 0.1;

        GRBEnv env;
        GRBModel model;

        int num_stations;
        int map_width;
        int map_height;

        double mean_station_energy_demand;

        std::vector<std::vector<int>> connectivity_map;
        std::vector<std::vector<double>> activity_map;
        std::vector<std::vector<double>> attractiveness_map;
        std::vector<std::vector<double>> energy_capacity_map;

        std::vector<std::vector<GRBVar>> grid_vars;
};