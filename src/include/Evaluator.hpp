#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>

#include "Maps.hpp"

class Evaluator {
public:
    Evaluator(
        const Maps& maps,
        std::pair<double, double> station_powers,
        std::pair<double, double> initial_costs,
        std::pair<double, double> maintenance_costs,
        double budget
    );

    double evaluate(
        const std::vector<std::vector<int>> &station_location,
        const std::vector<std::vector<int>> &l2_station_location,
        const std::vector<std::vector<int>> &l3_station_location,
        const std::vector<std::vector<std::vector<std::vector<double>>>> &demand_allocation_map
    ) const;

private:
    const Maps& maps;
    std::pair<double, double> station_powers;
    std::pair<double, double> initial_costs;
    std::pair<double, double> maintenance_costs;
    double budget;

    static constexpr double EPS = 1e-9;
    static constexpr double OVERLOAD_WEIGHT = 10000.0;
    static constexpr double UNMET_DEMAND_WEIGHT = 20000.0;
    static constexpr double BUDGET_WEIGHT = 500.0;
};