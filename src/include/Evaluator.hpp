#pragma once

#include <vector>
#include <algorithm>
#include <cmath>

#include "EvaluationContext.hpp"

class Evaluator {
public:
    explicit Evaluator(const EvaluationContext& ctx) : context(ctx) {}

    double evaluate(
        const std::vector<std::vector<int>> &station_location,
        const std::vector<std::vector<int>> &l2_station_location,
        const std::vector<std::vector<int>> &l3_station_location,
        const std::vector<std::vector<std::vector<std::vector<double>>>> &demand_allocation_map
    );

private:
    const EvaluationContext& context;

    static constexpr double EPS = 1e-9;
    static constexpr double OVERLOAD_WEIGHT = 10000.0;
    static constexpr double UNMET_DEMAND_WEIGHT = 20000.0;
    static constexpr double BUDGET_WEIGHT = 500.0;
};