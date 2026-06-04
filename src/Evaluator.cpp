#include "Evaluator.hpp"
#include <cmath>
#include <algorithm>

double Evaluator::evaluate(
    const std::vector<std::vector<int>>& station_location,
    const std::vector<std::vector<int>>& l2_station_location,
    const std::vector<std::vector<int>>& l3_station_location,
    const std::vector<std::vector<std::vector<std::vector<double>>>>& demand_allocation_map
) {
    double operator_cost = 0.0;
    double users_cost = 0.0;
    double penalty = 0.0;
    double current_investment = 0.0;

    int width = context.getWidth();
    int height = context.getHeight();

    auto powers = context.getStationPowers();
    auto init_costs = context.getInitialCosts();
    auto maint_costs = context.getMaintenanceCosts();

    std::vector<std::vector<double>> station_capacities(width, std::vector<double>(height, 0.0));
    std::vector<std::vector<double>> station_load(width, std::vector<double>(height, 0.0));
    std::vector<std::vector<double>> demand_satisfied(width, std::vector<double>(height, 0.0));

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            int current_l2 = l2_station_location[i][j];
            int current_l3 = l3_station_location[i][j];
            int xi = station_location[i][j];

            if (current_l2 + current_l3 > 0) {
                station_capacities[i][j] = (current_l2 * powers.first + current_l3 * powers.second) * 24.0;
            }

            double cell_investment = (current_l2 * (init_costs.first + maint_costs.first)) +
                                     (current_l3 * (init_costs.second + maint_costs.second)) +
                                     (context.getRentalCostAt(i, j) * xi);
            current_investment += cell_investment;
            operator_cost += cell_investment;
        }
    }

    if (current_investment > context.getBudget()) {
        penalty += (current_investment - context.getBudget()) * BUDGET_WEIGHT;
    }

    for (int m = 0; m < width; ++m) {
        for (int n = 0; n < height; ++n) {
            double demand = context.getDemandAt(m, n);
            if (demand <= 0) continue;

            for (int i = 0; i < width; ++i) {
                for (int j = 0; j < height; ++j) {
                    double ratio = demand_allocation_map[m][n][i][j];
                    if (ratio > EPS) {
                        double amount = ratio * demand;
                        users_cost += amount * context.getDistance(m, n, i, j);
                        station_load[i][j] += amount;
                        demand_satisfied[m][n] += amount;
                    }
                }
            }

            if (demand_satisfied[m][n] < demand - EPS) {
                penalty += (demand - demand_satisfied[m][n]) * UNMET_DEMAND_WEIGHT;
            }
        }
    }

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            if (station_load[i][j] > station_capacities[i][j] + EPS) {
                penalty += (station_load[i][j] - station_capacities[i][j]) * OVERLOAD_WEIGHT;
            }
        }
    }

    return operator_cost + users_cost + penalty;
}