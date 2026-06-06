#include "Evaluator.hpp"

Evaluator::Evaluator(
    const Maps& maps,
    std::pair<double, double> station_powers,
    std::pair<double, double> initial_costs,
    std::pair<double, double> maintenance_costs,
    int max_stations_per_cell,
    double budget
) : maps{maps},
    station_powers{std::move(station_powers)},
    initial_costs{std::move(initial_costs)},
    maintenance_costs{std::move(maintenance_costs)},
    max_stations_per_cell{max_stations_per_cell},
    budget{budget} 
{}

double Evaluator::evaluate(
    const std::vector<std::vector<int>>& station_location,
    const std::vector<std::vector<int>>& l2_station_location,
    const std::vector<std::vector<int>>& l3_station_location,
    const std::vector<std::vector<std::vector<std::vector<double>>>>& demand_allocation_map
) const {
    double operator_cost = 0.0;
    double users_cost = 0.0;
    double penalty = 0.0;
    double current_investment = 0.0;

    int width = maps.get_width();
    int height = maps.get_height();

    std::vector<std::vector<double>> station_capacities(width, std::vector<double>(height, 0.0));
    std::vector<std::vector<double>> station_load(width, std::vector<double>(height, 0.0));
    std::vector<std::vector<double>> demand_satisfied(width, std::vector<double>(height, 0.0));

    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            int current_l2 = l2_station_location[i][j];
            int current_l3 = l3_station_location[i][j];
            int xi = station_location[i][j];

            int total_chargers = current_l2 + current_l3;
            
            if (total_chargers > xi * max_stations_per_cell)
                penalty += 
                    HARD_CONSTRAINT_OFFSET + (total_chargers - (xi * max_stations_per_cell)) * CAPACITY_VIOLATION_WEIGHT;

            if(xi > 1)
                penalty += HARD_CONSTRAINT_OFFSET * xi;

            if (xi > 0 && maps.get_poi_at(i, j) == 0)
                penalty += HARD_CONSTRAINT_OFFSET;

            if (total_chargers > 0)
                station_capacities[i][j] =
                    (current_l2 * station_powers.first + current_l3 * station_powers.second) * 24.0;

            double cell_investment = (current_l2 * (initial_costs.first + maintenance_costs.first)) +
                                     (current_l3 * (initial_costs.second + maintenance_costs.second)) +
                                     (maps.get_rental_cost_at(i, j) * xi);
            current_investment += cell_investment;
            operator_cost += cell_investment;
        }
    }

    if (current_investment > budget)
        penalty += HARD_CONSTRAINT_OFFSET + (current_investment - budget) * BUDGET_WEIGHT;

    for (int m = 0; m < width; ++m)
    {
        for (int n = 0; n < height; ++n)
        {
            double demand = maps.get_demand_at(m, n);
            if (demand <= 0) continue;

            for (int i = 0; i < width; ++i)
            {
                for (int j = 0; j < height; ++j)
                {
                    double ratio = demand_allocation_map[m][n][i][j];
                    if (ratio > EPS) {
                        double amount = ratio * demand;
                        users_cost += amount * maps.get_distance(m, n, i, j) * USER_DISTANCE_MULTIPLIER; 
                        station_load[i][j] += amount;
                        demand_satisfied[m][n] += amount;
                    }
                }
            }

            if (demand_satisfied[m][n] < demand - EPS)
                penalty += (demand - demand_satisfied[m][n]) * UNMET_DEMAND_WEIGHT; 
        }
    }

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
            if (station_load[i][j] > station_capacities[i][j] + EPS)
                penalty += (station_load[i][j] - station_capacities[i][j]) * OVERLOAD_WEIGHT;

    return operator_cost + users_cost + penalty;
}