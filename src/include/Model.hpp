#pragma once

class Model
{
public:
    virtual void operator()(
        std::vector<std::vector<double>> _distances_costs_map,
        std::vector<std::vector<int>> _poi_map,
        std::vector<std::vector<double>> _demand_map,
        std::vector<std::vector<double>> _land_rental_cost_map,
        std::pair<double, double> _stations_powers,
        std::pair<double, double> _initial_costs,
        std::pair<double, double> _maintenance_costs
    ) = 0;

    virtual void print_solution() = 0;

    virtual void print_demand_distribution() = 0;

    virtual ~Model() = default;
};
