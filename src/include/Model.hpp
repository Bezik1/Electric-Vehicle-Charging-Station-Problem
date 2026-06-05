#pragma once

#include "Maps.hpp"

class Model
{
public:
    virtual void operator()(
        const Maps& input_maps,
        std::pair<double, double> _stations_powers,
        std::pair<double, double> _initial_costs,
        std::pair<double, double> _maintenance_costs
    ) = 0;

    virtual void print_solution() = 0;

    virtual void print_demand_distribution() = 0;

    virtual ~Model() = default;
};
