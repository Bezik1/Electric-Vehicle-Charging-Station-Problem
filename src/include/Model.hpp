#pragma once

#include "Maps.hpp"
#include "Solution.hpp"

/**
 * @brief Represents interface for implementing EVCSP solving model.
 * 
 */
class Model
{
public:
    /**
     * @brief This operator should be used, as function to call
     * to run the model.
     * 
     * @param input_maps Maps describing the problem infrastucture.
     * @param _stations_powers Powers of each station in kW.
     * @param _initial_costs Cost of building each station.
     * @param _maintenance_costs Cost of maintaining each station.
     */
    virtual void operator()(
        const Maps& input_maps,
        std::pair<double, double> _stations_powers,
        std::pair<double, double> _initial_costs,
        std::pair<double, double> _maintenance_costs
    ) = 0;

    /**
     * @brief Gets the solution of the model.
     * 
     * @return Solution Locations of potential stations, l2 and l3 maps,
     * with user demand allocation map.
     */
    virtual Solution get_solution() const = 0;

    virtual ~Model() = default;
};
