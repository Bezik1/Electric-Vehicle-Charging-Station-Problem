#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>

#include "Maps.hpp"

/**
 * @brief Evaluator class responsible for calculating the fitness of a Bacterium solution.
 *
 */
class Evaluator
{
public:
    explicit Evaluator(
        const Maps& maps,
        std::pair<double, double> station_powers,
        std::pair<double, double> initial_costs,
        std::pair<double, double> maintenance_costs,
        int max_stations_per_cell,
        double budget
    );

    /**
     * @brief Calculates the total cost of the solution, which is the sum of
     * operator expenses, user costs, and constraint violation penalties.
     *
     * The total cost C_total is calculated as:
     * 
     * C_total = C_operator + C_users + sum(P)
     *
     * Where:
     *
     * - Operator Cost C_operator: Sum of the initial purchase costs (Initial),
     * ongoing maintenance (Maintenance), and cell rental fees (Rental).
     *
     * - User Cost C_users: Sum of products of allocated demand,
     * travel distance to the station, and the USER_DISTANCE_MULTIPLIER.
     *
     * * * PENALTY CASES (P):
     *
     * 1. Charger Limit Per Cell Violation:
     *
     * - Condition: Number of chargers > xi * max_stations_per_cell
     * 
     * - Amount: HARD_CONSTRAINT_OFFSET + (excess * CAPACITY_VIOLATION_WEIGHT)
     *
     * 2. Station Location Not Binary Structure:
     *
     * - Condition: xi > 1 (more than one station allowed per grid cell)
     * 
     * - Amount: HARD_CONSTRAINT_OFFSET * xi
     *
     * 3. Invalid Station Location:
     *
     * - Condition: Station placed in a cell where POI == 0
     * 
     * - Amount: HARD_CONSTRAINT_OFFSET
     *
     * 4. Budget Exceeded:
     *
     * - Condition: Total investment > budget
     * 
     * - Amount: HARD_CONSTRAINT_OFFSET + (excess * BUDGET_WEIGHT)
     *
     * 5. Unmet Demand:
     *
     * - Condition: Sum of allocation for a demand point < total demand
     * 
     * - Amount: (unmet_demand * UNMET_DEMAND_WEIGHT)
     *
     * 6. Station Overload:
     *
     * - Condition: Assigned load (kWh) > station capacity (kWh)
     * 
     * - Amount: (overload_amount * OVERLOAD_WEIGHT)
     * 
     * @param station_location station_location Binary representation of station locations.
     * 
     * @param station_location l2_station_location L2 stations count at given cell.
     * 
     * @param station_location l3_station_location L3 stations count at given cell.
     * 
     * @param station_location demand_allocation_map Electric vehicle user demands allocation map.
     *
     * @return double The objective function value (the lower the value, the better the solution).
     */
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
    int max_stations_per_cell;
    double budget;

    static constexpr double EPS = 1e-6;
    static constexpr double HARD_CONSTRAINT_OFFSET = 1e7;
    static constexpr double OVERLOAD_WEIGHT = 5000.0;
    static constexpr double UNMET_DEMAND_WEIGHT = 10000.0;
    static constexpr double CAPACITY_VIOLATION_WEIGHT = 500000.0;
    static constexpr double BUDGET_WEIGHT = 10.0;
    static constexpr double USER_DISTANCE_MULTIPLIER = 2.0;
};