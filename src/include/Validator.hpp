#pragma once

#include "Maps.hpp"
#include "Solution.hpp"

/**
 * @brief Responsible for verifying if a Solution adheres to all hard constraints.
 * Ensures that a solution is technically and financially feasible before it can be
 * considered a valid proposal for the infrastructure layout.
 */
class Validator {
public:
    explicit Validator(
        const Maps &maps,
        int max_stations_per_cell,
        std::pair<double, double> stations_powers,
        std::pair<double, double> initial_costs,
        std::pair<double, double> maintenance_costs,
        double budget
    );

/**
     * @brief Checks if the given Solution is feasible (valid).
     * * A solution is considered valid only if it passes all spatial, structural, 
     * and financial checks. If any single criterion is violated, the method returns false.
     * 
     * Investment Calculation:
     * 
     * The total investment is the sum of costs for all grid cells, including:
     * 
     * 1. L2 Charger costs: count * (initial + maintenance)
     * 
     * 2. L3 Charger costs: count * (initial + maintenance)
     * 
     * 3. Operational costs: station_existence * rental_cost_at_location
     * 
     * Validation Criteria:
     * 
     * 1. POI Check
     * 
     * - Condition: if xi > 0 then POI must be 1
     * 
     * - Description: Stations cannot be placed on non-POI (Point of Interest) cells.
     * 
     * 2. Binary Integrity
     * 
     * - Condition: xi <= 1
     * 
     * - Description: Only one station (infrastructure unit) is allowed per cell.
     * 
     * 3. Charger Capacity
     * 
     * - Condition: (L2 + L3) <= (max_stations_per_cell * xi)
     * 
     * - Description: Total number of chargers is limited by the cell's capacity.
     * 
     * 4. Financial Limit
     * 
     * - Condition: Total_Investment <= budget + EPS
     * 
     * - Description: Total deployment and maintenance costs must not exceed the budget.
     * 
     * 5. Load Balance
     * 
     * - Condition: Station_Load <= Max_Capacity + EPS
     * 
     * - Description: The total demand assigned to a station cannot exceed its power output.
     * 
     * @param solution The complete solution object to be verified.
     * 
     * @return bool True if all criteria are met, false otherwise.
     */
    bool validate(const Solution& solution) const;

private:
    const Maps& maps;

    double budget;
    int max_stations_per_cell;

    std::pair<double, double> stations_powers;
    std::pair<double, double> initial_costs;
    std::pair<double, double> maintenance_costs;

    static constexpr double EPS = 1e-6;
};