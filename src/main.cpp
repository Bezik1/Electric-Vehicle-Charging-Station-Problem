#include "EVCSP.hpp"

int main() 
{
    int grid_width = 5;
    int grid_height = 5;
    int max_stations_per_cell = 3;
    double budget = 60000.0;

    std::vector<std::vector<int>> poi_map = {
        {1, 1, 0, 1, 1},
        {1, 0, 1, 1, 0},
        {1, 1, 1, 1, 1},
        {0, 1, 0, 1, 1},
        {1, 1, 1, 0, 1}
    };

    std::vector<std::vector<double>> demand_map = {
        {120.0, 80.0, 10.0, 50.0, 200.0},
        {90.0, 0.0, 150.0, 70.0, 30.0},
        {40.0, 110.0, 300.0, 85.0, 60.0},
        {5.0, 95.0, 0.0, 130.0, 45.0},
        {0.0, 60.0, 75.0, 20.0, 115.0}
    };

    std::vector<std::vector<double>> land_rental_cost_map = {
        {1200.0, 1100.0,  800.0, 1500.0, 2200.0},
        {1000.0,  700.0, 1400.0, 1300.0,  900.0},
        { 950.0, 1350.0, 2500.0, 1600.0, 1050.0},
        { 600.0, 1150.0,  500.0, 1200.0,  850.0},
        { 800.0,  900.0, 1100.0,  700.0, 1300.0}
    };

    std::vector<std::vector<double>> distances_costs_map(
        grid_width, std::vector<double>(grid_height, 1.0)
    );

    for (int i = 0; i < grid_width; ++i)
    {
        for (int j = 0; j < grid_height; ++j)
        {
            if(poi_map[i][j] == 0)
                distances_costs_map[i][j] = std::numeric_limits<double>::infinity();
        }
    }

    std::pair<double, double> stations_powers = {22.0, 150.0};
    std::pair<double, double> initial_costs = {3500.0, 15000.0};
    std::pair<double, double> maintenance_costs = {200.0, 800.0};

    EVCSP solver(grid_width, grid_height, max_stations_per_cell, budget);

    solver(
        std::move(distances_costs_map),
        std::move(poi_map),
        std::move(demand_map),
        std::move(land_rental_cost_map),
        std::move(stations_powers),
        std::move(initial_costs),
        std::move(maintenance_costs)
    );

    solver.print_solution();
    solver.print_demand_distribution();

    return 0;
}