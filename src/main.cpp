#include "EVCSP.hpp"
#include <nlohmann/json.hpp>
#include <fstream>

int main() 
{
    std::ifstream file("data/maps.json");
    nlohmann::json data = nlohmann::json::parse(file);

    int width = data["width"];
    int height = data["height"];

    auto distances_costs_map = data["distances_map"].get<std::vector<std::vector<double>>>();
    auto demand_map = data["demand_map"].get<std::vector<std::vector<double>>>();
    auto poi_map = data["poi_map"].get<std::vector<std::vector<int>>>();
    auto land_rental_cost_map = data["land_rental_cost_map"].get<std::vector<std::vector<double>>>();

    int max_stations_per_cell = 1;
    double budget = 400'000.0;

    std::pair<double, double> stations_powers = {25.0, 250.0}; 
    std::pair<double, double> initial_costs = {50'000.0, 150'000.0}; 
    std::pair<double, double> maintenance_costs = {2500.0, 9000.0};

    double mip_gap = 0.007;
    EVCSP solver(width, height, max_stations_per_cell, budget, mip_gap);

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