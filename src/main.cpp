#include "EVCSP.hpp"
#include <nlohmann/json.hpp>
#include <fstream>

int main() 
{
    std::ifstream mapFile("data/maps.json");
    std::ifstream hyperparametersFile("data/hyperparameters.json");

    nlohmann::json maps = nlohmann::json::parse(mapFile);
    nlohmann::json hyperparameters = nlohmann::json::parse(hyperparametersFile);

    int width = maps["width"];
    int height = maps["height"];

    auto distances_costs_map = maps["distances_map"].get<std::vector<std::vector<double>>>();
    auto demand_map = maps["demand_map"].get<std::vector<std::vector<double>>>();
    auto poi_map = maps["poi_map"].get<std::vector<std::vector<int>>>();
    auto land_rental_cost_map = maps["land_rental_cost_map"].get<std::vector<std::vector<double>>>();

    int max_stations_per_cell = hyperparameters["max_stations_per_cell"].get<int>();
    double budget = hyperparameters["budget"].get<double>();

    auto stations_powers = hyperparameters["stations_powers"].get<std::pair<double, double>>();
    auto initial_costs = hyperparameters["initial_costs"].get<std::pair<double, double>>();
    auto maintenance_costs = hyperparameters["maintenance_costs"].get<std::pair<double, double>>();

    double mip_gap = hyperparameters["mip_gap"].get<double>();
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