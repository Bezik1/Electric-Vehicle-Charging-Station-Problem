#include "LinearProgrammingModel.hpp"
#include "MetaheuristicModel.hpp"
#include "Evaluator.hpp"
#include "Validator.hpp"
#include "Maps.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

int main() 
{
    std::ifstream map_file("data/maps/30x30_map.json");
    std::ifstream hyperparameters_file("data/parameters/hyperparameters.json");
    std::ifstream problem_parameters_file("data/parameters/problem_parameters.json");
    
    nlohmann::json maps = nlohmann::json::parse(map_file);
    nlohmann::json hyperparameters = nlohmann::json::parse(hyperparameters_file);
    nlohmann::json problem_parameters = nlohmann::json::parse(problem_parameters_file);

    int width = maps["width"];
    int height = maps["height"];

    auto distances_costs_map = maps["distances_map"].get<std::vector<std::vector<double>>>();
    auto demand_map = maps["demand_map"].get<std::vector<std::vector<double>>>();
    auto poi_map = maps["poi_map"].get<std::vector<std::vector<int>>>();
    auto land_rental_cost_map = maps["land_rental_cost_map"].get<std::vector<std::vector<double>>>();

    int max_stations_per_cell = problem_parameters["max_stations_per_cell"].get<int>();
    double budget = problem_parameters["budget"].get<double>();

    auto stations_powers = problem_parameters["stations_powers"].get<std::pair<double, double>>();
    auto initial_costs = problem_parameters["initial_costs"].get<std::pair<double, double>>();
    auto maintenance_costs = problem_parameters["maintenance_costs"].get<std::pair<double, double>>();

    double mip_gap = problem_parameters["mip_gap"].get<double>();

    int bacteria_count = hyperparameters["bacteria_count"].get<int>();
    int hemotaxis_steps = hyperparameters["hemotaxis_steps"].get<int>();
    int swimming_steps = hyperparameters["swimming_steps"].get<int>();
    int reproduction_steps = hyperparameters["reproduction_steps"].get<int>();

    int elimination_dispersal_steps = hyperparameters["elimination_dispersal_steps"].get<int>();
    double elimination_dispersal_prob = hyperparameters["elimination_dispersal_prob"].get<double>();

    double d_attract = hyperparameters["d_attract"].get<double>();
    double w_attract = hyperparameters["w_attract"].get<double>();
    double h_repellant = hyperparameters["h_repellant"].get<double>();
    double w_repellant = hyperparameters["w_repellant"].get<double>();

    Maps map_data{
        width, 
        height, 
        distances_costs_map, 
        poi_map, 
        demand_map, 
        land_rental_cost_map
    };

    LinearProgrammingModel solver{width, height, max_stations_per_cell, budget, mip_gap};

    MetaheuristicModel meta_solver{
        max_stations_per_cell, bacteria_count,
        hemotaxis_steps, swimming_steps, reproduction_steps,
        elimination_dispersal_steps, elimination_dispersal_prob, 
        d_attract, w_attract, h_repellant, w_repellant, budget
    };

    Evaluator evaluator{
        map_data,
        stations_powers,
        initial_costs,
        maintenance_costs,
        max_stations_per_cell,
        budget
    };

    Validator validator{
        map_data,
        max_stations_per_cell,
        stations_powers,
        initial_costs,
        maintenance_costs,
        budget
    };

    meta_solver(
        map_data,
        stations_powers,
        initial_costs,
        maintenance_costs
    );

    const auto& meta_solver_solution = meta_solver.get_solution();
    meta_solver.printer->print_map(meta_solver_solution);
    meta_solver.printer->print_demand_distribution(meta_solver_solution);

    solver(
        map_data,
        stations_powers,
        initial_costs,
        maintenance_costs
    );

    const auto& solver_solution = solver.get_solution();
    solver.printer->print_map(solver_solution);
    solver.printer->print_demand_distribution(solver_solution);

    double linear_programming_cost = evaluator.evaluate(
        solver_solution.station_location,
        solver_solution.l2_station_location,
        solver_solution.l3_station_location,
        solver_solution.demand_allocation_map
    );

    std::println("Linear Programming Model: Cost: {:.2f}, Valid: {}.",
        linear_programming_cost,
        validator.validate(solver_solution)
    );

    double meta_heuristic_cost = evaluator.evaluate(
        meta_solver_solution.station_location,
        meta_solver_solution.l2_station_location,
        meta_solver_solution.l3_station_location,
        meta_solver_solution.demand_allocation_map
    );

    std::println(
        "Meta Heuristic Model: Cost: {:.2f}, Valid: {}.",
        meta_heuristic_cost,
        validator.validate(meta_solver_solution)
    );

    double optimality_gap = (meta_heuristic_cost - linear_programming_cost) / meta_heuristic_cost;

    std::println(
        "Optimality Gap: Cost: {:.2f}%",
        optimality_gap * 100
    );

    return 0;
}