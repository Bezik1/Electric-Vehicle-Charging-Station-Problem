#include "LinearProgrammingModel.hpp"
#include "MetaheuristicModel.hpp"
#include "Evaluator.hpp"
#include "Maps.hpp"

#include <nlohmann/json.hpp>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <print>

struct Context {
    int width;
    int height;
    int max_stations_per_cell;
    double budget;
    std::pair<double, double> stations_powers;
    std::pair<double, double> initial_costs;
    std::pair<double, double> maintenance_costs;
    const Maps* map_data;
};

double objective_function(const dlib::matrix<double, 0, 1>& params, const Context& ctx) {
    int bacteria_count = static_cast<int>(params(0));
    int hemotaxis_steps = static_cast<int>(params(1));
    int swimming_steps = static_cast<int>(params(2));
    int reproduction_steps = static_cast<int>(params(3));
    int elimination_steps = static_cast<int>(params(4));
    double elimination_prob = params(5);

    MetaheuristicModel model(
        ctx.max_stations_per_cell,
        bacteria_count, hemotaxis_steps, swimming_steps,
        reproduction_steps, elimination_steps, elimination_prob, ctx.budget
    );

    model(
        *ctx.map_data,
        ctx.stations_powers,
        ctx.initial_costs, 
        ctx.maintenance_costs
    );

    return model.best_cost;
}

int main() {
    std::ifstream map_file("data/maps/maps.json");
    std::ifstream hyperparameters_file("data/parameters/hyperparameters.json");
    std::ifstream problem_parameters_file("data/parameters/problem_parameters.json");
    
    nlohmann::json maps = nlohmann::json::parse(map_file);
    nlohmann::json hyperparameters = nlohmann::json::parse(hyperparameters_file);
    nlohmann::json problem_parameters = nlohmann::json::parse(problem_parameters_file);

    Context ctx;
    ctx.width = maps["width"];
    ctx.height = maps["height"];
    
    auto distances_costs_map = maps["distances_map"].get<std::vector<std::vector<double>>>();
    auto demand_map = maps["demand_map"].get<std::vector<std::vector<double>>>();
    auto poi_map = maps["poi_map"].get<std::vector<std::vector<int>>>();
    auto land_rental_cost_map = maps["land_rental_cost_map"].get<std::vector<std::vector<double>>>();

    ctx.max_stations_per_cell = problem_parameters["max_stations_per_cell"].get<int>();
    ctx.budget = problem_parameters["budget"].get<double>();
    ctx.stations_powers = problem_parameters["stations_powers"].get<std::pair<double, double>>();
    ctx.initial_costs = problem_parameters["initial_costs"].get<std::pair<double, double>>();
    ctx.maintenance_costs = problem_parameters["maintenance_costs"].get<std::pair<double, double>>();

    double mip_gap = problem_parameters["mip_gap"].get<double>();

    Maps map_data(
        ctx.width, 
        ctx.height, 
        distances_costs_map, 
        poi_map, 
        demand_map, 
        land_rental_cost_map
    );
    ctx.map_data = &map_data;

    dlib::matrix<double, 0, 1> lower_bound(6), upper_bound(6);
    lower_bound = 10, 5, 1, 2, 1, 0.01;
    upper_bound = 100, 50, 10, 10, 5, 0.50;

    std::println("Starting hyperparameter optimization...");

    auto result = dlib::find_min_global(
        [&ctx](const dlib::matrix<double, 0, 1>& p) { return objective_function(p, ctx); },
        lower_bound,
        upper_bound,
        dlib::max_function_calls(50)
    );

    MetaheuristicModel meta_solver(
        ctx.max_stations_per_cell,
        static_cast<int>(result.x(0)), static_cast<int>(result.x(1)),
        static_cast<int>(result.x(2)), static_cast<int>(result.x(3)),
        static_cast<int>(result.x(4)), result.x(5), ctx.budget
    );

    meta_solver(
        map_data,
        ctx.stations_powers,
        ctx.initial_costs,
        ctx.maintenance_costs
    );

    LinearProgrammingModel solver(ctx.width, ctx.height, ctx.max_stations_per_cell, ctx.budget, mip_gap);
    solver(
        map_data,
        ctx.stations_powers,
        ctx.initial_costs,
        ctx.maintenance_costs
    );

    Evaluator evaluator(map_data, ctx.stations_powers, ctx.initial_costs, ctx.maintenance_costs, ctx.budget);

    const auto& solver_solution = solver.get_solution();
    const auto& meta_solver_solution = *meta_solver.best_bacterium;

    std::println("--- Optimization Finished ---");
    std::println("Best Hyperparameters found: B:{} H:{} S:{} R:{} E:{} Prob:{:.2f}", 
        static_cast<int>(result.x(0)), static_cast<int>(result.x(1)), 
        static_cast<int>(result.x(2)), static_cast<int>(result.x(3)), 
        static_cast<int>(result.x(4)), result.x(5));

    double lp_cost = evaluator.evaluate(
        solver_solution.station_location, 
        solver_solution.l2_station_location, 
        solver_solution.l3_station_location, 
        solver_solution.demand_allocation_map
    );

    double meta_cost = evaluator.evaluate(
        meta_solver_solution.station_location, 
        meta_solver_solution.l2_station_location, 
        meta_solver_solution.l3_station_location, 
        meta_solver_solution.demand_allocation_map
    );

    std::println("Linear Programming Model: {:.2f}", lp_cost);
    std::println("Meta Heuristic Model (Tuned): {:.2f}", meta_cost);

    return 0;
}