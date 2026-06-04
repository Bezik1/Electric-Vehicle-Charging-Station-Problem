#include "LinearProgrammingModel.hpp"
#include "MetaheuristicModel.hpp"
#include "Evaluator.hpp"

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
    std::vector<std::vector<double>> distances_costs_map;
    std::vector<std::vector<int>> poi_map;
    std::vector<std::vector<double>> demand_map;
    std::vector<std::vector<double>> land_rental_cost_map;
    std::pair<double, double> stations_powers;
    std::pair<double, double> initial_costs;
    std::pair<double, double> maintenance_costs;
};

double objective_function(const dlib::matrix<double, 0, 1>& params, const Context& ctx) {
    int bacteria_count = static_cast<int>(params(0));
    int hemotaxis_steps = static_cast<int>(params(1));
    int swimming_steps = static_cast<int>(params(2));
    int reproduction_steps = static_cast<int>(params(3));
    int elimination_steps = static_cast<int>(params(4));
    double elimination_prob = params(5);

    MetaheuristicModel model(
        ctx.width, ctx.height, ctx.max_stations_per_cell,
        bacteria_count, hemotaxis_steps, swimming_steps,
        reproduction_steps, elimination_steps, elimination_prob, ctx.budget
    );

    model(
        ctx.distances_costs_map, ctx.poi_map, ctx.demand_map,
        ctx.land_rental_cost_map, ctx.stations_powers,
        ctx.initial_costs, ctx.maintenance_costs
    );

    return model.best_cost;
}

int main() {
    std::ifstream mapFile("data/maps.json");
    std::ifstream hyperparametersFile("data/hyperparameters.json");
    
    nlohmann::json maps = nlohmann::json::parse(mapFile);
    nlohmann::json hyperparameters = nlohmann::json::parse(hyperparametersFile);

    Context ctx;
    ctx.width = maps["width"];
    ctx.height = maps["height"];
    ctx.distances_costs_map = maps["distances_map"].get<std::vector<std::vector<double>>>();
    ctx.demand_map = maps["demand_map"].get<std::vector<std::vector<double>>>();
    ctx.poi_map = maps["poi_map"].get<std::vector<std::vector<int>>>();
    ctx.land_rental_cost_map = maps["land_rental_cost_map"].get<std::vector<std::vector<double>>>();
    ctx.max_stations_per_cell = hyperparameters["max_stations_per_cell"].get<int>();
    ctx.budget = hyperparameters["budget"].get<double>();
    ctx.stations_powers = hyperparameters["stations_powers"].get<std::pair<double, double>>();
    ctx.initial_costs = hyperparameters["initial_costs"].get<std::pair<double, double>>();
    ctx.maintenance_costs = hyperparameters["maintenance_costs"].get<std::pair<double, double>>();

    dlib::matrix<double, 0, 1> lower_bound(6), upper_bound(6);
    lower_bound = 50, 20, 2, 2, 2, 0.05;
    upper_bound = 250, 120, 10, 8, 8, 0.40;

    auto result = dlib::find_min_global(
        [&ctx](const dlib::matrix<double, 0, 1>& p) { return objective_function(p, ctx); },
        lower_bound,
        upper_bound,
        dlib::max_function_calls(60)
    );

    MetaheuristicModel meta_solver(
        ctx.width, ctx.height, ctx.max_stations_per_cell,
        static_cast<int>(result.x(0)), static_cast<int>(result.x(1)),
        static_cast<int>(result.x(2)), static_cast<int>(result.x(3)),
        static_cast<int>(result.x(4)), result.x(5), ctx.budget
    );

    meta_solver(
        ctx.distances_costs_map, ctx.poi_map, ctx.demand_map,
        ctx.land_rental_cost_map, ctx.stations_powers,
        ctx.initial_costs, ctx.maintenance_costs
    );

    LinearProgrammingModel lp_solver(ctx.width, ctx.height, ctx.max_stations_per_cell, ctx.budget, hyperparameters["mip_gap"].get<double>());
    lp_solver(ctx.distances_costs_map, ctx.poi_map, ctx.demand_map, ctx.land_rental_cost_map, ctx.stations_powers, ctx.initial_costs, ctx.maintenance_costs);

    Evaluator evaluator(meta_solver);
    const auto& lp_sol = lp_solver.get_solution();
    const auto& meta_sol = *meta_solver.best_bacterium;

    std::println("Tuned Parameters: B:{} H:{} S:{} R:{} E:{} Prob:{:.2f}", 
        static_cast<int>(result.x(0)), static_cast<int>(result.x(1)), 
        static_cast<int>(result.x(2)), static_cast<int>(result.x(3)), 
        static_cast<int>(result.x(4)), result.x(5));

    std::println("Linear Programming Model: {:.2f}", evaluator.evaluate(lp_sol.station_location, lp_sol.l2_station_location, lp_sol.l3_station_location, lp_sol.demand_allocation_map));
    std::println("Meta Heuristic Model: {:.2f}", evaluator.evaluate(meta_sol.station_location, meta_sol.l2_station_location, meta_sol.l3_station_location, meta_sol.demand_allocation_map));

    return 0;
}