#pragma once

#include <array>
#include <vector>
#include <memory>
#include <limits>
#include <variant>

#include "Maps.hpp"
#include "Model.hpp"
#include "Evaluator.hpp"
#include "Bacterium.hpp"

class MetaheuristicModel : public Model
{
public:
    MetaheuristicModel(
        int max_stations_per_cell,
        int bacteria_count,
        int hemotaxis_steps,
        int swimming_steps,
        int reproduction_steps,
        int elimination_dispersal_steps,
        double elimination_dispersal_prob,
        double budget
    );

    void operator()(
        const Maps& input_maps,
        std::pair<double, double> _stations_powers,
        std::pair<double, double> _initial_costs,
        std::pair<double, double> _maintenance_costs
    ) override;

    void optimization_loop();
    
    void print_solution() override;
    
    void print_demand_distribution() override;
    
    std::unique_ptr<Bacterium> best_bacterium = nullptr;
    double best_cost = std::numeric_limits<double>::infinity();

private:
    double calculate_swarming_effect(const Bacterium &current_bacterium) const;

    bool is_valid(const std::unique_ptr<Bacterium>& bacterium) const;

    const Maps* maps = nullptr;
    
    std::vector<std::unique_ptr<Bacterium>> colony;
    
    std::unique_ptr<Evaluator> evaluator = nullptr;
    
    std::pair<double, double> stations_powers;
    std::pair<double, double> initial_costs;
    std::pair<double, double> maintenance_costs;

    int max_stations_per_cell;
    double budget;
    int bacteria_count; 
    int hemotaxis_steps;
    int swimming_steps;
    int reproduction_steps;
    int elimination_dispersal_steps;
    double elimination_dispersal_prob;
};