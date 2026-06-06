#include "MetaheuristicModel.hpp"
#include <queue>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

MetaheuristicModel::MetaheuristicModel(
        int max_stations_per_cell,
        int bacteria_count,
        int hemotaxis_steps,
        int swimming_steps,
        int reproduction_steps,
        int elimination_dispersal_steps,
        double elimination_dispersal_prob,
        double d_attract,
        double w_attract,
        double h_repellant,
        double w_repellant,
        double budget
) : max_stations_per_cell{max_stations_per_cell},
    bacteria_count{bacteria_count},
    hemotaxis_steps{hemotaxis_steps},
    swimming_steps{swimming_steps},
    reproduction_steps{reproduction_steps},
    elimination_dispersal_steps{elimination_dispersal_steps},
    elimination_dispersal_prob{elimination_dispersal_prob},
    d_attract{d_attract},
    w_attract{w_attract},
    h_repellant{h_repellant},
    w_repellant{w_repellant},
    budget{budget}
{}

void MetaheuristicModel::operator()(
            const Maps& input_maps,
            std::pair<double, double> _stations_powers,
            std::pair<double, double> _initial_costs,
            std::pair<double, double> _maintenance_costs
        ) 
{
    maps = &input_maps;
    stations_powers = std::move(_stations_powers);
    initial_costs = std::move(_initial_costs);
    maintenance_costs = std::move(_maintenance_costs);

    evaluator = std::make_unique<Evaluator>(
        *maps,
        stations_powers,
        initial_costs,
        maintenance_costs,
        max_stations_per_cell,
        budget
    );
    validator = std::make_unique<Validator>(
        *maps,
        max_stations_per_cell,
        stations_powers,
        initial_costs,
        maintenance_costs,
        budget
    );

    printer = std::make_unique<SolutionPrinter>(*maps);

    optimization_loop();
}

double MetaheuristicModel::calculate_swarming_effect(const Bacterium& current_bacterium) const
{
    double swarming_cost = 0.0;

    for (const auto& other : colony)
    {
        if (!other) continue;

        double distance_sq = 0.0;
        for (int i = 0; i < maps->get_width(); ++i)
        {
            for (int j = 0; j < maps->get_height(); ++j)
            {
                double diff_l2 = current_bacterium.l2_station_location[i][j] - other->l2_station_location[i][j];
                double diff_l3 = current_bacterium.l3_station_location[i][j] - other->l3_station_location[i][j];
                distance_sq += (diff_l2 * diff_l2) + (diff_l3 * diff_l3);
            }
        }
        swarming_cost += -d_attract * std::exp(-w_attract * distance_sq);
        swarming_cost += h_repellant * std::exp(-w_repellant * distance_sq);
    }
    return swarming_cost;
}

void MetaheuristicModel::optimization_loop()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    best_cost = std::numeric_limits<double>::infinity();
    best_bacterium = nullptr;

    colony.clear();
    for(int i = 0; i < bacteria_count; ++i)
    {
        auto bacterium = std::make_unique<Bacterium>(*maps, max_stations_per_cell, stations_powers);
        bacterium->initialize();
        bacterium->update(*evaluator);
        colony.push_back(std::move(bacterium));
    }

    for (int l = 0; l < elimination_dispersal_steps; ++l)
    {
        for (int k = 0; k < reproduction_steps; ++k)
        {
            for (int j = 0; j < hemotaxis_steps; ++j)
            {
                for (auto& bacterium : colony) {
                    bacterium->update(*evaluator);
                    double old_fitness = bacterium->current_cost + calculate_swarming_effect(*bacterium);
                    
                    bacterium->tumble();
                    bacterium->swim();
                    bacterium->update(*evaluator);
                    
                    double current_fitness = bacterium->current_cost + calculate_swarming_effect(*bacterium);
                    
                    int m = 0;
                    while (m < swimming_steps) {
                        if (bacterium->current_cost < best_cost) {
                            best_cost = bacterium->current_cost;
                            best_bacterium = std::make_unique<Bacterium>(*bacterium);
                        }

                        if (current_fitness >= old_fitness)
                            break;

                        old_fitness = current_fitness;
                        bacterium->swim();
                        bacterium->update(*evaluator);
                        current_fitness = bacterium->current_cost + calculate_swarming_effect(*bacterium);

                        m++;
                    }
                    bacterium->health += current_fitness;
                }

                const auto& solution = get_solution();
                printer->print_iteration(
                    solution,
                    validator->validate(solution),
                    l,
                    k,
                    j
                );
            }
            
            std::sort(colony.begin(), colony.end(), [](const auto& a, const auto& b) { 
                return a->health < b->health; 
            });
            
            int half_size = bacteria_count / 2;
            for (int i = 0; i < half_size; ++i)
            {
                colony[i + half_size] = std::make_unique<Bacterium>(*colony[i]);
                colony[i]->health = 0;
                colony[i + half_size]->health = 0;
            }
        }
        
        for (auto& bacterium : colony)
        {
            if (dis(gen) < elimination_dispersal_prob)
            {
                bacterium = std::make_unique<Bacterium>(*maps, max_stations_per_cell, stations_powers);
                bacterium->initialize();
                bacterium->update(*evaluator);
            }
        }
    }
}

Solution MetaheuristicModel::get_solution() const
{
    Solution solution;

    solution.station_location.assign(maps->get_width(), std::vector<int>(maps->get_height()));
    solution.l2_station_location.assign(maps->get_width(), std::vector<int>(maps->get_height()));
    solution.l3_station_location.assign(maps->get_width(), std::vector<int>(maps->get_height()));
    solution.demand_allocation_map.assign(
        maps->get_width(), 
        std::vector<std::vector<std::vector<double>>>(maps->get_height(), 
            std::vector<std::vector<double>>(maps->get_width(), 
                std::vector<double>(maps->get_height(), 0.0)
            )
        )
    );
    

    solution.total_cost = best_cost;
    if (best_bacterium) {
        for (int i = 0; i < maps->get_width(); ++i)
        {
            for (int j = 0; j < maps->get_height(); ++j)
            {
                solution.station_location[i][j] = best_bacterium->station_location[i][j];
                solution.l2_station_location[i][j] = best_bacterium->l2_station_location[i][j];
                solution.l3_station_location[i][j] = best_bacterium->l3_station_location[i][j];
                
                for (int m = 0; m < maps->get_width(); ++m) {
                    for (int n = 0; n < maps->get_height(); ++n) {
                        if (maps->get_demand_at(m, n) > 0 && !std::isinf(maps->get_distance(m, n, i, j))) {
                            solution.demand_allocation_map[m][n][i][j] = 
                                best_bacterium->demand_allocation_map[m][n][i][j];
                        }
                    }
                }
            }
        }
    }
    return solution;
}