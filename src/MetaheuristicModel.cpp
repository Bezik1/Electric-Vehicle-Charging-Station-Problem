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
        double budget
) : max_stations_per_cell{max_stations_per_cell},
    bacteria_count{bacteria_count},
    hemotaxis_steps{hemotaxis_steps},
    swimming_steps{swimming_steps},
    reproduction_steps{reproduction_steps},
    elimination_dispersal_steps{elimination_dispersal_steps},
    elimination_dispersal_prob{elimination_dispersal_prob},
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

    evaluator = std::make_unique<Evaluator>(*maps, stations_powers, initial_costs, maintenance_costs, budget);

    optimization_loop();
}

double MetaheuristicModel::calculate_swarming_effect(const Bacterium& current_bacterium) const {
    double swarming_cost = 0.0;
    const double d_attract = 0.1;
    const double w_attract = 0.2;
    const double h_repellant = 0.1;
    const double w_repellant = 10.0;

    for (const auto& other : colony) {
        if (!other) continue;
        double distance_sq = 0.0;
        for (int i = 0; i < maps->get_width(); ++i) {
            for (int j = 0; j < maps->get_height(); ++j) {
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
    for(int i = 0; i < bacteria_count; ++i) {
        auto b = std::make_unique<Bacterium>(*maps, max_stations_per_cell, stations_powers);
        for (int x = 0; x < maps->get_width(); ++x) {
            for (int y = 0; y < maps->get_height(); ++y) {
                if (maps->get_poi_at(x, y) == 1 && dis(gen) < 0.4) {
                    b->station_location[x][y] = 1;

                    int start_l2 = std::min(1, max_stations_per_cell);
                    int start_l3 = std::min(1, max_stations_per_cell - start_l2);

                    b->l2_station_location[x][y] = start_l2;
                    b->l3_station_location[x][y] = start_l3;
                }
            }
        }
        
        double repair_penalty = 0.0;
        b->repair_and_validate_structure(repair_penalty);
        b->recalculate_demand_allocation();
        
        b->current_cost = evaluator->evaluate(
            b->station_location,
            b->l2_station_location,
            b->l3_station_location,
            b->demand_allocation_map) + repair_penalty;
        
        colony.push_back(std::move(b));
    }

    for (int l = 0; l < elimination_dispersal_steps; ++l) {
        for (int k = 0; k < reproduction_steps; ++k) {
            for (int j = 0; j < hemotaxis_steps; ++j) {
                for (auto& bacterium : colony) {
                    
                    double repair_p_old = 0.0;
                    bacterium->repair_and_validate_structure(repair_p_old);
                    bacterium->recalculate_demand_allocation();

                    double bacterium_fitness = evaluator->evaluate(
                        bacterium->station_location,
                        bacterium->l2_station_location,
                        bacterium->l3_station_location,
                        bacterium->demand_allocation_map) + repair_p_old;
                    
                    bacterium->current_cost = bacterium_fitness;
                    double old_fitness = bacterium_fitness + calculate_swarming_effect(*bacterium);
                    
                    bacterium->tumble();
                    bacterium->swim();

                    double repair_p_current = 0.0;
                    bacterium->repair_and_validate_structure(repair_p_current);
                    bacterium->recalculate_demand_allocation();

                    bacterium_fitness = evaluator->evaluate(
                        bacterium->station_location,
                        bacterium->l2_station_location,
                        bacterium->l3_station_location,
                        bacterium->demand_allocation_map) + repair_p_current;
                    
                    bacterium->current_cost = bacterium_fitness;
                    double current_fitness = bacterium_fitness + calculate_swarming_effect(*bacterium);

                    if (bacterium->current_cost < best_cost) {
                        best_cost = bacterium->current_cost;
                        best_bacterium = std::make_unique<Bacterium>(*bacterium);
                    }

                    int m = 0;
                    while (m < swimming_steps) {
                        if (current_fitness < old_fitness) {
                            old_fitness = current_fitness;
                            bacterium->swim();

                            double repair_p_swim = 0.0;
                            bacterium->repair_and_validate_structure(repair_p_swim);
                            bacterium->recalculate_demand_allocation();

                            bacterium_fitness = evaluator->evaluate(
                                bacterium->station_location,
                                bacterium->l2_station_location,
                                bacterium->l3_station_location,
                                bacterium->demand_allocation_map) + repair_p_swim;
                            
                            bacterium->current_cost = bacterium_fitness;
                            current_fitness = bacterium_fitness + calculate_swarming_effect(*bacterium);

                            if (bacterium->current_cost < best_cost) {
                                best_cost = bacterium->current_cost;
                                best_bacterium = std::make_unique<Bacterium>(*bacterium);
                            }
                            m++;
                        } else {
                            break;
                        }
                    }
                    bacterium->health += current_fitness;
                }
                
                std::cout << "Step L:" << l << " K:" << k << " J:" << j 
                          << " | Best Cost: " << best_cost 
                          << " | Valid: " << (best_bacterium ? is_valid(best_bacterium) : 0) << std::endl;
            }

            std::sort(colony.begin(), colony.end(), [](const auto& a, const auto& b) {
                return a->health < b->health;
            });

            int half_size = bacteria_count / 2;
            for (int i = 0; i < half_size; ++i) {
                colony[i + half_size] = std::make_unique<Bacterium>(*colony[i]);
                colony[i]->health = 0;
                colony[i + half_size]->health = 0;
            }
        }

        for (auto& bacterium : colony) {
            if (dis(gen) < elimination_dispersal_prob) {
                bacterium = std::make_unique<Bacterium>(*maps, max_stations_per_cell, stations_powers);
                for (int x = 0; x < maps->get_width(); ++x) {
                    for (int y = 0; y < maps->get_height(); ++y) {
                        if (maps->get_poi_at(x, y) == 1 && dis(gen) < 0.2) {
                            bacterium->station_location[x][y] = 1;
                            bacterium->l2_station_location[x][y] = std::min(1, max_stations_per_cell);
                        }
                    }
                }
                bacterium->recalculate_demand_allocation();
            }
        }
    }
}

bool MetaheuristicModel::is_valid(const std::unique_ptr<Bacterium>& bacterium) const {
    const double EPS = 1e-6;
    double total_investment = 0.0;

    for (int i = 0; i < maps->get_width(); ++i) {
        for (int j = 0; j < maps->get_height(); ++j) {
            int sL2 = bacterium->l2_station_location[i][j];
            int sL3 = bacterium->l3_station_location[i][j];
            int xi = bacterium->station_location[i][j];

            if (xi > 0 && maps->get_poi_at(i, j) == 0) return false;

            if ((sL2 + sL3) > (max_stations_per_cell * xi)) return false;

            total_investment += (sL2 * (initial_costs.first + maintenance_costs.first)) + 
                                (sL3 * (initial_costs.second + maintenance_costs.second)) + 
                                (xi * maps->get_rental_cost_at(i, j));
        }
    }

    if (total_investment > budget + EPS) return false;

    std::vector<std::vector<double>> station_load(maps->get_width(), std::vector<double>(maps->get_height(), 0.0));
    std::vector<std::vector<double>> demand_satisfied(maps->get_width(), std::vector<double>(maps->get_height(), 0.0));

    for (int m = 0; m < maps->get_width(); ++m) {
        for (int n = 0; n < maps->get_height(); ++n) {
            for (int i = 0; i < maps->get_width(); ++i) {
                for (int j = 0; j < maps->get_height(); ++j) {
                    double ratio = bacterium->demand_allocation_map[m][n][i][j];
                    if (ratio > 0) {
                        double amount = ratio * maps->get_demand_at(m, n);
                        station_load[i][j] += amount;
                        demand_satisfied[m][n] += amount;
                    }
                }
            }
        }
    }

    for (int i = 0; i < maps->get_width(); ++i) {
        for (int j = 0; j < maps->get_height(); ++j) {
            double max_cap = (bacterium->l2_station_location[i][j] * stations_powers.first +
                              bacterium->l3_station_location[i][j] * stations_powers.second) * 24.0;
            if (station_load[i][j] > max_cap + EPS) return false;
        }
    }

    for (int m = 0; m < maps->get_width(); ++m) {
        for (int n = 0; n < maps->get_height(); ++n) {
            if (maps->get_demand_at(m, n) > 0 && demand_satisfied[m][n] < maps->get_demand_at(m, n) - EPS) return false;
        }
    }

    return true;
}

void MetaheuristicModel::print_solution()
{
    std::cout << "Best Cost: " << best_cost << "\n";

    for (int n = 0; n < maps->get_width(); ++n)
    {
        for (int m = 0; m < maps->get_height(); ++m)
        {
            if (!maps->isPassable(n, m))
            {
                std::cout << " X ";
            }
            else
            {
                if (best_bacterium->station_location[n][m] > 0)
                {
                    int l2_count = best_bacterium->l2_station_location[n][m];
                    int l3_count = best_bacterium->l3_station_location[n][m];

                    if (l2_count > 0 && l3_count > 0)
                        std::cout << l2_count << ":" << l3_count;
                    else if (l2_count > 0)
                        std::cout << "2:" << l2_count;
                    else if (l3_count > 0)
                        std::cout << "3:" << l3_count;
                    else
                        std::cout << " 0 ";
                }
                else
                {
                    if (maps->get_poi_at(n, m) == 0)
                        std::cout << " - ";
                    else
                        std::cout << " . ";
                }
            }
            std::cout << " ";
        }
        std::cout << "\n";
    }
}

void MetaheuristicModel::print_demand_distribution()
{
    constexpr double EPSILON = 1e-4; 

    for (int n = 0; n < maps->get_width(); ++n)
    {
        for (int m = 0; m < maps->get_height(); ++m)
        {
            double current_demand = maps->get_demand_at(n, m);
            
            if (current_demand > 0.0)
            {
                std::cout << "Point (" << m << ", " << n << ") "
                          << "[Full-Demand: " << std::fixed << std::setprecision(1) << current_demand << " kWh]:\n";
                
                for (int i = 0; i < maps->get_width(); ++i)
                {
                    for (int j = 0; j < maps->get_height(); ++j)
                    {
                        if (!std::isinf(maps->get_distance(n, m, i, j)))
                        {
                            double allocation_ratio = best_bacterium->demand_allocation_map[n][m][i][j];

                            if (allocation_ratio > EPSILON)
                            {
                                double allocated_kwh = current_demand * allocation_ratio;
                                
                                std::cout << "  -> " << std::setprecision(1) << (allocation_ratio * 100.0) << "% demand "
                                          << "(" << allocated_kwh << " kWh) "
                                          << "goes to station: (" << i << ", " << j << ") "
                                          << "| distance: " << std::setprecision(2) << maps->get_distance(n, m, i, j) << "\n";
                            }
                        }
                    }
                }
                std::cout << "----------------------------------------------------------------------\n";
            }
        }
    }
    std::cout << "======================================================================\n\n";
}