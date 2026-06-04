#include "MetaheuristicModel.hpp"

#include <queue>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

MetaheuristicModel::MetaheuristicModel(
        int grid_width,
        int grid_height,
        int max_stations_per_cell,
        int bacteria_count,
        int hemotaxis_steps,
        int swimming_steps,
        int reproduction_steps,
        int elimination_dispersal_steps,
        double elimination_dispersal_prob,
        double budget
) : grid_width{grid_width},
    grid_height{grid_height},
    max_stations_per_cell{max_stations_per_cell},
    bacteria_count{bacteria_count},
    hemotaxis_steps{hemotaxis_steps},
    swimming_steps{swimming_steps},
    reproduction_steps{reproduction_steps},
    elimination_dispersal_steps{elimination_dispersal_steps},
    elimination_dispersal_prob{elimination_dispersal_prob},
    budget{budget},
    evaluator{std::make_unique<Evaluator>(*this)}
{}

void MetaheuristicModel::operator()(
            std::vector<std::vector<double>> _distances_costs_map,
            std::vector<std::vector<int>> _poi_map,
            std::vector<std::vector<double>> _demand_map,
            std::vector<std::vector<double>> _land_rental_cost_map,
            std::pair<double, double> _stations_powers,
            std::pair<double, double> _initial_costs,
            std::pair<double, double> _maintenance_costs
        ) 
{
    distances_costs_map = std::move(_distances_costs_map);
    poi_map = std::move(_poi_map);
    demand_map = std::move(_demand_map);
    land_rental_cost_map = std::move(_land_rental_cost_map);
    stations_powers = std::move(_stations_powers);
    initial_costs = std::move(_initial_costs);
    maintenance_costs = std::move(_maintenance_costs);

    distances_from_to.assign(grid_width, std::vector<std::vector<std::vector<double>>>(
        grid_height, 
        std::vector<std::vector<double>>(grid_width, std::vector<double>(
                                                        grid_height, std::numeric_limits<double>::infinity()))
    ));

    for (int m = 0; m < grid_width; ++m)
    {
        for (int n = 0; n < grid_height; ++n)
        {
            if (demand_map[m][n] <= 0) continue;

            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    distances_from_to[m][n][i][j] = calculate_distance(m, n, i, j);
                }
            }
        }
    }

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
        for (int i = 0; i < grid_width; ++i) {
            for (int j = 0; j < grid_height; ++j) {
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

void MetaheuristicModel::Bacterium::repair_and_validate_structure(const MetaheuristicModel& model, double& penalty) {
    const double MISC_PENALTY = 1e6;
    for (int i = 0; i < model.grid_width; ++i) {
        for (int j = 0; j < model.grid_height; ++j) {
            if (model.poi_map[i][j] == 0) {
                station_location[i][j] = 0;
                l2_station_location[i][j] = 0;
                l3_station_location[i][j] = 0;
                continue;
            }
            int current_l2 = l2_station_location[i][j];
            int current_l3 = l3_station_location[i][j];
            int total_chargers = current_l2 + current_l3;
            
            if (total_chargers > 0 && station_location[i][j] == 0) {
                station_location[i][j] = 1; 
            }
            int allowed = model.max_stations_per_cell * station_location[i][j];
            if (total_chargers > allowed) {
                penalty += MISC_PENALTY * (total_chargers - allowed);
                while (l2_station_location[i][j] + l3_station_location[i][j] > allowed) {
                    if (l2_station_location[i][j] > 0) l2_station_location[i][j]--;
                    else if (l3_station_location[i][j] > 0) l3_station_location[i][j]--;
                }
            }
        }
    }
}

void MetaheuristicModel::Bacterium::recalculate_demand_allocation(const MetaheuristicModel& model) {
    int width = model.getWidth();
    int height = model.getHeight();
    auto powers = model.getStationPowers();

    for (int m = 0; m < width; ++m) {
        for (int n = 0; n < height; ++n) {
            for (int i = 0; i < width; ++i) {
                std::fill(demand_allocation_map[m][n][i].begin(), demand_allocation_map[m][n][i].end(), 0.0);
            }
        }
    }

    std::vector<std::vector<double>> station_capacities(width, std::vector<double>(height, 0.0));
    std::vector<std::pair<int, int>> active_stations;

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            int total_chargers = l2_station_location[i][j] + l3_station_location[i][j];
            if (total_chargers > 0) {
                station_capacities[i][j] = (l2_station_location[i][j] * powers.first + 
                                           l3_station_location[i][j] * powers.second) * 24.0;
                active_stations.push_back({i, j});
            }
        }
    }

    std::vector<std::vector<double>> station_load(width, std::vector<double>(height, 0.0));

    for (int m = 0; m < width; ++m) {
        for (int n = 0; n < height; ++n) {
            double demand = model.getDemandAt(m, n);
            if (demand <= 0 || active_stations.empty()) continue;

            std::vector<std::pair<double, std::pair<int, int>>> sorted_stations;
            for (const auto& pos : active_stations) {
                double d = model.getDistance(m, n, pos.first, pos.second);
                if (!std::isinf(d)) sorted_stations.push_back({d, pos});
            }

            if (sorted_stations.empty()) continue;
            std::sort(sorted_stations.begin(), sorted_stations.end());

            double remaining_demand = demand;
            for (const auto& entry : sorted_stations) {
                int si = entry.second.first;
                int sj = entry.second.second;
                double capacity_left = std::max(0.0, station_capacities[si][sj] - station_load[si][sj]);
                double to_allocate = std::min(remaining_demand, capacity_left);

                if (to_allocate > 1e-9) {
                    demand_allocation_map[m][n][si][sj] = to_allocate / demand;
                    station_load[si][sj] += to_allocate;
                    remaining_demand -= to_allocate;
                }
                if (remaining_demand <= 1e-9) break;
            }
        }
    }
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
        auto b = std::make_unique<MetaheuristicModel::Bacterium>(grid_width, grid_height);
        for (int x = 0; x < grid_width; ++x) {
            for (int y = 0; y < grid_height; ++y) {
                if (poi_map[x][y] == 1 && dis(gen) < 0.4) {
                    b->station_location[x][y] = 1;
                    int start_l2 = std::min(1, max_stations_per_cell);
                    int start_l3 = std::min(1, max_stations_per_cell - start_l2);
                    b->l2_station_location[x][y] = start_l2;
                    b->l3_station_location[x][y] = start_l3;
                }
            }
        }
        
        double repair_penalty = 0.0;
        b->repair_and_validate_structure(*this, repair_penalty);
        
        b->recalculate_demand_allocation(*this);
        
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
                    bacterium->repair_and_validate_structure(*this, repair_p_old);
                    bacterium->recalculate_demand_allocation(*this);

                    double bacterium_fitness = evaluator->evaluate(
                        bacterium->station_location,
                        bacterium->l2_station_location,
                        bacterium->l3_station_location,
                        bacterium->demand_allocation_map) + repair_p_old;
                    
                    bacterium->current_cost = bacterium_fitness;
                    double old_fitness = bacterium_fitness + calculate_swarming_effect(*bacterium);
                    
                    bacterium->tumble(*this);
                    bacterium->swim(*this);

                    double repair_p_current = 0.0;
                    bacterium->repair_and_validate_structure(*this, repair_p_current);
                    
                    bacterium->recalculate_demand_allocation(*this);

                    bacterium_fitness = evaluator->evaluate(
                        bacterium->station_location,
                        bacterium->l2_station_location,
                        bacterium->l3_station_location,
                        bacterium->demand_allocation_map) + repair_p_current;
                    
                    bacterium->current_cost = bacterium_fitness;
                    double current_fitness = bacterium_fitness + calculate_swarming_effect(*bacterium);

                    if (bacterium->current_cost < best_cost) {
                        best_cost = bacterium->current_cost;
                        best_bacterium = std::make_unique<MetaheuristicModel::Bacterium>(*bacterium);
                    }

                    int m = 0;
                    while (m < swimming_steps) {
                        if (current_fitness < old_fitness) {
                            old_fitness = current_fitness;
                            bacterium->swim(*this);

                            double repair_p_swim = 0.0;
                            bacterium->repair_and_validate_structure(*this, repair_p_swim);
                            bacterium->recalculate_demand_allocation(*this);

                            bacterium_fitness = evaluator->evaluate(
                                bacterium->station_location,
                                bacterium->l2_station_location,
                                bacterium->l3_station_location,
                                bacterium->demand_allocation_map) + repair_p_swim;
                            
                            bacterium->current_cost = bacterium_fitness;
                            current_fitness = bacterium_fitness + calculate_swarming_effect(*bacterium);

                            if (bacterium->current_cost < best_cost) {
                                best_cost = bacterium->current_cost;
                                best_bacterium = std::make_unique<MetaheuristicModel::Bacterium>(*bacterium);
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
                colony[i + half_size] = std::make_unique<MetaheuristicModel::Bacterium>(*colony[i]);
                colony[i]->health = 0;
                colony[i + half_size]->health = 0;
            }
        }

        for (auto& bacterium : colony) {
            if (dis(gen) < elimination_dispersal_prob) {
                bacterium = std::make_unique<MetaheuristicModel::Bacterium>(grid_width, grid_height);
                for (int x = 0; x < grid_width; ++x) {
                    for (int y = 0; y < grid_height; ++y) {
                        if (poi_map[x][y] == 1 && dis(gen) < 0.2) {
                            bacterium->station_location[x][y] = 1;
                            bacterium->l2_station_location[x][y] = std::min(1, max_stations_per_cell);
                        }
                    }
                }
                bacterium->recalculate_demand_allocation(*this);
            }
        }
    }
}

std::pair<double, std::pair<int, int>> MetaheuristicModel::find_nearest_valid(int start_x, int start_y) const
{
    double distance_to_nearest_valid = 0;

    if (distances_costs_map[start_x][start_y] != INF_VAL)
        return { distance_to_nearest_valid, {start_x, start_y} };

    std::queue<std::pair<int, int>> q;
    std::vector<std::vector<bool>> visited(grid_width, std::vector<bool>(grid_height, false));

    q.push({start_x, start_y});
    visited[start_x][start_y] = true;

    while(!q.empty())
    {
        distance_to_nearest_valid += 1;
        auto [cx, cy] = q.front();
        q.pop();

        for(const auto& [dx, dy] : DIRECTIONS)
        {
            int nx = cx + dx;
            int ny = cy + dy;

            if(nx >= 0 && nx < grid_width && ny >= 0 && ny < grid_height && !visited[nx][ny])
            {
                if (distances_costs_map[nx][ny] != INF_VAL)
                    return { distance_to_nearest_valid, {nx, ny} };

                visited[nx][ny] = true;
                q.push({nx, ny});
            }
        }
    }
    return { distance_to_nearest_valid, {-1, -1} };
}

double MetaheuristicModel::calculate_distance(int m, int n, int i, int j) const
{
    using Node = std::pair<double, std::pair<int, int>>;

    if (m == i && n == j)
        return 0.0;

    if (distances_costs_map[i][j] == INF_VAL)
        return std::numeric_limits<double>::infinity();

    auto [extra_d, actual_start] = find_nearest_valid(m, n);
    if (actual_start.first == -1)
        return std::numeric_limits<double>::infinity();

    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> queue;
    std::vector<std::vector<double>> dist(
        grid_width, 
        std::vector<double>(grid_height, std::numeric_limits<double>::infinity())
    );

    dist[actual_start.first][actual_start.second] = extra_d;
    queue.push({extra_d, actual_start});

    while (!queue.empty())
    {
        auto [current_d, coord] = queue.top();
        auto [x, y] = coord;
        queue.pop();

        if (current_d > dist[x][y]) continue;
        if (x == i && y == j) return current_d;

        for (const auto& [dx, dy] : DIRECTIONS)
        {
            int nx = x + dx; int ny = y + dy;

            if (nx < 0 || nx >= grid_width || ny < 0 || ny >= grid_height || distances_costs_map[nx][ny] == INF_VAL)
                continue;

            double weight = (distances_costs_map[nx][ny] <= 0) ? 1.0 : distances_costs_map[nx][ny];
            if (current_d + weight < dist[nx][ny]) {
                dist[nx][ny] = current_d + weight;
                queue.push({dist[nx][ny], {nx, ny}});
            }
        }
    }
    return std::numeric_limits<double>::infinity();
}

bool MetaheuristicModel::is_valid(const std::unique_ptr<MetaheuristicModel::Bacterium>& bacterium) const {
    const double EPS = 1e-6;
    double total_investment = 0.0;

    for (int i = 0; i < grid_width; ++i) {
        for (int j = 0; j < grid_height; ++j) {
            int sL2 = bacterium->l2_station_location[i][j];
            int sL3 = bacterium->l3_station_location[i][j];
            int xi = bacterium->station_location[i][j];

            if (xi > 0 && poi_map[i][j] == 0) return false;

            if ((sL2 + sL3) > (max_stations_per_cell * xi)) return false;

            total_investment += (sL2 * (initial_costs.first + maintenance_costs.first)) + 
                                (sL3 * (initial_costs.second + maintenance_costs.second)) + 
                                (xi * land_rental_cost_map[i][j]);
        }
    }

    if (total_investment > budget + EPS) return false;

    std::vector<std::vector<double>> station_load(grid_width, std::vector<double>(grid_height, 0.0));
    std::vector<std::vector<double>> demand_satisfied(grid_width, std::vector<double>(grid_height, 0.0));

    for (int m = 0; m < grid_width; ++m) {
        for (int n = 0; n < grid_height; ++n) {
            for (int i = 0; i < grid_width; ++i) {
                for (int j = 0; j < grid_height; ++j) {
                    double ratio = bacterium->demand_allocation_map[m][n][i][j];
                    if (ratio > 0) {
                        double amount = ratio * demand_map[m][n];
                        station_load[i][j] += amount;
                        demand_satisfied[m][n] += amount;
                    }
                }
            }
        }
    }

    for (int i = 0; i < grid_width; ++i) {
        for (int j = 0; j < grid_height; ++j) {
            double max_cap = (bacterium->l2_station_location[i][j] * stations_powers.first +
                              bacterium->l3_station_location[i][j] * stations_powers.second) * 24.0;
            if (station_load[i][j] > max_cap + EPS) return false;
        }
    }

    for (int m = 0; m < grid_width; ++m) {
        for (int n = 0; n < grid_height; ++n) {
            if (demand_map[m][n] > 0 && demand_satisfied[m][n] < demand_map[m][n] - EPS) return false;
        }
    }

    return true;
}

void MetaheuristicModel::print_solution()
{
    std::cout << "Best Cost: " << best_cost << "\n";

    for (int n = 0; n < grid_width; ++n)
    {
        for (int m = 0; m < grid_height; ++m)
        {
            if (distances_costs_map[n][m] == INF_VAL)
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
                    if (poi_map[n][m] == 0)
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

    for (int n = 0; n < grid_width; ++n)
    {
        for (int m = 0; m < grid_height; ++m)
        {
            double current_demand = demand_map[n][m];
            
            if (current_demand > 0.0)
            {
                std::cout << "Point (" << m << ", " << n << ") "
                          << "[Full-Demand: " << std::fixed << std::setprecision(1) << current_demand << " kWh]:\n";
                
                for (int i = 0; i < grid_width; ++i)
                {
                    for (int j = 0; j < grid_height; ++j)
                    {
                        if (!std::isinf(distances_from_to[n][m][i][j]))
                        {
                            double allocation_ratio = best_bacterium->demand_allocation_map[n][m][i][j];

                            if (allocation_ratio > EPSILON)
                            {
                                double allocated_kwh = current_demand * allocation_ratio;
                                
                                std::cout << "  -> " << std::setprecision(1) << (allocation_ratio * 100.0) << "% demand "
                                          << "(" << allocated_kwh << " kWh) "
                                          << "goes to station: (" << i << ", " << j << ") "
                                          << "| distance: " << std::setprecision(2) << distances_from_to[n][m][i][j] << "\n";
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

MetaheuristicModel::Bacterium::Bacterium(int grid_width, int grid_height)
{
    station_location.assign(grid_width, std::vector<int>(grid_height, 0));
    l2_station_location.assign(grid_width, std::vector<int>(grid_height, 0));
    l3_station_location.assign(grid_width, std::vector<int>(grid_height, 0));

    demand_allocation_map.assign(grid_width, std::vector<std::vector<std::vector<double>>>(
        grid_height, std::vector<std::vector<double>>(
            grid_width, std::vector<double>(grid_height, 0.0)
        )
    ));
    
    current_cost = 0.0;
    health = 0.0;
}

MetaheuristicModel::Bacterium::Bacterium(const MetaheuristicModel::Bacterium& other) 
    : current_cost(other.current_cost),
      health(other.health),
      current_move(other.current_move),
      station_location(other.station_location),
      l2_station_location(other.l2_station_location),
      l3_station_location(other.l3_station_location),
      demand_allocation_map(other.demand_allocation_map)
{}

void MetaheuristicModel::Bacterium::tumble(const MetaheuristicModel& model)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> prob_dist(0.0, 1.0);
    
    std::uniform_int_distribution<> action_dist(0, 2);
    int action = action_dist(gen);
    
    const double change_probability = 0.15; 
    std::uniform_int_distribution<> move_dist(-2, 2); 

    this->current_move.map = static_cast<Map>(action);
    std::vector<std::vector<int>> s_move;
    s_move.assign(model.grid_width, std::vector<int>(model.grid_height, 0));

    bool consolidation_move = (prob_dist(gen) < 0.2); 

    for (int i = 0; i < model.grid_width; ++i)
    {
        for (int j = 0; j < model.grid_height; ++j)
        {
            if (model.poi_map[i][j] == 0) continue;
            
            if (consolidation_move) {
                if (station_location[i][j] > 0 && (l2_station_location[i][j] + l3_station_location[i][j] < 2)) {
                    if (prob_dist(gen) < 0.3) s_move[i][j] = -1;
                }
            } else if (prob_dist(gen) < change_probability) {
                s_move[i][j] = move_dist(gen);
            }
        }
    }
    this->current_move.data = std::move(s_move);
}

void MetaheuristicModel::Bacterium::swim(const MetaheuristicModel& model)
{
    auto& update = current_move.data;
    for (int i = 0; i < model.grid_width; ++i) {
        for (int j = 0; j < model.grid_height; ++j) {
            if (update[i][j] == 0) continue;

            if (current_move.map == Map::STATION_LOCATION) {
                station_location[i][j] = std::max(0, station_location[i][j] + update[i][j]);
                if (station_location[i][j] == 0) {
                    l2_station_location[i][j] = 0;
                    l3_station_location[i][j] = 0;
                }
            } else {
                int current_l2 = l2_station_location[i][j];
                int current_l3 = l3_station_location[i][j];
                int limit = station_location[i][j] * model.max_stations_per_cell;

                if (current_move.map == Map::L2_STATION_LOCATION) {
                    int possible_l2 = std::max(0, current_l2 + update[i][j]);
                    if (possible_l2 + current_l3 > limit) {
                        l2_station_location[i][j] = std::max(0, limit - current_l3);
                    } else {
                        l2_station_location[i][j] = possible_l2;
                    }
                } else if (current_move.map == Map::L3_STATION_LOCATION) {
                    int possible_l3 = std::max(0, current_l3 + update[i][j]);
                    if (possible_l3 + current_l2 > limit) {
                        l3_station_location[i][j] = std::max(0, limit - current_l2);
                    } else {
                        l3_station_location[i][j] = possible_l3;
                    }
                }
                
                if ((l2_station_location[i][j] > 0 || l3_station_location[i][j] > 0) && station_location[i][j] == 0) {
                    station_location[i][j] = 1;
                }
            }
        }
    }
}