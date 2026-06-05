#include "Bacterium.hpp"

#include <vector>    
#include <random>
#include <algorithm> 
#include <utility>   
#include <cmath>     

Bacterium::Bacterium(const Maps& maps, int max_stations_per_cell, std::pair<double, double> station_powers) 
    : maps(maps),
      max_stations_per_cell(max_stations_per_cell),
      station_powers(station_powers)
{
    int grid_width = maps.get_width();
    int grid_height = maps.get_height();

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

Bacterium::Bacterium(const Bacterium& other) 
    : maps(other.maps),
      max_stations_per_cell(other.max_stations_per_cell),
      station_powers(other.station_powers),
      current_cost(other.current_cost),
      health(other.health),
      current_move(other.current_move),
      station_location(other.station_location),
      l2_station_location(other.l2_station_location),
      l3_station_location(other.l3_station_location),
      demand_allocation_map(other.demand_allocation_map)
{}

void Bacterium::tumble()
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
    int width = maps.get_width();
    int height = maps.get_height();
    s_move.assign(width, std::vector<int>(height, 0));

    bool consolidation_move = (prob_dist(gen) < 0.2); 

    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            if (maps.get_poi_at(i, j) == 0) continue;
            
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

void Bacterium::swim()
{
    auto& update = current_move.data;
    int width = maps.get_width();
    int height = maps.get_height();
    int max_per_cell = max_stations_per_cell;

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
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
                int limit = station_location[i][j] * max_per_cell;

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

void Bacterium::recalculate_demand_allocation() {
    int width = maps.get_width();
    int height = maps.get_height();
    auto powers = station_powers;

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
            double demand = maps.get_demand_at(m, n);
            if (demand <= 0 || active_stations.empty()) continue;

            std::vector<std::pair<double, std::pair<int, int>>> sorted_stations;
            for (const auto& pos : active_stations) {
                double d = maps.get_distance(m, n, pos.first, pos.second);
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

void Bacterium::repair_and_validate_structure(double& penalty) {
    const double MISC_PENALTY = 1e6;
    int width = maps.get_width();
    int height = maps.get_height();
    int max_per_cell = max_stations_per_cell;

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            if (maps.get_poi_at(i, j) == 0) {
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
            int allowed = max_per_cell * station_location[i][j];
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