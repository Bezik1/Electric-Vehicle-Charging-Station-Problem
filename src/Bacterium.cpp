#include "Bacterium.hpp"

#include <random>
#include <vector>    
#include <algorithm> 
#include <utility>   
#include <cmath>     

Bacterium::Bacterium(const Maps& maps, int max_stations_per_cell, std::pair<double, double> station_powers) 
    : maps{maps},
      max_stations_per_cell{max_stations_per_cell},
      station_powers{station_powers}
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

void Bacterium::initialize()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int x = 0; x < maps.get_width(); ++x) {
        for (int y = 0; y < maps.get_height(); ++y) {
            if (maps.get_poi_at(x, y) == 1 && dis(gen) < INITIALIZATION_PROB) {
                l2_station_location[x][y] = 1;
                station_location[x][y] = 1;
            }
        }
    }
}

void Bacterium::tumble()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> prob_dist(0.0, 1.0);
    
    int width = maps.get_width();
    int height = maps.get_height();
    std::vector<std::vector<int>> move(width, std::vector<int>(height, 0));

    if (this->current_cost > HARD_THRESHOLD)
    {
        std::pair<int, int> worst_station = {-1, -1};
        double min_load = std::numeric_limits<double>::max();
        bool found_active = false;

        for (int i = 0; i < width; ++i)
        {
            for (int j = 0; j < height; ++j)
            {
                if (station_location[i][j] > 0)
                {
                    double total_load = 0;
                    for (int m = 0; m < width; ++m)
                        for (int n = 0; n < height; ++n) 
                            total_load += demand_allocation_map[m][n][i][j] * maps.get_demand_at(m, n);

                    if (total_load < min_load) {
                        min_load = total_load;
                        worst_station = {i, j};
                        found_active = true;
                    }
                }
            }
        }

        if (found_active)
        {
            this->current_move.map = l3_station_location[worst_station.first][worst_station.second] > 0
                ? Map::L3_STATION_LOCATION
                : Map::L2_STATION_LOCATION;

            move[worst_station.first][worst_station.second] = -1;
            this->current_move.move = std::move(move);
            return;
        }
    }

    std::uniform_int_distribution<> action_dist(1, 2);
    this->current_move.map = static_cast<Map>(action_dist(gen));
    
    std::uniform_int_distribution<> move_dist(-1, 1);

    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            if (maps.get_poi_at(i, j) == 0 || prob_dist(gen) >= CHANGE_PROB)
                continue;
            
            move[i][j] = move_dist(gen);
        }
    }
    this->current_move.move = std::move(move);
}

void Bacterium::swim()
{
    auto& update_move = current_move.move;
    int width = maps.get_width();
    int height = maps.get_height();

    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            if (update_move[i][j] == 0)
                continue;

            switch (current_move.map)
            {
            case Map::L2_STATION_LOCATION:
                l2_station_location[i][j] = std::max(0, l2_station_location[i][j] + update_move[i][j]);
                break;
            case Map::L3_STATION_LOCATION:
                l3_station_location[i][j] = std::max(0, l3_station_location[i][j] + update_move[i][j]);
                break;
            }

            station_location[i][j] = l2_station_location[i][j] + l3_station_location[i][j] > 0
                ? 1
                : 0;
        }
    }
}

void Bacterium::update(const Evaluator& evaluator)
{
    double repair_penalty = repair_and_validate_structure();
    recalculate_demand_allocation();
    current_cost = evaluator.evaluate(
        station_location,
        l2_station_location,
        l3_station_location,
        demand_allocation_map
    ) + repair_penalty;
}

void Bacterium::recalculate_demand_allocation()
{
    int width = maps.get_width();
    int height = maps.get_height();
    auto powers = station_powers;

    for (int m = 0; m < width; ++m)
        for (int n = 0; n < height; ++n)
            for (int i = 0; i < width; ++i)
                std::fill(demand_allocation_map[m][n][i].begin(), demand_allocation_map[m][n][i].end(), 0.0);

    std::vector<std::vector<double>> station_capacities(width, std::vector<double>(height, 0.0));
    std::vector<std::pair<int, int>> active_stations;
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            int total_chargers = l2_station_location[i][j] + l3_station_location[i][j];
            if (total_chargers > 0) {
                station_capacities[i][j] = (l2_station_location[i][j] * powers.first + 
                                           l3_station_location[i][j] * powers.second) * 24.0;
                active_stations.push_back({i, j});
            }
        }
    }
    std::vector<std::vector<double>> station_load(width, std::vector<double>(height, 0.0));
    for (int m = 0; m < width; ++m)
    {
        for (int n = 0; n < height; ++n)
        {
            double demand = maps.get_demand_at(m, n);

            if (demand <= 0 || active_stations.empty())
                continue;

            std::vector<std::pair<double, std::pair<int, int>>> sorted_stations;
            for (const auto& pos : active_stations)
            {
                double distance = maps.get_distance(m, n, pos.first, pos.second);
                if (!std::isinf(distance)) sorted_stations.push_back({distance, pos});
            }

            if (sorted_stations.empty())
                continue;

            std::sort(sorted_stations.begin(), sorted_stations.end());

            double remaining_demand = demand;
            for (const auto& entry : sorted_stations)
            {
                int si = entry.second.first;
                int sj = entry.second.second;
                double capacity_left = std::max(0.0, station_capacities[si][sj] - station_load[si][sj]);
                double to_allocate = std::min(remaining_demand, capacity_left);

                if (to_allocate > EPS) {
                    demand_allocation_map[m][n][si][sj] = to_allocate / demand;
                    station_load[si][sj] += to_allocate;
                    remaining_demand -= to_allocate;
                }

                if (remaining_demand <= EPS)
                    break;
            }
        }
    }
}

double Bacterium::repair_and_validate_structure()
{
    double penalty = 0.0;
    int width = maps.get_width();
    int height = maps.get_height();

    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            if (maps.get_poi_at(i, j) == 0)
            {
                if (station_location[i][j] > 0 || l2_station_location[i][j] > 0 || l3_station_location[i][j] > 0)
                {
                    penalty += MISC_PENALTY;
                    station_location[i][j] = 0;
                    l2_station_location[i][j] = 0;
                    l3_station_location[i][j] = 0;
                }
                continue;
            }

            station_location[i][j] = (l2_station_location[i][j] + l3_station_location[i][j] > 0) 
                ? 1 
                : 0;
            
            int allowed = station_location[i][j] * max_stations_per_cell;
            int current_total = l2_station_location[i][j] + l3_station_location[i][j];
            if (current_total > allowed)
            {
                penalty += MISC_PENALTY * (current_total - allowed);

                while (l2_station_location[i][j] + l3_station_location[i][j] > allowed)
                {
                    if (l3_station_location[i][j] > 0) l3_station_location[i][j]--;
                    else if (l2_station_location[i][j] > 0) l2_station_location[i][j]--;
                    else break;
                }
            }
        }
    }

    return penalty;
}