#pragma once

#include <vector>
#include <utility>
#include "Maps.hpp"

enum class Map {
    STATION_LOCATION,
    L2_STATION_LOCATION,
    L3_STATION_LOCATION,
};

struct Move {
    Map map;
    std::vector<std::vector<int>> data;
};

struct Bacterium
{
    Bacterium(const Maps& maps, int max_stations_per_cell, std::pair<double, double> station_powers);
    Bacterium(const Bacterium& other);
    
    void tumble();
    void swim();
    void repair_and_validate_structure(double& penalty);
    void recalculate_demand_allocation();

    const Maps& maps;
    int max_stations_per_cell;
    std::pair<double, double> station_powers;

    Move current_move;
    
    std::vector<std::vector<int>> station_location;
    std::vector<std::vector<int>> l2_station_location;
    std::vector<std::vector<int>> l3_station_location;
    std::vector<std::vector<std::vector<std::vector<double>>>> demand_allocation_map;

    double current_cost = 0.0;
    double health = 0.0;
};