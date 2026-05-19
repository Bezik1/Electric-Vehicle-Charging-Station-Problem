#include <vector>
#include <print>
#include <cmath>

#include "EVCSP.hpp"

int main() 
{
    double mean_station_energy_demand = 150.0;
    int num_stations = 5;
    int map_width = 20;
    int map_height = 20;
    
    std::vector<std::vector<int>> connectivity_map(map_width, std::vector<int>(map_height, 0));
    std::vector<std::vector<double>> activity_map(map_width, std::vector<double>(map_height, 10.0));
    std::vector<std::vector<double>> attractiveness_map(map_width, std::vector<double>(map_height, 0.2));
    std::vector<std::vector<double>> energy_capacity_map(map_width, std::vector<double>(map_height, 300.0));

    for (int x = 0; x < map_width; ++x) 
    {
        for (int y = 0; y < map_height; ++y) 
        {
            if (x == 10 || y == 10) 
            {
                activity_map[x][y] = 90.0;
                attractiveness_map[x][y] = 0.5;
                energy_capacity_map[x][y] = 400.0;
            }
            
            if (std::abs(x - 10) <= 3 && std::abs(y - 10) <= 3) 
            {
                activity_map[x][y] += 20.0;
                attractiveness_map[x][y] = 0.8;
            }

            if (x == 9 && y != 10) 
            {
                connectivity_map[x][y] = 1;
            }

            if (x < 4 && y < 4) 
            {
                energy_capacity_map[x][y] = 800.0;
                activity_map[x][y] = 5.0;
                attractiveness_map[x][y] = 0.1;
            }

            if (x > 15 && y > 15) 
            {
                energy_capacity_map[x][y] = 30.0;
            }
        }
    }

    attractiveness_map[10][4] = 1.0; activity_map[10][4] = 95.0;
    attractiveness_map[10][16] = 1.0; activity_map[10][16] = 95.0;
    attractiveness_map[4][10] = 1.0; activity_map[4][10] = 95.0;
    attractiveness_map[16][10] = 1.0; activity_map[16][10] = 95.0;

    EVCSP problem{mean_station_energy_demand, num_stations, map_width, map_height};

    std::printf("Starting EVCSP problem with a 20x20 Grid Roadmap!\n");

    problem(
        std::move(connectivity_map),
        std::move(activity_map),
        std::move(attractiveness_map),
        std::move(energy_capacity_map)
    );

    return 0;
}