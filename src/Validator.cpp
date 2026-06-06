#include "Validator.hpp"

Validator::Validator(
        const Maps &maps,
        int max_stations_per_cell,
        std::pair<double, double> stations_powers,
        std::pair<double, double> initial_costs,
        std::pair<double, double> maintenance_costs,
        double budget
) : maps{maps},
    max_stations_per_cell{max_stations_per_cell},
    stations_powers{std::move(stations_powers)},
    initial_costs{std::move(initial_costs)},
    maintenance_costs{std::move(maintenance_costs)},
    budget{budget}
{}

bool Validator::validate(const Solution& solution) const {
    double total_investment = 0.0;
    int width = maps.get_width();
    int height = maps.get_height();

    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            int sL2 = solution.l2_station_location[i][j];
            int sL3 = solution.l3_station_location[i][j];
            int xi = solution.station_location[i][j];

            if (xi > 0 && maps.get_poi_at(i, j) == 0)
                return false;
            
            if (xi > 1)
                return false;
            
            if ((sL2 + sL3) > (max_stations_per_cell * xi))
                return false;
            

            total_investment += (sL2 * (initial_costs.first + maintenance_costs.first)) + 
                                (sL3 * (initial_costs.second + maintenance_costs.second)) + 
                                (xi * maps.get_rental_cost_at(i, j));
        }
    }
    if (total_investment > budget + EPS)
        return false;

    std::vector<std::vector<double>> station_load(width, std::vector<double>(height, 0.0));
    std::vector<std::vector<double>> demand_satisfied(width, std::vector<double>(height, 0.0));
    for (int m = 0; m < width; ++m)
    {
        for (int n = 0; n < height; ++n)
        {
            double demand = maps.get_demand_at(m, n);
            if (demand <= 0) continue;
            
            for (int i = 0; i < width; ++i)
            {
                for (int j = 0; j < height; ++j)
                {
                    double ratio = solution.demand_allocation_map[m][n][i][j];
                    if (ratio > 0) {
                        double amount = ratio * demand;
                        station_load[i][j] += amount;
                        demand_satisfied[m][n] += amount;
                    }
                }
            }
        }
    }
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            double max_cap = (solution.l2_station_location[i][j] * stations_powers.first +
                              solution.l3_station_location[i][j] * stations_powers.second) * 24.0;

            if (station_load[i][j] > max_cap + EPS)
                return false;
        }
    }
    return true;
}
