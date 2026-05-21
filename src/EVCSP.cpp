#include "EVCSP.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <queue>

EVCSP::EVCSP(
            int _grid_width,
            int _grid_height,
            int _max_stations_per_cell,
            double _budget
) : grid_width{_grid_width},
    grid_height{_grid_height},
    max_stations_per_cell{_max_stations_per_cell},
    budget{_budget},
    env{},
    model{env}
{
    model.set(GRB_StringAttr_ModelName, MODEL_NAME);
}

void EVCSP::operator()(
            std::vector<std::vector<double>>&& _distances_costs_map,
            std::vector<std::vector<int>>&& _poi_map,
            std::vector<std::vector<double>>&& _demand_map,
            std::vector<std::vector<double>>&& _land_rental_cost_map,
            std::pair<double, double>&& _stations_powers,
            std::pair<double, double>&& _initial_costs,
            std::pair<double, double>&& _maintenance_costs
)
{
    distances_costs_map = std::move(_distances_costs_map);
    poi_map = std::move(_poi_map);
    demand_map = std::move(_demand_map);
    land_rental_cost_map = std::move(_land_rental_cost_map);
    stations_powers = std::move(_stations_powers);
    initial_costs = std::move(_initial_costs);
    maintenance_costs = std::move(_maintenance_costs);

    demand_allocation_map.resize(grid_width);
    for (int m = 0; m < grid_width; ++m)
    {
        demand_allocation_map[m].resize(grid_height);
        for (int n = 0; n < grid_height; ++n)
        {
            demand_allocation_map[m][n].resize(grid_width);
            for (int i = 0; i < grid_width; ++i)
            {
                demand_allocation_map[m][n][i].resize(grid_height);
            }
        }
    }

    distances_from_to.resize(grid_width);
    for (int m = 0; m < grid_width; ++m)
    {
        distances_from_to[m].resize(grid_height);
        for (int n = 0; n < grid_height; ++n)
        {
            distances_from_to[m][n].resize(grid_width);
            for (int i = 0; i < grid_width; ++i)
            {
                distances_from_to[m][n][i].resize(grid_height);
                for (int j = 0; j < grid_height; ++j)
                {
                    distances_from_to[m][n][i][j] = calculate_distance(m, n, i, j);
                }
            }
        }
    }

    build_variables();
    build_constraints();
    build_criterion();

    model.optimize();
}

double EVCSP::calculate_distance(int m, int n, int i, int j) const
{
    constexpr double INF = std::numeric_limits<double>::infinity();

    if (m == i && n == j) return 0.0;

    using Node = std::pair<double, std::pair<int, int>>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> queue;
    std::vector<std::vector<double>> dist(grid_width, std::vector<double>(grid_height, INF));

    dist[m][n] = 0.0;
    queue.push({0.0, {m, n}});

    while (!queue.empty())
    {
        auto& [current_d, coord] = queue.top();
        auto& [x, y] = coord;
        queue.pop();

        if (current_d > dist[x][y])
            continue;

        if (x == i && y == j)
            return dist[x][y];

        for (const auto& [dx, dy] : DIRECTIONS)
        {
            int next_x = x + dx;
            int next_y = y + dy;

            if (next_x < 0 || next_x >= grid_width || next_y < 0 || next_y >= grid_height)
                continue;

            if (distances_costs_map[next_x][next_y] == INF)
                continue;

            double new_dist = dist[x][y] + distances_costs_map[next_x][next_y];

            if (new_dist < dist[next_x][next_y])
            {
                dist[next_x][next_y] = new_dist;
                queue.push({new_dist, {next_x, next_y}});
            }
        }
    }
    return INF;
}

void EVCSP::build_variables()
{
    station_location.resize(grid_width);
    l2_station_location.resize(grid_width);
    l3_station_location.resize(grid_width);

    for (int m = 0; m < grid_width; ++m)
    {
        station_location[m].resize(grid_height);
        l2_station_location[m].resize(grid_height);
        l3_station_location[m].resize(grid_height);

        for (int n = 0; n < grid_height; ++n)
        {
            station_location[m][n] = model.addVar(
                0.0,
                1.0,
                0.0,
                GRB_BINARY
            );
            
            l2_station_location[m][n] = model.addVar(
                0.0,
                max_stations_per_cell,
                0.0,
                GRB_INTEGER
            );

            l3_station_location[m][n] = model.addVar(
                0.0,
                max_stations_per_cell,
                0.0,
                GRB_INTEGER
            );

            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    demand_allocation_map[m][n][i][j] = model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_CONTINUOUS);
                }
            }
        }
    }
    model.update();
}

void EVCSP::build_constraints()
{
    GRBLinExpr total_investment_cost = 0.0;

    for (int i = 0; i < grid_width; ++i)
    {
        for (int j = 0; j < grid_height; ++j)
        {
            // Station Construction Constraint
            model.addConstr(
                l2_station_location[i][j] + l3_station_location[i][j] <= 
                    max_stations_per_cell * station_location[i][j]);
            
            // POI Constraint
            if(poi_map[i][j] == 0)
                model.addConstr(station_location[i][j] == 0);

            total_investment_cost += (initial_costs.first * l2_station_location[i][j] + 
                                      initial_costs.second * l3_station_location[i][j] + 
                                      land_rental_cost_map[i][j] * station_location[i][j]);
        }
    }
    
    // Budget Constraint
    model.addConstr(total_investment_cost <= budget);

    for (int m = 0; m < grid_width; ++m)
    {
        for (int n = 0; n < grid_height; ++n)
        {
            GRBLinExpr sum = 0.0;
            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    // Obstacle Constraint
                    if(std::isinf(distances_from_to[m][n][i][j]))
                        model.addConstr(demand_allocation_map[m][n][i][j] == 0.0);

                    sum += demand_allocation_map[m][n][i][j];
                }
            }
            // Demand Satisfaction Constraint
            model.addConstr(sum == 1.0);
        }
    }

    for (int i = 0; i < grid_width; ++i)
    {
        for (int j = 0; j < grid_height; ++j)
        {
            GRBLinExpr total_demand_at_station = 0.0;
            for (int m = 0; m < grid_width; ++m)
            {
                for (int n = 0; n < grid_height; ++n)
                {
                    total_demand_at_station += demand_map[m][n] * demand_allocation_map[m][n][i][j];
                }
            }

            // Capacity Constraint
            model.addConstr(
                total_demand_at_station <=
                    (l2_station_location[i][j] * stations_powers.first * 24 * 0.4) +
                    (l3_station_location[i][j] * stations_powers.second * 24 * 0.4)
            );
        }
    }
    
    model.update();
}

void EVCSP::build_criterion()
{
    GRBLinExpr users_cost = 0.0;
    GRBLinExpr land_rental_cost = 0.0;
    
    GRBLinExpr sum_l2_station_location  = 0.0;
    GRBLinExpr sum_l3_station_location  = 0.0;
    
    for (int i = 0; i < grid_width; ++i)
    {
        for (int j = 0; j < grid_height; ++j)
        {
            sum_l2_station_location += l2_station_location[i][j];
            sum_l3_station_location += l3_station_location[i][j];
            land_rental_cost += land_rental_cost_map[i][j] * station_location[i][j];
        }
    }

    for (int m = 0; m < grid_width; ++m)
    {
        for (int n = 0; n < grid_height; ++n)
        {
            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    if(!std::isinf(distances_from_to[m][n][i][j]))
                        users_cost += demand_map[m][n] * demand_allocation_map[m][n][i][j] * distances_from_to[m][n][i][j];
                }
            }
        }
    }
    
    GRBLinExpr investment_cost = sum_l2_station_location * initial_costs.first + 
                                 sum_l3_station_location * initial_costs.second;

    GRBLinExpr maintenance_cost = sum_l2_station_location * maintenance_costs.first +
                                  sum_l3_station_location * maintenance_costs.second;

    GRBLinExpr operator_cost = investment_cost + land_rental_cost + maintenance_cost;
    GRBLinExpr total_min_cost = operator_cost + users_cost;

    model.setObjective(total_min_cost, GRB_MINIMIZE);
    model.update();
}

void EVCSP::print_solution()
{
    constexpr double INF = std::numeric_limits<double>::infinity();

    for (int n = 0; n < grid_height; ++n)
    {
        for (int m = 0; m < grid_width; ++m)
        {
            if (distances_costs_map[m][n] == INF)
            {
                std::cout << " X ";
            }
            else
            {
                if (station_location[m][n].get(GRB_DoubleAttr_X) > 0.5)
                {
                    int l2_count = std::round(l2_station_location[m][n].get(GRB_DoubleAttr_X));
                    int l3_count = std::round(l3_station_location[m][n].get(GRB_DoubleAttr_X));

                    if (l2_count > 0 && l3_count > 0)
                    {
                        std::cout << l2_count << ":" << l3_count;
                    }
                    else if (l2_count > 0)
                    {
                        std::cout << "2:" << l2_count;
                    }
                    else if (l3_count > 0)
                    {
                        std::cout << "3:" << l3_count;
                    }
                    else
                    {
                        std::cout << " 0 ";
                    }
                }
                else
                {
                    if (poi_map[m][n] == 0)
                    {
                        std::cout << " - ";
                    }
                    else
                    {
                        std::cout << " . ";
                    }
                }
            }
            std::cout << " ";
        }
        std::cout << "\n";
    }
}

void EVCSP::print_demand_distribution()
{
    constexpr double EPSILON = 1e-4; 

    for (int n = 0; n < grid_height; ++n)
    {
        for (int m = 0; m < grid_width; ++m)
        {
            double current_demand = demand_map[m][n];
            
            if (current_demand > 0.0)
            {
                std::cout << "Point (" << m << ", " << n << ") "
                          << "[Full-Demand: " << std::fixed << std::setprecision(1) << current_demand << " kWh]:\n";
                
                double total_allocated_ratio = 0.0;

                for (int i = 0; i < grid_width; ++i)
                {
                    for (int j = 0; j < grid_height; ++j)
                    {
                        double allocation_ratio = demand_allocation_map[m][n][i][j].get(GRB_DoubleAttr_X);

                        if (allocation_ratio > EPSILON)
                        {
                            double allocated_kwh = current_demand * allocation_ratio;
                            
                            std::cout << "  -> " << std::setprecision(1) << (allocation_ratio * 100.0) << "% demand "
                                      << "(" << allocated_kwh << " kWh) "
                                      << "goes to station: (" << i << ", " << j << ") "
                                      << "| distance: " << std::setprecision(2) << distances_from_to[m][n][i][j] << "\n";
                            
                            total_allocated_ratio += allocation_ratio;
                        }
                    }
                }

                if (total_allocated_ratio < (1.0 - EPSILON))
                {
                    std::cout << "  --> BŁĄD: Obsłużono tylko " << (total_allocated_ratio * 100.0) 
                              << "% popytu dla tego punktu!\n";
                }
                std::cout << "----------------------------------------------------------------------\n";
            }
        }
    }
    std::cout << "======================================================================\n\n";
}