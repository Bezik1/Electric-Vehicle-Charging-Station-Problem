#include "EVCSP.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <queue>
#include <limits>

EVCSP::EVCSP(int _grid_width, int _grid_height, int _max_stations_per_cell, double _budget)
    : grid_width{_grid_width}, grid_height{_grid_height}, max_stations_per_cell{_max_stations_per_cell}, budget{_budget}, env{}, model{env}
{
    model.set(GRB_StringAttr_ModelName, MODEL_NAME);
}

EVCSP::EVCSP(int _grid_width, int _grid_height, int _max_stations_per_cell, double _budget, double mip_gap)
    : grid_width{_grid_width}, grid_height{_grid_height}, max_stations_per_cell{_max_stations_per_cell}, budget{_budget}, env{}, model{env}
{
    model.set(GRB_StringAttr_ModelName, MODEL_NAME);
    model.set(GRB_DoubleParam_MIPGap, mip_gap);
}

void EVCSP::operator()(
    std::vector<std::vector<double>>&& _distances_costs_map,
    std::vector<std::vector<int>>&& _poi_map,
    std::vector<std::vector<double>>&& _demand_map,
    std::vector<std::vector<double>>&& _land_rental_cost_map,
    std::pair<double, double>&& _stations_powers,
    std::pair<double, double>&& _initial_costs,
    std::pair<double, double>&& _maintenance_costs)
{
    distances_costs_map = std::move(_distances_costs_map);
    poi_map = std::move(_poi_map);
    demand_map = std::move(_demand_map);
    land_rental_cost_map = std::move(_land_rental_cost_map);
    stations_powers = std::move(_stations_powers);
    initial_costs = std::move(_initial_costs);
    maintenance_costs = std::move(_maintenance_costs);

    demand_allocation_map.assign(grid_width, std::vector<std::vector<std::vector<GRBVar>>>(grid_height, 
        std::vector<std::vector<GRBVar>>(grid_width, std::vector<GRBVar>(grid_height))));

    distances_from_to.assign(grid_width, std::vector<std::vector<std::vector<double>>>(grid_height, 
        std::vector<std::vector<double>>(grid_width, std::vector<double>(grid_height, std::numeric_limits<double>::infinity()))));

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

    build_variables();
    build_constraints();
    build_criterion();

    model.optimize();
}

double EVCSP::calculate_distance(int m, int n, int i, int j) const
{
    using Node = std::pair<double, std::pair<int, int>>;
    constexpr double INF_VAL = -1.0;

    if (m == i && n == j)
        return 0.0;

    if (distances_costs_map[i][j] == INF_VAL)
        return std::numeric_limits<double>::infinity();

    double distanceToNearestValid = 0;
    auto find_nearest_valid = [&](int start_x, int start_y) -> std::pair<int, int>
    {
        if (distances_costs_map[start_x][start_y] != INF_VAL)
            return {start_x, start_y};

        std::queue<std::pair<int, int>> q;
        std::vector<std::vector<bool>> visited(grid_width, std::vector<bool>(grid_height, false));

        q.push({start_x, start_y});
        visited[start_x][start_y] = true;

        while(!q.empty())
        {
            distanceToNearestValid += 1;
            auto [cx, cy] = q.front();
            q.pop();

            for(const auto& [dx, dy] : DIRECTIONS)
            {
                int nx = cx + dx;
                int ny = cy + dy;

                if(nx >= 0 && nx < grid_width && ny >= 0 && ny < grid_height && !visited[nx][ny])
                {
                    if (distances_costs_map[nx][ny] != INF_VAL)
                        return {nx, ny};

                    visited[nx][ny] = true;
                    q.push({nx, ny});
                }
            }
        }
        return {-1, -1};
    };

    std::pair<int, int> actual_start = find_nearest_valid(m, n);
    if (actual_start.first == -1) return std::numeric_limits<double>::infinity();

    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> queue;
    std::vector<std::vector<double>> dist(grid_width, std::vector<double>(grid_height, std::numeric_limits<double>::infinity()));

    dist[actual_start.first][actual_start.second] = distanceToNearestValid;
    queue.push({distanceToNearestValid, actual_start});

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

void EVCSP::build_variables()
{
    station_location.assign(grid_width, std::vector<GRBVar>(grid_height));
    l2_station_location.assign(grid_width, std::vector<GRBVar>(grid_height));
    l3_station_location.assign(grid_width, std::vector<GRBVar>(grid_height));

    for (int i = 0; i < grid_width; ++i)
    {
        for (int j = 0; j < grid_height; ++j)
        {
            station_location[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            l2_station_location[i][j] = model.addVar(0.0, max_stations_per_cell, 0.0, GRB_INTEGER);
            l3_station_location[i][j] = model.addVar(0.0, max_stations_per_cell, 0.0, GRB_INTEGER);
        }
    }

    for (int m = 0; m < grid_width; ++m)
    {
        for (int n = 0; n < grid_height; ++n)
        {
            if (demand_map[m][n] <= 0) continue;
            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    if (!std::isinf(distances_from_to[m][n][i][j])) {
                        demand_allocation_map[m][n][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                    }
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
            if(poi_map[i][j] == 0) {
                model.addConstr(station_location[i][j] == 0);
                model.addConstr(l2_station_location[i][j] == 0);
                model.addConstr(l3_station_location[i][j] == 0);
            }

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
            if (demand_map[m][n] <= 0) continue;

            GRBLinExpr sum = 0.0;
            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    // Obstacle Constraint
                    if(!std::isinf(distances_from_to[m][n][i][j])) {
                        sum += demand_allocation_map[m][n][i][j];
                    }
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
                    if (demand_map[m][n] > 0 && !std::isinf(distances_from_to[m][n][i][j])) {
                        total_demand_at_station += demand_map[m][n] * demand_allocation_map[m][n][i][j];
                    }
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
    GRBLinExpr total_cost = 0.0;

    // Operator Cost: Investment + Maintenance + Land Rental
    for (int i = 0; i < grid_width; ++i)
    {
        for (int j = 0; j < grid_height; ++j)
        {
            total_cost += l2_station_location[i][j] * (initial_costs.first + maintenance_costs.first);
            total_cost += l3_station_location[i][j] * (initial_costs.second + maintenance_costs.second);
            total_cost += land_rental_cost_map[i][j] * station_location[i][j];
        }
    }

    // Users Cost: Travel costs
    for (int m = 0; m < grid_width; ++m)
    {
        for (int n = 0; n < grid_height; ++n)
        {
            if (demand_map[m][n] <= 0) continue;
            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    if(!std::isinf(distances_from_to[m][n][i][j])) {
                        total_cost += demand_map[m][n] * demand_allocation_map[m][n][i][j] * distances_from_to[m][n][i][j];
                    }
                }
            }
        }
    }

    model.setObjective(total_cost, GRB_MINIMIZE);
    model.update();
}

void EVCSP::print_solution()
{
    for (int n = 0; n < grid_width; ++n)
    {
        for (int m = 0; m < grid_height; ++m)
        {
            if (distances_costs_map[n][m] == -1.0)
            {
                std::cout << " X ";
            }
            else
            {
                if (station_location[n][m].get(GRB_DoubleAttr_X) > 0.5)
                {
                    int l2_count = std::round(l2_station_location[n][m].get(GRB_DoubleAttr_X));
                    int l3_count = std::round(l3_station_location[n][m].get(GRB_DoubleAttr_X));

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
                    if (poi_map[n][m] == 0)
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
                            double allocation_ratio = demand_allocation_map[n][m][i][j].get(GRB_DoubleAttr_X);

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