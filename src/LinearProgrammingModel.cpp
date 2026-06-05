#include "LinearProgrammingModel.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <queue>
#include <limits>

LinearProgrammingModel::LinearProgrammingModel(int grid_width, int grid_height, int max_stations_per_cell, double budget)
    : grid_width{grid_width}, grid_height{grid_height}, max_stations_per_cell{max_stations_per_cell}, budget{budget}, env{}, model{env}, map_data{nullptr}
{
    model.set(GRB_StringAttr_ModelName, MODEL_NAME);
    model.set(GRB_DoubleParam_MIPGap, DEFAULT_MIP_GAP);
}

LinearProgrammingModel::LinearProgrammingModel(int grid_width, int grid_height, int max_stations_per_cell, double budget, double mip_gap)
    : grid_width{grid_width}, grid_height{grid_height}, max_stations_per_cell{max_stations_per_cell}, budget{budget}, env{}, model{env}, map_data{nullptr}
{
    model.set(GRB_StringAttr_ModelName, MODEL_NAME);
    model.set(GRB_DoubleParam_MIPGap, mip_gap);
}

void LinearProgrammingModel::operator()(
    const Maps& _map_data,
    std::pair<double, double> _stations_powers,
    std::pair<double, double> _initial_costs,
    std::pair<double, double> _maintenance_costs
) 
{
    map_data = &_map_data;
    stations_powers = std::move(_stations_powers);
    initial_costs = std::move(_initial_costs);
    maintenance_costs = std::move(_maintenance_costs);

    demand_allocation_map.assign(
        grid_width, std::vector<std::vector<std::vector<GRBVar>>>(
            grid_height, 
            std::vector<std::vector<GRBVar>>(grid_width, std::vector<GRBVar>(grid_height))
    ));

    build_variables();
    build_constraints();
    build_criterion();

    model.optimize();
}

void LinearProgrammingModel::build_variables()
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
            if (map_data->get_demand_at(m, n) <= 0) continue;
            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    if (!std::isinf(map_data->get_distance(m, n, i, j))) {
                        demand_allocation_map[m][n][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                    }
                }
            }
        }
    }
    model.update();
}

void LinearProgrammingModel::build_constraints()
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
            if(map_data->get_poi_at(i, j) == 0) {
                model.addConstr(station_location[i][j] == 0);
                model.addConstr(l2_station_location[i][j] == 0);
                model.addConstr(l3_station_location[i][j] == 0);
            }

            total_investment_cost += (initial_costs.first * l2_station_location[i][j] + 
                                      initial_costs.second * l3_station_location[i][j] + 
                                      map_data->get_rental_cost_at(i, j) * station_location[i][j]);
        }
    }
    
    // Budget Constraint
    model.addConstr(total_investment_cost <= budget);

    for (int m = 0; m < grid_width; ++m)
    {
        for (int n = 0; n < grid_height; ++n)
        {
            double d_val = map_data->get_demand_at(m, n);
            if (d_val <= 0) continue;

            GRBLinExpr sum = 0.0;
            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    // Obstacle Constraint
                    if(!std::isinf(map_data->get_distance(m, n, i, j))) {
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
                    double d_val = map_data->get_demand_at(m, n);
                    if (d_val > 0 && !std::isinf(map_data->get_distance(m, n, i, j))) {
                        total_demand_at_station += d_val * demand_allocation_map[m][n][i][j];
                    }
                }
            }

            // Capacity Constraint
            model.addConstr(
                total_demand_at_station <=
                    (l2_station_location[i][j] * stations_powers.first +
                    l3_station_location[i][j] * stations_powers.second) * 24
            );
        }
    }
    
    model.update();
}

void LinearProgrammingModel::build_criterion()
{
    GRBLinExpr total_cost = 0.0;

    // Operator Cost: Investment + Maintenance + Land Rental
    for (int i = 0; i < grid_width; ++i)
    {
        for (int j = 0; j < grid_height; ++j)
        {
            total_cost += l2_station_location[i][j] * (initial_costs.first + maintenance_costs.first);
            total_cost += l3_station_location[i][j] * (initial_costs.second + maintenance_costs.second);
            total_cost += map_data->get_rental_cost_at(i, j) * station_location[i][j];
        }
    }

    // Users Cost: Travel costs
    for (int m = 0; m < grid_width; ++m)
    {
        for (int n = 0; n < grid_height; ++n)
        {
            double d_val = map_data->get_demand_at(m, n);
            if (d_val <= 0) continue;
            for (int i = 0; i < grid_width; ++i)
            {
                for (int j = 0; j < grid_height; ++j)
                {
                    double dist = map_data->get_distance(m, n, i, j);
                    if(!std::isinf(dist)) {
                        total_cost += d_val * demand_allocation_map[m][n][i][j] * dist;
                    }
                }
            }
        }
    }

    model.setObjective(total_cost, GRB_MINIMIZE);
    model.update();
}

LinearProgrammingModel::Solution LinearProgrammingModel::get_solution() const {
    LinearProgrammingModel::Solution sol;
    
    sol.station_location.assign(grid_width, std::vector<int>(grid_height));
    sol.l2_station_location.assign(grid_width, std::vector<int>(grid_height));
    sol.l3_station_location.assign(grid_width, std::vector<int>(grid_height));
    sol.demand_allocation_map.assign(grid_width, 
        std::vector<std::vector<std::vector<double>>>(grid_height, 
            std::vector<std::vector<double>>(grid_width, 
                std::vector<double>(grid_height, 0.0))));

    try {
        int optim_status = model.get(GRB_IntAttr_Status);

        if (optim_status == GRB_OPTIMAL || optim_status == GRB_SUBOPTIMAL) {
            sol.total_cost = model.get(GRB_DoubleAttr_ObjVal);

            for (int i = 0; i < grid_width; ++i) {
                for (int j = 0; j < grid_height; ++j) {
                    sol.station_location[i][j] = (int)(station_location[i][j].get(GRB_DoubleAttr_X) + 0.5);
                    sol.l2_station_location[i][j] = (int)(l2_station_location[i][j].get(GRB_DoubleAttr_X) + 0.5);
                    sol.l3_station_location[i][j] = (int)(l3_station_location[i][j].get(GRB_DoubleAttr_X) + 0.5);

                    for (int m = 0; m < grid_width; ++m) {
                        for (int n = 0; n < grid_height; ++n) {
                            if (map_data->get_demand_at(m, n) > 0 && !std::isinf(map_data->get_distance(m, n, i, j))) {
                                sol.demand_allocation_map[m][n][i][j] = demand_allocation_map[m][n][i][j].get(GRB_DoubleAttr_X);
                            }
                        }
                    }
                }
            }
        } else {
            std::cerr << "Model has not been optimized or no solution found!" << std::endl;
            sol.total_cost = std::numeric_limits<double>::infinity();
        }
    } catch (GRBException& e) {
        std::cerr << "Error during getting solution: " << e.getMessage() << std::endl;
    }

    return sol;
}

void LinearProgrammingModel::print_solution()
{
    try {
        int optim_status = model.get(GRB_IntAttr_Status);
        if (optim_status != GRB_OPTIMAL && optim_status != GRB_SUBOPTIMAL) return;

        for (int n = 0; n < grid_width; ++n)
        {
            for (int m = 0; m < grid_height; ++m)
            {
                if (!map_data->isPassable(n, m))
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
                        if (map_data->get_poi_at(n, m) == 0)
                            std::cout << " - ";
                        else
                            std::cout << " . ";
                    }
                }
                std::cout << " ";
            }
            std::cout << "\n";
        }
    } catch (...) {}
}

void LinearProgrammingModel::print_demand_distribution()
{
    constexpr double EPSILON = 1e-4; 
    
    try {
        int optim_status = model.get(GRB_IntAttr_Status);
        if (optim_status != GRB_OPTIMAL && optim_status != GRB_SUBOPTIMAL) return;

        for (int n = 0; n < grid_width; ++n)
        {
            for (int m = 0; m < grid_height; ++m)
            {
                double current_demand = map_data->get_demand_at(n, m);
                
                if (current_demand > 0.0)
                {
                    std::cout << "Point (" << m << ", " << n << ") "
                              << "[Full-Demand: " << std::fixed << std::setprecision(1) << current_demand << " kWh]:\n";
                    
                    for (int i = 0; i < grid_width; ++i)
                    {
                        for (int j = 0; j < grid_height; ++j)
                        {
                            double dist = map_data->get_distance(n, m, i, j);
                            if (!std::isinf(dist))
                            {
                                double allocation_ratio = demand_allocation_map[n][m][i][j].get(GRB_DoubleAttr_X);

                                if (allocation_ratio > EPSILON)
                                {
                                    double allocated_kwh = current_demand * allocation_ratio;
                                    
                                    std::cout << "  -> " << std::setprecision(1) << (allocation_ratio * 100.0) << "% demand "
                                              << "(" << allocated_kwh << " kWh) "
                                              << "goes to station: (" << i << ", " << j << ") "
                                              << "| distance: " << std::setprecision(2) << dist << "\n";
                                }
                            }
                        }
                    }
                    std::cout << "----------------------------------------------------------------------\n";
                }
            }
        }
    } catch (...) {}
    std::cout << "======================================================================\n\n";
}