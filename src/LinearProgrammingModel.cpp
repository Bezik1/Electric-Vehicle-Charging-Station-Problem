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

    printer = std::make_unique<SolutionPrinter>(*map_data);

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

Solution LinearProgrammingModel::get_solution() const {
    Solution solution;
    
    solution.station_location.assign(grid_width, std::vector<int>(grid_height));
    solution.l2_station_location.assign(grid_width, std::vector<int>(grid_height));
    solution.l3_station_location.assign(grid_width, std::vector<int>(grid_height));
    solution.demand_allocation_map.assign(grid_width, 
        std::vector<std::vector<std::vector<double>>>(grid_height, 
            std::vector<std::vector<double>>(grid_width, 
                std::vector<double>(grid_height, 0.0))));

    solution.total_cost = model.get(GRB_DoubleAttr_ObjVal);

    for (int i = 0; i < grid_width; ++i) {
        for (int j = 0; j < grid_height; ++j) {
            solution.station_location[i][j] = (int)(station_location[i][j].get(GRB_DoubleAttr_X) + 0.5);
            solution.l2_station_location[i][j] = (int)(l2_station_location[i][j].get(GRB_DoubleAttr_X) + 0.5);
            solution.l3_station_location[i][j] = (int)(l3_station_location[i][j].get(GRB_DoubleAttr_X) + 0.5);

            for (int m = 0; m < grid_width; ++m) {
                for (int n = 0; n < grid_height; ++n) {
                    if (map_data->get_demand_at(m, n) > 0 && !std::isinf(map_data->get_distance(m, n, i, j))) {
                        solution.demand_allocation_map[m][n][i][j] = demand_allocation_map[m][n][i][j].get(GRB_DoubleAttr_X);
                    }
                }
            }
        }
    }

    return solution;
}