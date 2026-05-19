#include "EVCSP.hpp"

EVCSP::EVCSP(
    double _mean_station_energy_demand,
    int _num_stations,
    int _map_width,
    int _map_height
):  env{},
    model{env},
    num_stations{_num_stations},
    map_width{_map_width},
    map_height{_map_height},
    mean_station_energy_demand{_mean_station_energy_demand}
{
    model.set(GRB_StringAttr_ModelName, MODEL_NAME);
}

void EVCSP::build_variables()
{
    grid_vars.resize(map_width, std::vector<GRBVar>(map_height));

    for (int x = 0; x < map_width; ++x)
    {
        for (int y = 0; y < map_height; ++y)
        {
            grid_vars[x][y] = model.addVar(
                0.0,
                1.0,
                0.0,
                GRB_BINARY
            );
        }
    }
    model.update();
}

void EVCSP::build_constraints()
{
    GRBLinExpr total_stations = 0;

    for (int x = 0; x < map_width; ++x)
    {
        for (int y = 0; y < map_height; ++y)
        {
            total_stations += grid_vars[x][y];

            if (connectivity_map[x][y] == 1)
                model.addConstr(grid_vars[x][y] == 0);
        }
    }

    model.addConstr(total_stations == num_stations);
    model.update();
}

void EVCSP::build_criterion()
{
    GRBLinExpr objective = 0;

    for (int x = 0; x < map_width; ++x)
    {
        for (int y = 0; y < map_height; ++y)
        {
            double commercial_potential = attractiveness_map[x][y] * activity_map[x][y];
            double energy_risk = mean_station_energy_demand / 
                                (energy_capacity_map[x][y] + MINIMAL_ENERGY_DEMAND);
            double score = (COMMERCIAL_POTENTIAL_COEFF * commercial_potential) 
                            - (NETWORK_SECURITY_COEFF * energy_risk);

            objective += score * grid_vars[x][y];
        }
    }

    model.setObjective(objective, GRB_MAXIMIZE);
    model.update();
}

void EVCSP::print_map()
{
    std::printf("\n=========================================\n");
    
    int status = model.get(GRB_IntAttr_Status);
    if (status != GRB_OPTIMAL) 
    {
        std::printf(" No optimal solution available to print. Status: %d\n", status);
        std::printf("=========================================\n");
        return;
    }

    std::printf(" Optimization Success!\n");
    std::printf(" Objective Value: %.2f\n", model.get(GRB_DoubleAttr_ObjVal));
    std::printf("=========================================\n\n");

    std::printf("Selected Charging Station Positions (X, Y):\n");
    for (int x = 0; x < map_width; ++x) 
    {
        for (int y = 0; y < map_height; ++y) 
        {
            if (grid_vars[x][y].get(GRB_DoubleAttr_X) > 0.5) 
            {
                std::printf(" -> Station located at: [%d, %d]\n", x, y);
            }
        }
    }

    std::printf("\nMap Visualization [ . Empty | X Obstacle | S Station ]:\n");
    for (int y = 0; y < map_height; ++y) 
    {
        for (int x = 0; x < map_width; ++x) 
        {
            if (connectivity_map[x][y] == 1) 
                std::printf(" X ");
            else if (grid_vars[x][y].get(GRB_DoubleAttr_X) > 0.5) 
                std::printf(" S ");
            else 
                std::printf(" . ");
        }
        std::printf("\n");
    }
    std::printf("=========================================\n");
}

void EVCSP::operator()(
            std::vector<std::vector<int>>&& _connectivity_map,
            std::vector<std::vector<double>>&& _activity_map,
            std::vector<std::vector<double>>&& _attractiveness_map,
            std::vector<std::vector<double>>&& _energy_capacity_map
        )
{
    connectivity_map = std::move(_connectivity_map);
    activity_map = std::move(_activity_map);
    attractiveness_map = std::move(_attractiveness_map);
    energy_capacity_map = std::move(_energy_capacity_map);

    build_variables();
    build_constraints();
    build_criterion();

    model.optimize();

    print_map();
}