#include "SolutionPrinter.hpp"

SolutionPrinter::SolutionPrinter(const Maps& maps): maps{maps} {}

void SolutionPrinter::print_iteration(
    const Solution& solution,
    bool validity,
    int elimination_step,
    int reproduction_step,
    int hemotaxis_step
)
{
    std::println("E: {:2}, R: {:2}, H: {:2} | Cost: {:12.2f} | Validity: {:d}",
        elimination_step,
        reproduction_step,
        hemotaxis_step,
        solution.total_cost,
        validity ? 1 : 0
    );
}

void SolutionPrinter::print_map(const Solution& solution)
{
    std::println("Best Cost: {:.2f}", solution.total_cost);

    for (int n = 0; n < maps.get_width(); ++n)
    {
        for (int m = 0; m < maps.get_height(); ++m)
        {
            std::string cell = "";

            if (!maps.isPassable(n, m))
            {
                cell = "X";
            }
            else if (solution.station_location[n][m] > 0)
            {
                int l2 = solution.l2_station_location[n][m];
                int l3 = solution.l3_station_location[n][m];

                if (l2 > 0 && l3 > 0)
                    cell = std::format("L2:{}L3:{}", l2, l3);
                else if (l2 > 0)
                    cell = std::format("L2:{}", l2);
                else if (l3 > 0)
                    cell = std::format("L3:{}", l3);
                else
                    cell = " ";
            }
            else
            {
                cell = (maps.get_poi_at(n, m) == 0) ? " " : ".";
            }

            std::print("{:^9}", cell);
        }
        std::println("");
    }
}

void SolutionPrinter::print_demand_distribution(const Solution& solution)
{
    for (int n = 0; n < maps.get_width(); ++n)
    {
        for (int m = 0; m < maps.get_height(); ++m)
        {
            double current_demand = maps.get_demand_at(n, m);
            
            if (current_demand > 0.0)
            {
                std::println("Point ({}, {}), [Full-Demand: {:.2f} kWh]:", n, m, current_demand);

                for (int i = 0; i < maps.get_width(); ++i)
                {
                    for (int j = 0; j < maps.get_height(); ++j)
                    {
                        double allocation_ratio = solution.demand_allocation_map[n][m][i][j];

                        if (allocation_ratio <= EPSILON)
                            continue;

                        double allocated_kwh = current_demand * allocation_ratio;
                        std::println(
                            " -> {:.1f}% demand ({:.2f} kWh) goes to station: ({}, {}) | distance: {:.2f}.",
                            allocation_ratio * 100.0,
                            allocated_kwh,
                            i,
                            j,
                            maps.get_distance(n, m, i, j) 
                        );
                    }
                }
                std::println("----------------------------------------------------------------------");
            }
        }
    }
    std::println("======================================================================\n");
}