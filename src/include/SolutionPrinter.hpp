#pragma once

#include "Maps.hpp"
#include "Solution.hpp"

#include <print>

class SolutionPrinter {
public:
    explicit SolutionPrinter(const Maps &maps);

    void print_iteration(
        const Solution& solution,
        bool validity,
        int elimination_step,
        int reproduction_step,
        int hemotaxis_step
    );

    void print_map(const Solution& solution);

    void print_demand_distribution(const Solution& solution);
    
private:
    const Maps &maps;

    static constexpr double EPSILON = 1e-9;
};