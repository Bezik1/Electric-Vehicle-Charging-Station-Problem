#pragma once

#include <array>
#include <vector>
#include <memory>
#include <limits>
#include <variant>
#include <mutex>
#include <thread>

#include "Maps.hpp"
#include "Model.hpp"
#include "Solution.hpp"
#include "Evaluator.hpp"
#include "Validator.hpp"
#include "Bacterium.hpp"
#include "SolutionPrinter.hpp"

/**
 * @brief Copy of the l2 and l3 station locations. It's purpose is to hold
 * information about station locations of each bacterium between each hemotaxis
 * step to pass this information for calculation of swarming effect.
 * 
 * This structure is necessary for the correct working mechanism of multi-threading
 * in this application.
 * 
 */
struct BacteriumSnapshot {
    std::vector<std::vector<int>> l2_station_location;
    std::vector<std::vector<int>> l3_station_location;
};

/**
 * @brief MetaheuristicModel implements BFO meta heuristic algorithm. It is used
 * to solve EVCSP. It creates and manages the lifecycle of bacterias in order to generate
 * optimal solution. BFOA concentrates on the movement of bacterium towards the best solution.
 * To do that bacterias are using tumble, swim and reproduction mechanism. MetaheuristicModel
 * additionally resets and eliminates some part of the population to avoid stagnation.
 *
 */
class MetaheuristicModel : public Model
{
public:
    explicit MetaheuristicModel(
        int max_stations_per_cell,
        int bacteria_count,
        int hemotaxis_steps,
        int swimming_steps,
        int reproduction_steps,
        int elimination_dispersal_steps,
        double elimination_dispersal_prob,
        double d_attract,
        double w_attract,
        double h_repellant,
        double w_repellant,
        double budget,
        int num_threads = DEFAULT_NUM_THREADS
    );

    void operator()(
        const Maps& input_maps,
        std::pair<double, double> _stations_powers,
        std::pair<double, double> _initial_costs,
        std::pair<double, double> _maintenance_costs
    ) override;

    /**
     * @brief Main method of the model. It manages the whole meta heuristic
     * optimization process, by creating new population and launching 4-layered loop.
     * 
     * Firstly model creates each Bacterium and calls `initialize` method to set their
     * initial solutions as the random ones.
     * 
     * Mentioned 4-layered loop contains following phases (from the most inner-loop 
     * to the most-outer one):
     * 
     * 1. Swimming Loop.
     * 
     * In this loop bacterias swim in the direction,of the last tumble
     * to the point, when the swimming direction is no longer optimal, or
     * the number of `swimming_steps` is exceeded. They use `swim` method
     * to achieve that.
     * 
     * 2. Hemotaxis Loop.
     * 
     * In this loop bacterias change their direction using `tumble` method.
     * They save their current_move and they swim one time in that direction,
     * after that Bacterium checks if the direction is the beneficial one and if
     * it is, when the swimming loop is initiated.
     * 
     * 3. Reproduction Loop.
     * 
     * In this loop program sorts the solutions based on health of each Bacterium.
     * The healthier half of the population reproduce themselves using cell divison
     * creating identical copies of themselves. The other half is replaced by the 
     * new clones of the healthier half of population and dies.
     * 
     * 4. Dispersal Loop.
     * Random number of the individuals from the population is being replaced by the 
     * completly random ones to avoid stagnation.
     * 
     *
     */
    void optimization_loop();
    

    Solution get_solution() const override;

    /**
     * @brief It's the pointer to the object that can print
     * solutions ot the model.
     *
     */
    std::unique_ptr<SolutionPrinter> printer = nullptr;

private:
    /**
     * @brief Calcualtes swarming effect, that affect other bacterias to either
     * attract them to the places, were solutions are great. Other use of this
     * mechanism is to repell bacterias from the overcrowded areas to avoid
     * collecting all population in some local minimum.
     * 
     * Formula for the swarming effect:
     * 
     * swarming_cost += -d_attract * exp(-w_attract * distance_sq);
     *   swarming_cost += h_repellant * exp(-w_repellant * distance_sq);
     * 
     * 
     * @param current_bacterium Bacterium that should be affected by the swarming effect.
     * @param snapshot Copy of the l2 and l3 station locations.
     * 
     * @return double 
     */
    double calculate_swarming_effect(
        const Bacterium& current_bacterium,
        const std::vector<BacteriumSnapshot>& snapshot
    ) const;

    
    /**
     * @brief Represents the entire population of solutions as the bacterias.
     * 
     */
    std::vector<std::unique_ptr<Bacterium>> colony;
    
    /**
     * @brief Pointer to the Bacterium that solution has the lowest cost.
     * 
     */
    std::unique_ptr<Bacterium> best_bacterium = nullptr;

    /**
     * @brief Value of the cost of `best_bacterium`.
     * 
     */
    double best_cost = std::numeric_limits<double>::infinity();

    const Maps* maps = nullptr;
    std::unique_ptr<Validator> validator = nullptr;
    std::unique_ptr<Evaluator> evaluator = nullptr;

    std::pair<double, double> stations_powers;
    std::pair<double, double> initial_costs;
    std::pair<double, double> maintenance_costs;

    double d_attract;
    double w_attract;
    double h_repellant;
    double w_repellant;

    int max_stations_per_cell;
    double budget;
    int bacteria_count; 
    int hemotaxis_steps;
    int swimming_steps;
    int reproduction_steps;
    int elimination_dispersal_steps;
    double elimination_dispersal_prob;

    int num_threads;
    std::mutex best_solution_mutex;

    static constexpr int DEFAULT_NUM_THREADS = 10;
};