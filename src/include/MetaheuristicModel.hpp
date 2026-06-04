#pragma once

#include <array>
#include <vector>
#include <memory>
#include <limits>
#include <variant>

#include "Model.hpp"
#include "EvaluationContext.hpp"
#include "Evaluator.hpp"

class MetaheuristicModel : public Model, public EvaluationContext
{
public:
    MetaheuristicModel(
        int grid_width,
        int grid_height,
        int max_stations_per_cell,
        int bacteria_count,
        int hemotaxis_steps,
        int swimming_steps,
        int reproduction_steps,
        int elimination_dispersal_steps,
        double elimination_dispersal_prob,
        double budget
    );

    void operator()(
        std::vector<std::vector<double>> _distances_costs_map,
        std::vector<std::vector<int>> _poi_map,
        std::vector<std::vector<double>> _demand_map,
        std::vector<std::vector<double>> _land_rental_cost_map,
        std::pair<double, double> _stations_powers,
        std::pair<double, double> _initial_costs,
        std::pair<double, double> _maintenance_costs
    ) override;

    void optimization_loop();
    
    void print_solution() override;
    
    void print_demand_distribution() override;

    int getWidth() const override { return grid_width; }
    int getHeight() const override { return grid_height; }

    int getPoiAt(int i, int j) const override { return poi_map[i][j]; }
    double getDemandAt(int i, int j) const override { return demand_map[i][j]; }
    double getRentalCostAt(int i, int j) const override { return land_rental_cost_map[i][j]; }
    double getDistance(int m, int n, int i, int j) const override { return distances_from_to[m][n][i][j]; }

    int getMaxStationsPerCell() const override { return max_stations_per_cell; }
    double getBudget() const override { return budget; }

    std::pair<double, double> getStationPowers() const override { return stations_powers; }
    std::pair<double, double> getInitialCosts() const override { return initial_costs; }
    std::pair<double, double> getMaintenanceCosts() const override { return maintenance_costs; }

    struct Bacterium
    {
        enum class Map {
            STATION_LOCATION,
            L2_STATION_LOCATION,
            L3_STATION_LOCATION,
        };
        
        struct Move {
            Bacterium::Map map;
            std::vector<std::vector<int>>  data;
        };
        
        Bacterium(int grid_width, int grid_height);
        Bacterium(const Bacterium& other);
        
        void tumble(const MetaheuristicModel& model);
        void swim(const MetaheuristicModel& model);
        void repair_and_validate_structure(const MetaheuristicModel& model, double& penalty);
        void recalculate_demand_allocation(const MetaheuristicModel& model);
        
        double current_cost = 0.0;
        double health = 0.0;
        Move current_move;
        
        std::vector<std::vector<int>> station_location;
        std::vector<std::vector<int>> l2_station_location;
        std::vector<std::vector<int>> l3_station_location;
        std::vector<std::vector<std::vector<std::vector<double>>>> demand_allocation_map;
    };
    
    std::unique_ptr<Bacterium> best_bacterium = nullptr;
    double best_cost = std::numeric_limits<double>::infinity();

private:
    double calculate_swarming_effect(const Bacterium &current_bacterium) const;
    double calculate_distance(int m, int n, int i, int j) const;
    std::pair<double, std::pair<int, int>> find_nearest_valid(int start_x, int start_y) const;
    bool is_valid(const std::unique_ptr<Bacterium>& bacterium) const;

    std::vector<std::vector<double>> distances_costs_map;
    std::vector<std::vector<int>> poi_map;
    std::vector<std::vector<double>> demand_map;
    std::vector<std::vector<double>> land_rental_cost_map;

    std::pair<double, double> stations_powers;
    std::pair<double, double> initial_costs;
    std::pair<double, double> maintenance_costs;

    std::vector<std::vector<std::vector<std::vector<double>>>> distances_from_to;
    std::vector<std::unique_ptr<Bacterium>> colony;
    std::unique_ptr<Evaluator> evaluator = nullptr;

    int grid_width;
    int grid_height;
    int max_stations_per_cell;
    double budget;
    int bacteria_count; 
    int hemotaxis_steps;
    int swimming_steps;
    int reproduction_steps;
    int elimination_dispersal_steps;
    double elimination_dispersal_prob;
    
    static constexpr double INF_VAL = -1.0;
    static constexpr std::array<std::pair<int, int>, 4> DIRECTIONS = {{ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }};
};