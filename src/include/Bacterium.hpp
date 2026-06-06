#pragma once

#include <vector>
#include <utility>
#include "Maps.hpp"
#include "Evaluator.hpp"

/**
 * @brief Represents types of station maps.
 * 
 */
enum class Map {
    L2_STATION_LOCATION,
    L3_STATION_LOCATION,
};

/**
 * @brief Represents `move` of the bacteria.
 * 
 */
struct Move {
    Map map;
    std::vector<std::vector<int>> move;
};

/**
 * @brief It's the representation of the solution for EVCSP. It was created based on
 * BFO meta heuristic algorithm entities. This class represents single individual solution
 * in the large population of bacterias. They can move forward, change direction and creaate
 * swarms to achieve best solutions.
 * 
 */
class Bacterium
{
public:
    Bacterium(
        const Maps& maps,
        int max_stations_per_cell,
        std::pair<double, double> station_powers
    );
    Bacterium(const Bacterium& other);

    /**
     * @brief Metaphoricly it's the random change of direction of the bacterium.
     * Bacteria cannot move in the direction entirely of their choice, so they
     * inequally move their locomotion cell-part called flagella's to change their
     * direction to explore more nutritient-rich regions.
     * 
     * 
     * In this implementation this function creates new move. Move is just the
     * change in the location of stations.
     *
     * This implementation of tumble is splited into two parts:
     *
     * 1. Repair.
     *
     * If solution have realy high cost model checks if there is an inefficient
     * station currently installed. If there is one program just creates new move
     * of removing the unncessary station and ends the current function call.
     *
     * This mechanism is necessary for performance of the algorithm because, eventually
     * there would be situation, where current bacterium would place station not in the
     * best locations, but the deletion of the station would mean not satisfying the demands,
     * so model just sticks to the incorrect solutions. This fact makes it necessary to implement
     * repair mechanism, that will definietly remove the station and will allow bacterium to develop
     * new station in the different better location.
     *
     * 2. Actual tumble.
     *
     * Model first randomly chooses type of station to try to modify and then
     * it updates cell of the given move map by the value of either -1 / 0 / 1.
     * It makes bacterium either increase, decrease or not modify at all any cell
     * that of course have poi_map value equal to 1.
     * 
     * After choosing mentioned parameters bacterium just sets current move as the
     * choosen one.
     *
     */
    void tumble();


    /**
     * @brief Well-placed bacterium can decide not to change direction and instead
     * swim in the direction that was proved to be lead to nutritient-rich regions.
     * 
     * In this implementation bacterium can also do `swim`, that mean re-create the move
     * that led to the more efficient solution of the model. Swim method just adds current
     * move map into the actual map of the choosen station.
     * 
     * @warning
     * After this operation user must ensure that bacterium is correct, by using `update method`.
     * 
     */
    void swim();

    /**
     * @brief Initialize bacterium, with the randomly placed locations, with the probability
     * of `INITIALIZATION_PROB` in each POI.
     *
     * @warning
     * After this operation user must ensure that bacterium is correct, by using `update method`.
     *
     */
    void initialize();

    /**
     * @brief  Update method recalculates bacterium `current_cost`, recalculates
     * demand_allocation_map and ensures that bacterium meets EVCSP basic constraints.
     * 
     * To achieve that it uses `repair_and_validate_structure` and `recalculate_demand_allocation`
     * functions. They with the cooperation with `Evaluator` allows bacterium to measure optional
     * penalty that will be added to the `current_cost` if needed.
     * 
     * @param evaluator Calculates cost of the current solution stored in Bacterium.
     */
    void update(const Evaluator& evaluator);

    const Maps& maps;
    int max_stations_per_cell;
    std::pair<double, double> station_powers;

    Move current_move;
    
    std::vector<std::vector<int>> station_location;
    std::vector<std::vector<int>> l2_station_location;
    std::vector<std::vector<int>> l3_station_location;
    std::vector<std::vector<std::vector<std::vector<double>>>> demand_allocation_map;

    double current_cost = 0.0;
    double health = 0.0;

private:
    /**
     * @brief Ensures the validity of the Bacterium and if neeeded repairs it
     * and returns the cost of the repairs.
     * 
     * It checks following constraints:
     * 
     * 1. That at POI = 0, there should be no station.
     * 
     * 2. It cut's the station_location value range to {0, 1}.
     * 
     * 3. It checks, if the number of stations at the given cell match
     * the number of `max_stations_per_cell`. In case it would not be 
     * true function tries to subtruct stations from the given cell. First
     * from L3 stations, then from L2 stations.
     * 
     * @return double Penalty for repairing the Bacterium.
     */
    double repair_and_validate_structure();

    /**
     * @brief Bacterium recalculates values of demand_allocation_map in order
     * to make them as efficient as they can be.
     * 
     * This method is necessary in our case because Bacterium do not engage
     * demand_allocaiton_map in the `tumble` process, so it should greedy
     * check the potential stations that will handle user demands.
     * 
     * This method tries to allocate demand of the user to the nearest station
     * and if it is overloaded, then it goes to the next station, based on distance.
     * If, there are no stations that can handle the rest of the demand, then program
     * breaks the loop and Bacterium has unmet demand. This situation is handled, by the
     * evaluator that account big penalty for unmet demands. In that way model can learn
     * that it should create, as many stations as the user demand require. 
     * 
     */
    void recalculate_demand_allocation();

    static constexpr double INITIALIZATION_PROB = 0.3;
    static constexpr double CHANGE_PROB = 0.25; 

    static constexpr double HARD_THRESHOLD = 1e6;
    static constexpr double MISC_PENALTY = 1e7;
    static constexpr double EPS = 1e-9;
};